'use strict';

/*
 * STRATA — Hardened Local Static Server (zero external dependencies)
 * Binds strictly to 127.0.0.1, enforces canonical root containment,
 * allow-lists GET/HEAD + static extensions, and applies strict security headers.
 */
const http = require('node:http');
const fs = require('node:fs');
const path = require('node:path');

const HOST = '127.0.0.1';
const PORT = Number(process.env.STRATA_PORT || 4174);
const ROOT_DIR = path.resolve(__dirname);

const MIME_ALLOWLIST = new Map([
  ['.html', 'text/html; charset=utf-8'],
  ['.css', 'text/css; charset=utf-8'],
  ['.js', 'application/javascript; charset=utf-8'],
  ['.png', 'image/png'],
  ['.mp4', 'video/mp4'],
  ['.svg', 'image/svg+xml'],
  ['.json', 'application/json; charset=utf-8']
]);

const SECURITY_HEADERS = {
  'Content-Security-Policy':
    "default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data: blob:; media-src 'self'; connect-src 'none'; object-src 'none'; base-uri 'none'; form-action 'none'; frame-ancestors 'none';",
  'X-Content-Type-Options': 'nosniff',
  'X-Frame-Options': 'DENY',
  'Referrer-Policy': 'no-referrer',
  'Permissions-Policy': 'camera=(), microphone=(), geolocation=(), payment=(), usb=()',
  'Cache-Control': 'no-store'
};

function sendError(res, status, message) {
  res.writeHead(status, Object.assign({ 'Content-Type': 'text/plain; charset=utf-8' }, SECURITY_HEADERS));
  res.end(message);
}

const server = http.createServer((req, res) => {
  if (req.method !== 'GET' && req.method !== 'HEAD') {
    res.setHeader('Allow', 'GET, HEAD');
    sendError(res, 405, 'Method Not Allowed');
    return;
  }

  let pathname;
  try {
    const parsed = new URL(req.url || '/', `http://${HOST}:${PORT}`);
    pathname = decodeURIComponent(parsed.pathname);
  } catch {
    sendError(res, 400, 'Bad Request');
    return;
  }

  if (pathname === '/') pathname = '/index.html';
  if (pathname.includes('\0')) {
    sendError(res, 400, 'Bad Request');
    return;
  }

  const resolved = path.resolve(ROOT_DIR, '.' + pathname);
  if (!resolved.startsWith(ROOT_DIR + path.sep)) {
    sendError(res, 403, 'Forbidden');
    return;
  }

  const rel = path.relative(ROOT_DIR, resolved);
  const segments = rel.split(path.sep);
  if (segments.some((s) => s.startsWith('.')) || segments[0] === 'reference' || rel === 'server.js') {
    sendError(res, 403, 'Forbidden');
    return;
  }

  const ext = path.extname(resolved).toLowerCase();
  const contentType = MIME_ALLOWLIST.get(ext);
  if (!contentType) {
    sendError(res, 403, 'Forbidden');
    return;
  }

  fs.stat(resolved, (statErr, stats) => {
    if (statErr || !stats.isFile()) {
      sendError(res, 404, 'Not Found');
      return;
    }

    res.writeHead(
      200,
      Object.assign({ 'Content-Type': contentType, 'Content-Length': stats.size }, SECURITY_HEADERS)
    );
    if (req.method === 'HEAD') {
      res.end();
      return;
    }
    const stream = fs.createReadStream(resolved);
    stream.on('error', () => {
      if (!res.headersSent) sendError(res, 500, 'Internal Server Error');
      else res.destroy();
    });
    stream.pipe(res);
  });
});

// TODO(security): Local single-user static visualization server over loopback (127.0.0.1);
// authentication, rate limiting, and TLS are intentionally omitted for local offline use.
server.listen(PORT, HOST, () => {
  console.log(`STRATA Geomechanics Lab running at http://${HOST}:${PORT}`);
});
