#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description='Watch and serve documentation.')
parser.add_argument('--host,h', dest='host', type=str, default=None)
parser.add_argument('--port,p', dest='port', type=str, default=5500)
parser.add_argument('--watch', dest='watch', default=False, action='store_true')

args = parser.parse_args()

if args.watch:
    from livereload import Server, shell
    server = Server()

    server.watch('source/*.rst', shell('make html', cwd='.'))

    server.serve(port=args.port, host=args.host, root='build/html')
else:
    import http.server
    import socketserver

    host = args.host if args.host is not None else "127.0.0.1"
    port = args.port

    Handler = lambda request, client_address, server: http.server.SimpleHTTPRequestHandler(
                request, client_address, server, directory='build/html')

    with socketserver.TCPServer((host, port), Handler) as httpd:
        print ("Serving on http://%s:%d" % (host, port))
        httpd.serve_forever()

