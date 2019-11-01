#!/usr/bin/env python
import argparse
from livereload import Server, shell

parser = argparse.ArgumentParser(description='Watch and serve documentation.')
parser.add_argument('--host,h', dest='host', type=str, default=None)
parser.add_argument('--port,p', dest='port', type=str, default=5500)

args = parser.parse_args()

server = Server()
server.watch('source/*.rst', shell('make html', cwd='.'))
server.serve(port=args.port, host=args.host, root='build/html')

