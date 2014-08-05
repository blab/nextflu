# Run when starting up app

import os
import rethinkdb as r
from rethinkdb.errors import RqlRuntimeError, RqlDriverError

RDB_HOST =  os.environ.get('RDB_HOST') or 'localhost'
RDB_PORT = os.environ.get('RDB_PORT') or 28015

def db_setup():
    connection = r.connect(host=RDB_HOST, port=RDB_PORT)
    try:
        r.db_create('augur').run(connection)
        print "augur db created"
    except RqlRuntimeError:
        print "augur db exists"
    try:
        r.db('augur').table_create('seq').run(connection) 
        print "seq table created"
    except RqlRuntimeError:
        print "seq table exists"               
    finally:
        connection.close()

def main():
	db_setup()

if __name__ == "__main__":
    main()
