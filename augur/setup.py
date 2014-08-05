# Run when starting up app

import os
import rethinkdb as r
from rethinkdb.errors import RqlRuntimeError, RqlDriverError

RDB_HOST =  os.environ.get('RDB_HOST') or 'localhost'
RDB_PORT = os.environ.get('RDB_PORT') or 28015
RDB_DB = 'augur'

def db_setup():

	try:
		connection = r.connect(host=RDB_HOST, port=RDB_PORT)
	except RqlDriverError:
		abort(503, "No database connection could be established.")

	try:
		r.db_create(RDB_DB).run(connection)
	except RqlRuntimeError:
		pass

	try:
		r.db(RDB_DB).table_create('seq').run(connection) 
	except RqlRuntimeError:
		pass               
	
	connection.close()

def main():
	db_setup()

if __name__ == "__main__":
	main()
