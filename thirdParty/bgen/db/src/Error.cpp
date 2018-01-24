
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include "db/SQLite3Statement.hpp"
#include "db/Error.hpp"

namespace db {
	Error::Error( std::string const& caller, std::string const& db_spec, int error, std::string const& sql ):
		m_spec( db_spec ),
		m_error( error ),
		m_sql( sql )
	{
		assert( m_error != SQLite3Statement::Error::OK ) ;
	}
	
	Error::~Error() throw() {}
	
	std::string Error::description() const {
		typedef db::SQLite3Statement::Error Error ;
		std::string result ;
		switch( m_error ) {
			case Error::ERROR: 		result = "SQL error or missing database"; break ;
			case Error::INTERNAL: 	result = "Internal logic error in SQLite"; break ;
			case Error::PERM: 		result = "Access permission denied"; break ;
			case Error::ABORT: 		result = "Callback routine requested an abort"; break ;
			case Error::BUSY: 		result = "The database file is locked"; break ;
			case Error::LOCKED: 	result = "A table in the database is locked"; break ;
			case Error::NOMEM: 		result = "A malloc() failed"; break ;
			case Error::READONLY: 	result = "Attempt to write a readonly database"; break ;
			case Error::INTERRUPT: 	result = "Operation terminated by sqlite3_interrupt()"; break ;
			case Error::IOERR: 		result = "Some kind of disk I/O error occurred"; break ;
			case Error::CORRUPT: 	result = "The database disk image is malformed"; break ;
			case Error::NOTFOUND: 	result = "NOT USED. Table or record not found"; break ;
			case Error::FULL: 		result = "Insertion failed because database is full"; break ;
			case Error::CANTOPEN: 	result = "Unable to open the database file"; break ;
			case Error::PROTOCOL: 	result = "Database lock protocol error"; break ;
			case Error::EMPTY: 		result = "Database is empty"; break ;
			case Error::SCHEMA: 	result = "The database schema changed"; break ;
			case Error::TOOBIG: 	result = "String or BLOB exceeds size limit"; break ;
			case Error::CONSTRAINT: result = "Abort due to constraint violation"; break ;
			case Error::MISMATCH: 	result = "Data type mismatch"; break ;
			case Error::MISUSE: 	result = "Library used incorrectly"; break ;
			case Error::NOLFS: 		result = "Uses OS features not supported on host"; break ;
			case Error::AUTH: 		result = "Authorization denied"; break ;
			case Error::FORMAT: 	result = "Auxiliary database format error"; break ;
			case Error::RANGE: 		result = "2nd parameter to sqlite3_bind out of range"; break ;
			case Error::NOTADB: 	result = "File opened that is not a database file"; break ;
			case Error::ROW: 		result = "sqlite3_step() has another row ready"; break ;
			case Error::DONE: 		result = "sqlite3_step() has finished executing"; break ;
			default: 				assert(0) ;
		}
		result += ", in statement \"" + m_sql + "\"" ;
		return result ;
	}
}

