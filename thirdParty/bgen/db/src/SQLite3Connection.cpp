
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <stdint.h>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/thread/thread_time.hpp>
#include <boost/thread/thread.hpp>
#include "db/SQLite3Connection.hpp"
#include "db/SQLite3Statement.hpp"
#include "db/Error.hpp"

extern "C" {
	void sqlite3_trace_callback( void* udp, const char* sql ) {
		std::cerr << "SQLite3 trace: SQL = \"" << sql << "\".\n" ;
	}
	
	int sqlite3_busy_callback( void*, int number_of_tries ) {
		if( number_of_tries > 10 ) {
			return 0 ;
		}
		boost::this_thread::sleep( boost::posix_time::milliseconds( 10 ) ) ;
		return 1 ;
	}
}

namespace db {
	SQLite3Connection::SQLite3Connection( std::string const& filename, bool overwrite, std::string const& mode ):
		m_filename( filename ),
		m_db_connection(0),
		m_managed( true )
	{
		open_db_connection( filename, overwrite, mode ) ;
	}

	SQLite3Connection::~SQLite3Connection() {
		close_db_connection_if_necessary() ;
	}

	
	SQLite3Connection::StatementPtr SQLite3Connection::get_statement( std::string const& SQL ) {
		return StatementPtr( new SQLite3Statement( this, SQL ) ) ;
	}	

	Connection::RowId SQLite3Connection::get_last_insert_row_id() const {
		uint64_t result = sqlite3_last_insert_rowid( m_db_connection ) ;
		return result ;
	}

	sqlite3_stmt* SQLite3Connection::prepare_sql( std::string const& SQL ) const {
		assert( m_db_connection != 0 ) ;
		sqlite3_stmt* statement ;
		int code = sqlite3_prepare_v2(
			m_db_connection,
			SQL.c_str(),
			static_cast< int >( SQL.size()+1 ), // +1 accounts for null terminating byte.
			&statement,
			0 // ignore pzTail
		) ;
		if( code != SQLITE_OK ) {
			throw StatementPreparationError( "SQLite3Connection::prepare_sql()", get_spec(), code, SQL ) ;
		}
		// SQLite might return 0 if the SQL consisted of comments and whitespace only
		// but I want to treat this as a programmer error.
		assert( statement != 0 ) ;
		return statement ;
	}

	int SQLite3Connection::finalise_statement( sqlite3_stmt* statement ) {
		assert( statement != 0 ) ;
		return sqlite3_finalize( statement ) ;
	}

	int SQLite3Connection::step_statement( sqlite3_stmt* statement ) {
		int code = sqlite3_step( statement ) ;
		if( code != SQLITE_ROW && code != SQLITE_DONE ) {
			throw StatementStepError( "SQLite3Connection::step_statement()", get_spec(), code, std::string( sqlite3_sql( statement ) ) ) ;
		}
		return (code == SQLITE_ROW);
	}
	
	void SQLite3Connection::open_db_connection( std::string const& filename, bool overwrite, std::string const& mode ) {
		int flags = 0 ;
		if( mode == "r" ) {
			flags |= SQLITE_OPEN_READONLY ;
		} else if( mode == "rw" ) {
			flags |= SQLITE_OPEN_READWRITE ;
			if( overwrite ) {
				flags |= SQLITE_OPEN_CREATE ;
			}
		} else {
			assert( 0 ) ;
		}
		int code = sqlite3_open_v2( filename.c_str(), &m_db_connection, flags, NULL ) ;
		if( code != SQLITE_OK ) {
			throw ConnectionError( "SQLite3Connection::open_db_connection()", get_spec(), code ) ;
		}
		// We add a busy handler.  This makes the database more robust by retrying failed transactions.
		sqlite3_busy_handler( m_db_connection, &sqlite3_busy_callback, NULL ) ;
		// Uncomment the next line to trace SQL statements executed.
		//sqlite3_trace( m_db_connection, &sqlite3_trace_callback, NULL ) ;
	}

	void SQLite3Connection::close_db_connection_if_necessary() {
		if( m_managed && m_db_connection != 0 ) {
			// According to SQLite docs, we must finalise any prepared statements
			// before we can close the db connection.
			finalise_prepared_statements() ;
			sqlite3_close( m_db_connection ) ;
			m_db_connection = 0 ;
		}
	}

	void SQLite3Connection::finalise_prepared_statements() {
		assert( m_db_connection != 0 ) ;
#if SQLITE_VERSION_NUMBER > 3006000
		sqlite3_stmt* statement ;
		while(( statement = sqlite3_next_stmt( m_db_connection, 0)) != 0 ) {
			finalise_statement( statement ) ;
		}
#endif
	}
	
	SQLite3Connection::ScopedTransactionPtr SQLite3Connection::open_transaction( double max_seconds_to_wait ) {		
		ScopedTransactionPtr transaction ;
		// we wait 10 milliseconds between attempts.
		
		for( std::size_t count = 0 ; true; ++count ) {
			try {
				transaction.reset( new SQLite3Connection::Transaction( *this ) ) ;
				break ;
			}
			catch( db::StatementStepError const& e ) {
				// Because of the busy handler (see top of this file)
				// each attempt takes ~0.1s anyway
				// We wait an additional 0.1s so that each attempt takes 0.2s in total.
				if( ( count * 0.2 ) > max_seconds_to_wait ) {
					std::cerr << "Open transaction: failure count=" << count << " (~" << count*0.2 << "s).  Bailing out.\n" ;
					boost::this_thread::sleep( boost::posix_time::milliseconds( 100 ) ) ;
					throw TransactionError( "SQLite3Connection::open_transaction()", get_spec(), e.error_code(), e.sql() ) ;
				}
			}
		}
		return transaction ;
	}
	
	SQLite3Connection::Transaction::Transaction( SQLite3Connection& connection ):
		m_connection( connection )
	{
		m_connection.run_statement( "BEGIN IMMEDIATE TRANSACTION" ) ;
	}
	
	SQLite3Connection::Transaction::~Transaction() {
		m_connection.run_statement( "COMMIT" ) ;
	}
}
