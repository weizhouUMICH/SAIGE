
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef DB_SQLITE3_STATEMENT_HPP
#define DB_SQLITE3_STATEMENT_HPP

#include <cassert>
#include <string>
#include <exception>

#include "sqlite3/sqlite3.h"
#include "db/SQLite3Connection.hpp"
#include "db/SQLStatement.hpp"

namespace db {
	class SQLite3Statement: public SQLStatement
	{
	public:
		SQLite3Statement( SQLite3Connection* connector, std::string const& SQL ) ;
		~SQLite3Statement() ;
		bool step() ;
		bool empty() const ;
		
		std::size_t get_number_of_columns() const ;
		std::string get_name_of_column( std::size_t i ) const ;
		bool is_null( int column_id ) const ;
		SQLite3Statement& bind( std::size_t i, int32_t value ) ;
		SQLite3Statement& bind( std::size_t i, uint32_t value ) ;
		SQLite3Statement& bind( std::size_t i, int64_t value ) ;
		SQLite3Statement& bind( std::size_t i, uint64_t value ) ;
		SQLite3Statement& bind( std::size_t i, double value ) ;
		SQLite3Statement& bind( std::size_t i, std::string const& value ) ;
		SQLite3Statement& bind( std::size_t i, char const* buffer, char const* const end ) ;
		SQLite3Statement& bind( std::size_t i, uint8_t const* buffer, uint8_t const* const end ) ;
		SQLite3Statement& bind_NULL( std::size_t i ) ;

		SQLite3Statement& reset() ;

		std::string get_sql() const ;

		struct Error {
			enum {
				OK			= SQLITE_OK,
				ERROR		= SQLITE_ERROR,			/* SQL error or missing database */
				INTERNAL	= SQLITE_INTERNAL,		/* Internal logic error in SQLite */
				PERM		= SQLITE_PERM,			/* Access permission denied */
				ABORT		= SQLITE_ABORT,			/* Callback routine requested an abort */
				BUSY		= SQLITE_BUSY,			/* The database file is locked */
				LOCKED		= SQLITE_LOCKED,		/* A table in the database is locked */
				NOMEM		= SQLITE_NOMEM,			/* A malloc() failed */
				READONLY	= SQLITE_READONLY,		/* Attempt to write a readonly database */
				INTERRUPT	= SQLITE_INTERRUPT,		/* Operation terminated by sqlite3_interrupt()*/
				IOERR		= SQLITE_IOERR,			/* Some kind of disk I/O error occurred */
				CORRUPT		= SQLITE_CORRUPT,		/* The database disk image is malformed */
				NOTFOUND	= SQLITE_NOTFOUND,		/* NOT USED. Table or record not found */
				FULL		= SQLITE_FULL,			/* Insertion failed because database is full */
				CANTOPEN	= SQLITE_CANTOPEN,		/* Unable to open the database file */
				PROTOCOL	= SQLITE_PROTOCOL,		/* Database lock protocol error */
				EMPTY		= SQLITE_EMPTY,			/* Database is empty */
				SCHEMA		= SQLITE_SCHEMA,		/* The database schema changed */
				TOOBIG		= SQLITE_TOOBIG,		/* String or BLOB exceeds size limit */
				CONSTRAINT	= SQLITE_CONSTRAINT,	/* Abort due to constraint violation */
				MISMATCH	= SQLITE_MISMATCH,		/* Data type mismatch */
				MISUSE		= SQLITE_MISUSE,		/* Library used incorrectly */
				NOLFS		= SQLITE_NOLFS,			/* Uses OS features not supported on host */
				AUTH		= SQLITE_AUTH,			/* Authorization denied */
				FORMAT		= SQLITE_FORMAT,		/* Auxiliary database format error */
				RANGE		= SQLITE_RANGE,			/* 2nd parameter to sqlite3_bind out of range */
				NOTADB		= SQLITE_NOTADB,		/* File opened that is not a database file */
				ROW			= SQLITE_ROW,			/* sqlite3_step() has another row ready */
				DONE		= SQLITE_DONE			/* sqlite3_step() has finished executing */
			} ;
		} ;
	protected:
		
		int64_t get_column_int64( int column_id ) const ;
		int get_column_int( int column_id ) const ;
		double get_column_double( int column_id ) const ;
		std::string get_column_string( int column_id ) const ;
		char get_column_char( int column_id ) const ;
		std::vector< uint8_t > get_column_blob( int column_id ) const ;
		
		int get_column_count() const ;

	private:
		sqlite3_stmt* m_statement ;
		SQLite3Connection* m_connection ;	
		bool m_have_results ;
	} ;

}

#endif
