
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <cassert>
#include <string>
#include <exception>
#include <boost/lexical_cast.hpp>
#include "sqlite3/sqlite3.h"
#include "db/SQLite3Connection.hpp"
#include "db/SQLStatement.hpp"
#include "db/SQLite3Statement.hpp"
#include "db/SQLite3Error.hpp"

namespace db {
	SQLite3Statement::SQLite3Statement( SQLite3Connection* connection, std::string const& SQL ):
		m_connection( connection ),
		m_have_results( false )
	{
		assert( m_connection ) ;
		m_statement = m_connection->prepare_sql( SQL ) ;
		//std::cerr << "SQLite3Statement::SQLite3Statement(): statement is \"" + SQL + "\".\n" ;
		assert( m_statement != 0 ) ;
	}
	
	SQLite3Statement::~SQLite3Statement() {
		m_connection->finalise_statement( m_statement ) ;
	}

	bool SQLite3Statement::step() {
		m_have_results = m_connection->step_statement( m_statement ) ;
		return m_have_results ;
	}

	bool SQLite3Statement::empty() const {
		return !m_have_results ;
	}

	std::size_t SQLite3Statement::get_number_of_columns() const {
		return std::size_t( get_column_count() ) ;
	}

	std::string SQLite3Statement::get_name_of_column( std::size_t i ) const {
		assert( m_statement != 0 ) ;
		assert( i < get_number_of_columns() ) ;
		return std::string( sqlite3_column_origin_name( m_statement, i )) ;
	}

	SQLite3Statement& SQLite3Statement::bind( std::size_t i, int32_t value ) {
		assert( m_statement != 0 ) ;
		int error = sqlite3_bind_int( m_statement, i, value ) ;
		if( error != SQLITE_OK ) {
			throw ValueBindError( "SQLite3Statement::bind()", m_connection->get_spec(), error, boost::lexical_cast<std::string>( i ) ) ;
		}
		return *this ;
	}

	SQLite3Statement& SQLite3Statement::bind( std::size_t i, uint32_t value ) {
		assert( m_statement != 0 ) ;
		int error = sqlite3_bind_int( m_statement, i, value ) ;
		if( error != SQLITE_OK ) {
			throw ValueBindError( "SQLite3Statement::bind()", m_connection->get_spec(), error, boost::lexical_cast<std::string>( i ) ) ;
		}
		return *this ;
	}

	SQLite3Statement& SQLite3Statement::bind( std::size_t i, int64_t value ) {
		assert( m_statement != 0 ) ;
		int error = sqlite3_bind_int64( m_statement, i, sqlite3_int64( value ) ) ;
		if( error != SQLITE_OK ) {
			throw ValueBindError( "SQLite3Statement::bind()", m_connection->get_spec(), error, boost::lexical_cast<std::string>( i ) ) ;
		}
		return *this ;
	}

	SQLite3Statement& SQLite3Statement::bind( std::size_t i, uint64_t value ) {
		assert( m_statement != 0 ) ;
		int error = sqlite3_bind_int64( m_statement, i, sqlite3_int64( value ) ) ;
		if( error != SQLITE_OK ) {
			throw ValueBindError( "SQLite3Statement::bind()", m_connection->get_spec(), error, boost::lexical_cast<std::string>( i ) ) ;
		}
		return *this ;
	}

	SQLite3Statement& SQLite3Statement::bind( std::size_t i, double value ) {
		assert( m_statement != 0 ) ;
		int error = sqlite3_bind_double( m_statement, i, value ) ;
		if( error != SQLITE_OK ) {
			throw ValueBindError( "SQLite3Statement::bind()", m_connection->get_spec(), error, boost::lexical_cast<std::string>( i ) ) ;
		}
		return *this ;
	}

	SQLite3Statement& SQLite3Statement::bind( std::size_t i, std::string const& value ) {
		assert( m_statement != 0 ) ;
		int error = sqlite3_bind_text( m_statement, i, value.c_str(), value.size(), SQLITE_TRANSIENT ) ;
		if( error != SQLITE_OK ) {
			throw ValueBindError( "SQLite3Statement::bind()", m_connection->get_spec(), error, boost::lexical_cast<std::string>( i ) ) ;
		}
		return *this ;
	}

	SQLite3Statement& SQLite3Statement::bind( std::size_t i, char const* buffer, char const* const end ) {
		assert( m_statement != 0 ) ;
		int error = sqlite3_bind_blob( m_statement, i, reinterpret_cast< void const* >( buffer ), int( end - buffer ), SQLITE_TRANSIENT ) ;
		if( error != SQLITE_OK ) {
			throw ValueBindError( "SQLite3Statement::bind()", m_connection->get_spec(), error, boost::lexical_cast<std::string>( i ) ) ;
		}
		return *this ;
	}

	SQLite3Statement& SQLite3Statement::bind( std::size_t i, uint8_t const* buffer, uint8_t const* const end ) {
		assert( m_statement != 0 ) ;
		int error = sqlite3_bind_blob( m_statement, i, reinterpret_cast< void const* >( buffer ), int( end - buffer ), SQLITE_TRANSIENT ) ;
		if( error != SQLITE_OK ) {
			throw ValueBindError( "SQLite3Statement::bind()", m_connection->get_spec(), error, boost::lexical_cast<std::string>( i ) ) ;
		}
		return *this ;
	}

	SQLite3Statement& SQLite3Statement::bind_NULL( std::size_t i ) {
		assert( m_statement != 0 ) ;
		int error = sqlite3_bind_null( m_statement, i ) ;
		if( error != SQLITE_OK ) {
			throw ValueBindError( "SQLite3Statement::bind()", m_connection->get_spec(), error, boost::lexical_cast<std::string>( i ) ) ;
		}
		return *this ;
	}

	SQLite3Statement& SQLite3Statement::reset() {
		assert( m_statement != 0 ) ;
		int error = sqlite3_reset( m_statement ) ;
		if( error != SQLITE_OK ) {
			throw db::Error( "SQLite3Statement::reset()", m_connection->get_spec(), error ) ;
		}
		return *this ;
	}

	std::string SQLite3Statement::get_sql() const {
		assert( m_statement != 0 ) ;
		return std::string( sqlite3_sql( m_statement )) ;
	}

	bool SQLite3Statement::is_null( int column_id ) const {
		return ( sqlite3_column_type( m_statement, column_id ) == SQLITE_NULL ) ;
	}

	int SQLite3Statement::get_column_int( int column_id ) const {
		assert( m_statement != 0 ) ;
		return sqlite3_column_int( m_statement, column_id ) ;
	}

	int64_t SQLite3Statement::get_column_int64( int column_id ) const {
		assert( m_statement != 0 ) ;
		return sqlite3_column_int64( m_statement, column_id ) ;
	}

	double SQLite3Statement::get_column_double( int column_id ) const {
		assert( m_statement != 0 ) ;
		return sqlite3_column_double( m_statement, column_id ) ;
	}
	
	std::string SQLite3Statement::get_column_string( int column_id ) const {
		assert( m_statement != 0 ) ;
		return reinterpret_cast< char const * >( sqlite3_column_text( m_statement, column_id ) ) ;
	}

	char SQLite3Statement::get_column_char( int column_id ) const {
		assert( m_statement != 0 ) ;
		char const* p = reinterpret_cast< char const * >( sqlite3_column_text( m_statement, column_id )) ;
		int bytes = sqlite3_column_bytes( m_statement, column_id ) ;
		assert( bytes == 1 ) ;
		return *p ;
	}
	
	std::vector< uint8_t > SQLite3Statement::get_column_blob( int column_id ) const {
		assert( m_statement != 0 ) ;
		uint8_t const* p = reinterpret_cast< uint8_t const* >( sqlite3_column_blob( m_statement, column_id )) ;
		int nBytes = sqlite3_column_bytes( m_statement, column_id ) ;
		return std::vector< uint8_t >( p, p+nBytes ) ;
	}
	
	int SQLite3Statement::get_column_count() const {
		assert( m_statement != 0 ) ;
		return sqlite3_column_count( m_statement ) ;
	}
}
