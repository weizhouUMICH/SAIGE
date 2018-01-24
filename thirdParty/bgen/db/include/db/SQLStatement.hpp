
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef DB_SQL_STATEMENT_HPP
#define DB_SQL_STATEMENT_HPP

#include <cassert>
#include <string>
#include <vector>
#include <exception>
#include <stdint.h>
#include "sqlite3/sqlite3.h"
#include "db/SQLite3Connection.hpp"

namespace db {
	
	struct SQLError: public std::exception
	{
		virtual ~SQLError() throw() {}
		char const* what() const throw() { return "db::SQLError" ; }
		virtual std::string description() const = 0 ;
	} ;
	
	class SQLStatement
	{
	public:
		typedef std::auto_ptr< SQLStatement > UniquePtr ;
	public:
		virtual ~SQLStatement() ;
		
		// Step to the next result row.  The entries can be accessed using get_column().
		// Return false if there is no next row, otherwise true.
		virtual bool step() = 0 ;
		// Return 
		virtual bool empty() const = 0 ;
		
		// Get the result for the given column.
		// Columns are 0-indexed.
		template< typename T > T get_column( int column_id ) const ;
		template< typename T > T get( int column_id ) const { return get_column< T >( column_id ) ; }
		virtual bool is_null( int column_id ) const = 0 ;

		virtual std::size_t get_number_of_columns() const = 0 ;
		virtual std::string get_name_of_column( std::size_t i ) const = 0 ;

		// For parameterised queries, bind an integer value to the ith placeholder.
		// Placeholders are indexed starting from 1 on the left.
		// Named placeholders with the same name have the same index.
		virtual SQLStatement& bind( std::size_t i, int32_t value ) = 0 ;
		virtual SQLStatement& bind( std::size_t i, uint32_t value ) = 0 ;
		virtual SQLStatement& bind( std::size_t i, int64_t value ) = 0 ;
		virtual SQLStatement& bind( std::size_t i, uint64_t value ) = 0 ;
		virtual SQLStatement& bind( std::size_t i, double value ) = 0 ;
		// For parameterised queries, bind a string value to the ith placeholder.
		// Placeholders are indexed starting from 1 on the left.
		// Named placeholders with the same name have the same index.
		// The string will be copied so the caller need not preserve it beyond the call site.
		virtual SQLStatement& bind( std::size_t i, std::string const& value ) = 0 ;
		SQLStatement& bind( std::size_t i, char const* value ) ;
		//SQLStatement& bind( std::size_t i, genfile::string_utils::slice const& value ) ;
		// For parameterised queries, bind a BLOB value (array of chars) to the ith placeholder.
		// Placeholders are indexed starting from 1 on the left.
		// Named placeholders with the same name have the same index.
		// The data will not be copied and so the caller must preserve the data until
		// such time as no further steps() are performed, or the parameter is re-bound.
		virtual SQLStatement& bind( std::size_t i, char const* buffer, char const* const end ) = 0 ;
		virtual SQLStatement& bind( std::size_t i, uint8_t const* buffer, uint8_t const* const end ) = 0 ;

		// bind NULL to a parameter
		virtual SQLStatement& bind_NULL( std::size_t i ) = 0 ;

		// Bind a genfile::VariantEntry
		// SQLStatement& bind( std::size_t, genfile::VariantEntry const& value ) ;

		// Reset the statement, ready to be re-executed.
		virtual SQLStatement& reset() = 0 ;
		
		// Return the SQL this statement contains.
		virtual std::string get_sql() const = 0 ;
		
	protected:
		virtual int get_column_int( int column_id ) const = 0 ;
		virtual int64_t get_column_int64( int column_id ) const = 0 ;
		virtual double get_column_double( int column_id ) const = 0 ;
		virtual std::string get_column_string( int column_id ) const = 0 ;
		virtual char get_column_char( int column_id ) const = 0 ;
		virtual std::vector< uint8_t > get_column_blob( int column_id ) const = 0 ;
	} ;
	
	template<> int SQLStatement::get_column< int >( int column_id ) const ;
	template<> int64_t SQLStatement::get_column< int64_t >( int column_id ) const ;
	template<> double SQLStatement::get_column< double >( int column_id ) const ;
	template<> std::string SQLStatement::get_column< std::string >( int column_id ) const ;
}

#endif
