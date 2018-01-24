
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <cassert>
#include <string>
#include <stdint.h>
#include "sqlite3/sqlite3.h"
#include "db/SQLStatement.hpp"

namespace db {
	SQLStatement::~SQLStatement() {}
	
	template<>
	int SQLStatement::get_column< int >( int column_id ) const {
		return this->get_column_int( column_id ) ;
	}

	template<>
	int64_t SQLStatement::get_column< int64_t >( int column_id ) const {
		return this->get_column_int64( column_id ) ;
	}

	template<>
	double SQLStatement::get_column< double >( int column_id ) const {
		return this->get_column_double( column_id ) ;
	}

	template<>
	std::string SQLStatement::get_column< std::string >( int column_id ) const {
		return this->get_column_string( column_id ) ;
	}

	template<>
	char SQLStatement::get_column< char >( int column_id ) const {
		return this->get_column_char( column_id ) ;
	}

	template<>
	std::vector< uint8_t > SQLStatement::get_column< std::vector< uint8_t > >( int column_id ) const {
		return this->get_column_blob( column_id ) ;
	}
	
	SQLStatement& SQLStatement::bind( std::size_t i, char const* value ) {
		bind( i, std::string( value )) ;
		return *this ;
	}
}
