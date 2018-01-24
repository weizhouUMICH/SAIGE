
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef DB_CONNECTION_HPP
#define DB_CONNECTION_HPP

#include <memory>
#include <string>
#include <stdint.h>
#include "db/Transaction.hpp"

namespace db {

	class SQLStatement ;
	
	class Connection
		// Base class for classes representing a connection to a database.
		// The only supported operation is getting a query representing some SQL.
	{
	public:
		typedef std::auto_ptr< Connection > UniquePtr ;
		typedef std::auto_ptr< SQLStatement > StatementPtr ;
		
		static UniquePtr create( std::string const& filename, std::string const& mode = "rw" ) ;
		
		virtual ~Connection() {}

		virtual StatementPtr get_statement( std::string const& SQL ) = 0 ;
		virtual void run_statement( std::string const& SQL ) ;
		virtual std::string get_spec() const = 0 ;
		typedef int64_t RowId ;
		virtual RowId get_last_insert_row_id() const = 0 ;
		typedef Transaction::UniquePtr ScopedTransactionPtr ;
		virtual ScopedTransactionPtr open_transaction( double max_seconds_to_wait = 0.1 ) = 0 ;
	} ;	
}

#endif
