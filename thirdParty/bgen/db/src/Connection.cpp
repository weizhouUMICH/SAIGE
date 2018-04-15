
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "db/Connection.hpp"
#include "db/SQLite3Connection.hpp"
#include "db/SQLStatement.hpp"

namespace db {
	Connection::UniquePtr Connection::create( std::string const& filename, std::string const& mode ) {
		return Connection::UniquePtr( new SQLite3Connection( filename, true, mode )) ;
	}

	void Connection::run_statement( std::string const& SQL ) {
		get_statement( SQL )->step() ;
	}
}
