
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <iomanip>
#include <iostream>
#include <sstream>
#include "test_utils.hpp"

std::string to_hex( std::string const& str ) {
	return to_hex( str.data(), str.data() + str.size() ) ;
}

