
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef APPCONTEXT_STRING_UTILS_HPP
#define APPCONTEXT_STRING_UTILS_HPP

#include <vector>
#include <string>
#include <sstream>
#include <exception>

namespace appcontext {
	namespace string_utils {
		std::string to_lower( std::string aString ) ;
		void to_lower( std::string* aString ) ;
		std::string to_upper( std::string aString ) ;
		void to_upper( std::string* aString ) ;

		std::string join( std::vector< std::string > const& strings, std::string const& joiner ) ;
		std::string wrap( std::string const& string_to_wrap, unsigned int wrap_column, unsigned int starting_column, std::size_t indent_amount ) ;
	}
}

#endif
