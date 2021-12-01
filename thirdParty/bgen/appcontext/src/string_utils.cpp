
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <limits>
#include <vector>
#include <string>
#include <sstream>
#include <exception>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include "appcontext/string_utils.hpp"

namespace appcontext {
	namespace string_utils {
		std::string to_lower( std::string aString ) {
			for( std::string::iterator i = aString.begin(); i != aString.end(); ++i ) {
				if( *i >= 'A' && *i <= 'Z' ) {
					*i += 32 ;
				}
			}
			return aString ;
		}

		void to_lower( std::string* aString ) {
			for( std::string::iterator i = aString->begin(); i != aString->end(); ++i ) {
				if( *i >= 'A' && *i <= 'Z' ) {
					*i += 32 ;
				}
			}
		}

		std::string to_upper( std::string aString ) {
			for( std::string::iterator i = aString.begin(); i != aString.end(); ++i ) {
				if( *i >= 'a' && *i <= 'z' ) {
					*i -= 32 ;
				}
			}
			return aString ;
		}

		void to_upper( std::string* aString ) {
			for( std::string::iterator i = aString->begin(); i != aString->end(); ++i ) {
				if( *i >= 'a' && *i <= 'z' ) {
					*i -= 32 ;
				}
			}
		}

		std::string join( std::vector< std::string > const& strings, std::string const& joiner ) {
			std::string result ;
			for( std::size_t i = 0; i < strings.size(); ++i ) {
					if( i > 0 ) {
							result += joiner ;
					}
					result += strings[i] ;
			}
			return result ;
		}

		namespace impl {
			bool is_alphabetic( char c) {
					return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ;
			}

			bool is_blank( char c) {
					return c == ' ' || c == '\n' ;
			}
		}
	
	

		std::string wrap( std::string const& string_to_wrap, unsigned int wrap_column, unsigned int starting_column, std::size_t indent_amount ) { 
			assert( wrap_column > starting_column ) ; 
			assert( wrap_column > indent_amount ) ; 

			if( string_to_wrap.size() < (wrap_column - starting_column) ) { 
					return string_to_wrap ;
			}	

			std::string result ;
			unsigned int current_column = starting_column ;
			std::string indent( indent_amount, ' ' ) ; 

			std::string::const_iterator
					this_char = string_to_wrap.begin(),
					next_char = string_to_wrap.begin() ;
			for( ++next_char; this_char < string_to_wrap.end(); ++this_char, ++next_char ) { 
					if( current_column == 0 ) { 
							result.append( indent ) ; 
							current_column = indent_amount ;
							// skip whitespace at beginning of line
							for( ; this_char < string_to_wrap.end() && impl::is_blank( *this_char ); ++this_char, ++next_char ) ; 
					}	

					if( this_char < string_to_wrap.end() ) { 
							result.push_back( *this_char ) ; 
							if( *this_char == '\n' ) { 
									current_column = 0 ; 
							}	
							else if( ++current_column < wrap_column || next_char == string_to_wrap.end() ) { 
							}	
							else {
									if( impl::is_alphabetic( *this_char ) && impl::is_alphabetic( *next_char )) {
											// we are in the middle of a word.
											result.push_back( '-' ) ;
											result.push_back( '\n' ) ;
											current_column = 0 ;
									}
									else if( impl::is_blank( *this_char ) || impl::is_blank( *next_char )) {
											result.push_back( '\n' ) ;
											current_column = 0 ;
									}
									else {
											// do nothing in case this is something which can't be split.
									}
							}
					}
			}
			return result ;
		}
	}
}
