
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef NULL_OSTREAM_HPP
#define NULL_OSTREAM_HPP

#include <iostream>

struct null_ostream: public std::ostream
{
	struct null_streambuf: public std::streambuf
	{
		int overflow( int c ) { return std::char_traits< char >::not_eof(c) ; }
	} ;
	
	null_ostream(): std::ios( &m_buf), std::ostream( &m_buf ) {}

private:
	null_streambuf m_buf ;
} ;

#endif
