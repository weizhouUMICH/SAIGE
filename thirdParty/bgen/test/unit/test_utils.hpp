
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <memory>

std::string to_hex( std::string const& str ) ;
template< typename Iterator, typename End >
std::string to_hex( Iterator begin, End end ) {
	std::ostringstream o ;
	for( Iterator p = begin; p != end; ++p ) {
		if( (p-begin) % 4 == 0 )
			o << "|" ;
		o << std::hex << std::setw(2) << std::setfill('0') << static_cast<int> ( static_cast<unsigned char>( *p ) ) ;
	}
	return o.str() ;
}

	

template<typename ... Args>
std::string format( std::string const& format, Args... args )
{
	size_t size = std::snprintf( nullptr, 0, format.c_str(), args... ) + 1; // Extra space for '\0'
	std::unique_ptr<char[]> buf( new char[ size ] ) ; 
	std::snprintf( buf.get(), size, format.c_str(), args... );
	return std::string( buf.get(), buf.get() + size - 1 );
}

