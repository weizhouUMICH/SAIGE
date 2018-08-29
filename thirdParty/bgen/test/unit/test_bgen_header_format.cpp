
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cassert>
#include "stdint.h"
#include "catch.hpp"
#include "test_utils.hpp"
#include "genfile/bgen/bgen.hpp"

namespace data {
	std::string construct_header_block(
		uint32_t number_of_snp_blocks,
		uint32_t number_of_samples,
		std::string free_data,
		uint32_t flags
	) {
		uint32_t reserved = 0x6e656762 ; // 'bgen'
		std::ostringstream oStream ;
		uint32_t header_length = free_data.size() + 20 ;
		genfile::bgen::write_little_endian_integer( oStream, header_length ) ;
		genfile::bgen::write_little_endian_integer( oStream, number_of_snp_blocks ) ;
		genfile::bgen::write_little_endian_integer( oStream, number_of_samples ) ;
		genfile::bgen::write_little_endian_integer( oStream, reserved ) ;
		oStream.write( free_data.data(), free_data.size() ) ;
		genfile::bgen::write_little_endian_integer( oStream, flags ) ;

		return oStream.str() ;
	}
}

template< typename T >
struct Setter
{
	Setter( T& field ): m_field( field ) {} ;
	void operator()( T const& value ) { m_field = value ; }
private:
	T& m_field ;
} ;

template< typename T >
Setter< T > make_setter( T& field ) { return Setter<T>( field ) ; }

void do_header_block_read_test(
	uint32_t number_of_snp_blocks,
	uint32_t number_of_samples,
	std::string free_data,
	uint32_t flags
) {
	std::istringstream inStream ;
	inStream.str(
		data::construct_header_block(
			number_of_snp_blocks,
			number_of_samples,
			free_data,
			flags
		)
	) ;

	uint32_t header_size, number_of_snp_blocks2, number_of_samples2, flags2 ;
	std::string free_data2, version ;

	genfile::bgen::Context context ;

	genfile::bgen::read_header_block(
		inStream,
		&context
	) ;
	
	REQUIRE( inStream ) ;
	REQUIRE( context.header_size() == 20 + free_data.size()) ;
	REQUIRE( context.number_of_variants == number_of_snp_blocks ) ;	
	REQUIRE( context.number_of_samples == number_of_samples ) ;	
	REQUIRE( context.free_data == free_data ) ;
	REQUIRE( context.flags == flags ) ;	
}

void do_header_block_write_test( 
	uint32_t number_of_snp_blocks,
	uint32_t number_of_samples,
	std::string free_data,
	uint32_t flags
) {
	std::ostringstream outStream ;
	genfile::bgen::Context context ;
	context.number_of_variants = number_of_snp_blocks ;
	context.number_of_samples = number_of_samples ;
	context.free_data = free_data ;
	context.flags = flags ;
	genfile::bgen::write_header_block( 
		outStream,
		context
	) ;

	std::string expected = data::construct_header_block(
		number_of_snp_blocks,
		number_of_samples,
		free_data,
		flags
	) ;
#if DEBUG
	std::cout << "\"" << to_hex( outStream.str() ) << "\"\n" ;
	std::cout << "\"" << to_hex( expected ) << "\"\n" ;
#endif
	
	REQUIRE( outStream.str() == expected ) ;
}

TEST_CASE( "Test that legal header blocks can be read", "[bgen]" ) {
	std::cout << "test_header_block_input\n" ;
	do_header_block_read_test( 100, 45, "Hi, this is some data", 1 ) ;
	do_header_block_read_test( 4294967295u, 4294967295u, "Ochen ochen, ochen ura. coheh   heheh  he  agg agg 767  $$%%$   **  ", 4294967295u ) ;
	do_header_block_read_test( 0, 0, "", 1 ) ;
	do_header_block_read_test( 0, 0, "", 0 ) ;
}

TEST_CASE( "Test that legal header blocks can be written", "[bgen]" ) {
	std::cout << "test_header_block_output\n" ;
	do_header_block_write_test( 6, 6, "This is some free text.", 1 ) ;
	do_header_block_write_test( 4294967295u, 4294967295u, "Some more free text here, you see &&&& $$$ $$$", 4294967295u ) ;
	do_header_block_write_test( 0, 0, "", 1 ) ;
}
