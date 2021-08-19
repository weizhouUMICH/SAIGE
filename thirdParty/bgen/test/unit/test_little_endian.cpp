
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

namespace {
	genfile::byte_t buffer[100] ;
	template< typename Integer >
	void test_rw( Integer value ) {
		genfile::byte_t const* end = genfile::bgen::write_little_endian_integer( buffer, buffer + 100, value ) ;

		// Test most and least significant bytes stored in the right place
		std::size_t const lsb = 0 ;
		std::size_t const msb = sizeof(Integer)-1 ;
		REQUIRE(( static_cast< unsigned char >( buffer[lsb] ) == static_cast< unsigned char>(value & 0xFF) )) ;
		REQUIRE(( static_cast< unsigned char >( buffer[msb] ) == static_cast< unsigned char >( (value>>(msb*8)) & 0xFF) )) ;

		// Test we can reconstruct it
		Integer result ;
		genfile::bgen::read_little_endian_integer( buffer, end, &result ) ;
		REQUIRE(( result == value )) ;
	}
}

TEST_CASE( "Test that integers can be written and recovered from a buffer." ) {
	std::cout << "test_little_endian_integer_io\n" ;
	genfile::byte_t buffer[100] ;

	test_rw<char>( 0 ) ;
	test_rw<char>( 1 ) ;
	test_rw<char>( 127 ) ;
	test_rw<char>( -128 ) ;
	test_rw<char>( std::numeric_limits<char>::max() ) ;
	test_rw<char>( std::numeric_limits<char>::min() ) ;

	test_rw<unsigned char>( 0 ) ;
	test_rw<unsigned char>( 1 ) ;
	test_rw<unsigned char>( 127 ) ;
	test_rw<unsigned char>( 255 ) ;
	test_rw<unsigned char>( std::numeric_limits<unsigned char>::max() ) ;
	test_rw<unsigned char>( std::numeric_limits<unsigned char>::min() ) ;

	test_rw<uint16_t>( 0 ) ;
	test_rw<uint16_t>( 1 ) ;
	test_rw<uint16_t>( 255 ) ;
	test_rw<uint16_t>( std::numeric_limits<uint16_t>::max() ) ;
	test_rw<uint16_t>( std::numeric_limits<uint16_t>::min() ) ;

	test_rw<int16_t>( 0 ) ;
	test_rw<int16_t>( 1 ) ;
	test_rw<int16_t>( -1 ) ;
	test_rw<int16_t>( 255 ) ;
	test_rw<int16_t>( std::numeric_limits<int16_t>::max() ) ;
	test_rw<int16_t>( std::numeric_limits<int16_t>::min() ) ;

	test_rw<uint32_t>( 0 ) ;
	test_rw<uint32_t>( 1 ) ;
	test_rw<uint32_t>( 255 ) ;
	test_rw<uint32_t>( std::numeric_limits<uint32_t>::max() ) ;
	test_rw<uint32_t>( std::numeric_limits<uint32_t>::min() ) ;

	test_rw<int32_t>( 0 ) ;
	test_rw<int32_t>( 1 ) ;
	test_rw<int32_t>( -1 ) ;
	test_rw<int32_t>( 255 ) ;
	test_rw<int32_t>( std::numeric_limits<int32_t>::max() ) ;
	test_rw<int32_t>( std::numeric_limits<int32_t>::min() ) ;

	test_rw<uint64_t>( 0 ) ;
	test_rw<uint64_t>( 1 ) ;
	test_rw<uint64_t>( 255 ) ;
	test_rw<uint64_t>( std::numeric_limits<uint64_t>::max() ) ;
	test_rw<uint64_t>( std::numeric_limits<uint64_t>::min() ) ;

	test_rw<int64_t>( 0 ) ;
	test_rw<int64_t>( 1 ) ;
	test_rw<int64_t>( -1 ) ;
	test_rw<int64_t>( 255 ) ;
	test_rw<int64_t>( std::numeric_limits<int64_t>::max() ) ;
	test_rw<int64_t>( std::numeric_limits<int64_t>::min() ) ;
}
