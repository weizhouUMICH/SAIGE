
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_ZLIB_HPP
#define GENFILE_ZLIB_HPP

#include <vector>
#include <stdint.h>
#include <cassert>
#include <zlib.h>
#include "zstd.h"
#include "genfile/types.hpp"

namespace genfile {

	// Compress the given data into the given destination buffer.  The destination will be resized
	// to fit the compressed data.  (Since the capacity of dest may be larger than its size,
	// to save memory you may need to copy the contents of dest elsewhere after calling
	// this function).
	//
	// If offset is nonzero, compressed data will be written starting at position [offset].
	// The first [offset] bytes will be untouched.
	void zlib_compress(
		byte_t const* buffer,
		byte_t const* const end,
		std::vector< byte_t >* dest,
		std::size_t const offset = 0,
		int const compressionLevel = Z_BEST_COMPRESSION
	) ;

	void zstd_compress(
		byte_t const* buffer,
		byte_t const* const end,
		std::vector< byte_t >* dest,
		std::size_t const offset = 0,
		int const compressionLevel = 22
	) ;

	// Compress the given data into the given destination buffer.  The destination will be resized
	// to fit the compressed data.  (Since the capacity of dest may be larger than its size,
	// to save memory you may need to copy the contents of dest elsewhere after calling
	// this function).
	template< typename T >
	void zlib_compress(
		std::vector< T > const& source,
		std::vector< byte_t >* dest,
		int const compressionLevel = Z_BEST_COMPRESSION
	 ) {
		byte_t const* begin = reinterpret_cast< byte_t const* >( &source[0] ) ;
		byte_t const* const end = reinterpret_cast< byte_t const* >( &source[0] + source.size() ) ;
		return zlib_compress( begin, end, dest, 0, compressionLevel ) ;
	}

	template< typename T >
	void zlib_uncompress(
		byte_t const* begin,
		byte_t const* const end,
		std::vector< T >* dest
	) {
		uLongf const source_size = ( end - begin ) ;
		uLongf dest_size = dest->size() * sizeof( T ) ;
		int result = uncompress(
			reinterpret_cast< Bytef* >( &dest->operator[]( 0 ) ),
			&dest_size,
			reinterpret_cast< Bytef const* >( begin ),
			source_size
		) ;
		assert( result == Z_OK ) ;
		assert( dest_size % sizeof( T ) == 0 ) ;
		dest->resize( dest_size / sizeof( T )) ;
	}

	template< typename T >
	void zstd_uncompress( byte_t const* begin, byte_t const* const end, std::vector< T >* dest ) {
		std::size_t const source_size = ( end - begin ) ;
		std::size_t const dest_size = dest->size() * sizeof( T ) ;
	    std::size_t const uncompressed_size = ZSTD_getDecompressedSize( reinterpret_cast< void const* >( begin ), source_size ) ;
		std::size_t const result = ZSTD_decompress(
			reinterpret_cast< void* >( &dest->operator[]( 0 ) ),
			dest_size,
			reinterpret_cast< void const* >( begin ),
			source_size
		) ;
		assert( result == uncompressed_size ) ;
		dest->resize( dest_size / sizeof( T )) ;
	}

	// Uncompress the given data, symmetric with zlib_compress.
	// The destination must be large enough to fit the uncompressed data,
	// and it will be resized to exactly fit the uncompressed data.
	template< typename T >
	void zlib_uncompress( std::vector< byte_t > const& source, std::vector< T >* dest ) {
		byte_t const* begin = &source[0] ;
		byte_t const* const end = &source[0] + source.size() ;
		zlib_uncompress( begin, end, dest ) ;
	}
}

#endif
