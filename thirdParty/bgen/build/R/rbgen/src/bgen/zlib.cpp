
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <cassert>
#include <zlib.h>
#include "genfile/zlib.hpp"

namespace genfile {
	void zlib_compress(
		uint8_t const* buffer,
		uint8_t const* const end,
		std::vector< uint8_t >* dest,
		std::size_t const offset,
		int const compressionLevel
	) {
		assert( dest != 0 ) ;
		assert( compressionLevel >= 0 && compressionLevel <= Z_BEST_COMPRESSION ) ;
		uLongf const source_size = ( end - buffer ) ;
		uLongf compressed_size = compressBound( source_size ) ;
		dest->resize( compressed_size + offset ) ;
		int result = compress2(
			reinterpret_cast< Bytef* >( const_cast< uint8_t* >( &( dest->operator[](0) ) + offset ) ),
			&compressed_size,
			reinterpret_cast< Bytef const* >( buffer ),
			source_size,
			compressionLevel
		) ;
		assert( result == Z_OK ) ;
		dest->resize( compressed_size + offset ) ;
	}

	void zstd_compress(
		uint8_t const* buffer,
		uint8_t const* const end,
		std::vector< uint8_t >* dest,
		std::size_t const offset,
		int const compressionLevel
	) {
		assert( dest != 0 ) ;
		assert( compressionLevel >= 1 && compressionLevel <= 22 ) ;
		std::size_t const source_size = ( end - buffer ) ;
		std::size_t compressed_size = ZSTD_compressBound( source_size ) ;
		dest->resize( compressed_size + offset ) ;
		compressed_size = ZSTD_compress(
			reinterpret_cast< void* >( &( dest->operator[](0) ) + offset ),
			compressed_size,
			reinterpret_cast< void const* >( buffer ),
			source_size,
			compressionLevel
		) ;
		assert( !ZSTD_isError( compressed_size )) ;
		dest->resize( compressed_size + offset ) ;
	}
}
