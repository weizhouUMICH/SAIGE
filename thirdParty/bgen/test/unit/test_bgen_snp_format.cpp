
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cassert>
#include <cstring>
#include <map>
#include <set>
#include <functional>
#include <cstdio>
#include <cstring>
#include "stdint.h"
#include "catch.hpp"
#include "test_utils.hpp"
#include "genfile/bgen/bgen.hpp"
#include "genfile/types.hpp"
#include "ProbSetCheck.hpp"

using namespace std::placeholders ;

// #define DEBUG 1

// The following section contains a simple snp block writer.
namespace data {
	std::string construct_snp_block_v10(
		uint32_t number_of_samples,
		unsigned char max_id_size,
		std::string SNPID,
		std::string RSID,
		std::string chromosome,
		uint32_t SNP_position,
		std::string a_allele,
		std::string b_allele,
		std::function< double ( std::size_t i, std::size_t g ) > get_probs
	) {
		std::ostringstream oStream ;
		genfile::bgen::write_little_endian_integer( oStream, number_of_samples ) ;
		genfile::bgen::write_little_endian_integer( oStream, max_id_size ) ;
		genfile::bgen::write_little_endian_integer( oStream, static_cast< char >( SNPID.size() )) ;
		oStream.write( SNPID.data(), SNPID.size() ) ;
		oStream.write( "                ", max_id_size - SNPID.size()) ;
		genfile::bgen::write_little_endian_integer( oStream, static_cast< char >( RSID.size() )) ;
		oStream.write( RSID.data(), RSID.size() ) ;
		oStream.write( "                ", max_id_size - RSID.size()) ;
		unsigned char chr = std::stoi( chromosome ) & 0xFF ;
		genfile::bgen::write_little_endian_integer( oStream, chr ) ;
		genfile::bgen::write_little_endian_integer( oStream, SNP_position ) ;

		assert( a_allele.size() == 1 ) ;
		assert( b_allele.size() == 1 ) ;
		unsigned char const a_allele_size = a_allele.size() ;
		unsigned char const b_allele_size = a_allele.size() ;
		oStream.put( a_allele_size ) ;
		oStream.write( a_allele.data(), a_allele.size() ) ;
		oStream.put( b_allele_size ) ;
		oStream.write( b_allele.data(), b_allele.size() ) ;
		
		for( std::size_t i = 0; i < number_of_samples; ++i ) {
			// We store the probability
			// ((i*100)+o) / 10000.0
			uint16_t
				AA = std::floor( 0.5 + get_probs( i, 0 ) * 10000 ),
				AB = std::floor( 0.5 + get_probs( i, 1 ) * 10000 ),
				BB = std::floor( 0.5 + get_probs( i, 2 ) * 10000 )
			;
			genfile::bgen::write_little_endian_integer( oStream, AA ) ;
			genfile::bgen::write_little_endian_integer( oStream, AB ) ;
			genfile::bgen::write_little_endian_integer( oStream, BB ) ;
		}

		return oStream.str() ;
	}

	std::string construct_snp_block_v11(
		uint32_t number_of_samples,
		std::string SNPID,
		std::string RSID,
		std::string chromosome,
		uint32_t SNP_position,
		std::string a_allele,
		std::string b_allele,
		std::function< double ( std::size_t i, std::size_t g ) > get_probs
	) {
		std::ostringstream oStream ;
		genfile::bgen::write_little_endian_integer( oStream, number_of_samples ) ;
		genfile::bgen::write_little_endian_integer( oStream, static_cast< uint16_t >( SNPID.size() )) ;
		oStream.write( SNPID.data(), SNPID.size() ) ;
		genfile::bgen::write_little_endian_integer( oStream, static_cast< uint16_t >( RSID.size() )) ;
		oStream.write( RSID.data(), RSID.size() ) ;
		genfile::bgen::write_little_endian_integer( oStream, static_cast< uint16_t >( chromosome.size() ) ) ;
		oStream.write( chromosome.data(), chromosome.size() ) ;
		genfile::bgen::write_little_endian_integer( oStream, SNP_position ) ;

		genfile::bgen::write_little_endian_integer( oStream, static_cast< uint32_t >( a_allele.size() )) ;
		oStream.write( a_allele.data(), a_allele.size() ) ;
		genfile::bgen::write_little_endian_integer( oStream, static_cast< uint32_t >( b_allele.size() )) ;
		oStream.write( b_allele.data(), b_allele.size() ) ;
		
		for( std::size_t i = 0; i < number_of_samples; ++i ) {
			uint16_t
				AA = std::floor( 0.5 + get_probs( i, 0 ) * 32768.0 ),
				AB = std::floor( 0.5 + get_probs( i, 1 ) * 32768.0 ),
				BB = std::floor( 0.5 + get_probs( i, 2 ) * 32768.0 )
			;
			genfile::bgen::write_little_endian_integer( oStream, AA ) ;
			genfile::bgen::write_little_endian_integer( oStream, AB ) ;
			genfile::bgen::write_little_endian_integer( oStream, BB ) ;
		}

		return oStream.str() ;
	}
	
	void round_probs_to_simplex( double *p, std::size_t n, std::size_t number_of_bits ) {
		double const scale = uint64_t( 0xFFFFFFFFFFFFFFFF ) >> ( 64 - number_of_bits ) ;
		std::multimap< double, double*, std::greater< double > > fractional_parts ;
		double total_fractional_part = 0.0 ;
		for( std::size_t i = 0; i < n; ++i ) {
			*(p+i) *= scale ;
			double const fractional_part = *(p+i) - std::floor(*(p+i)) ;
			fractional_parts.insert( std::make_pair( fractional_part, (p+i) ) ) ;
			total_fractional_part += fractional_part ;
		}
		std::size_t const upper = std::floor( 0.5 + total_fractional_part ) ;
#if DEBUG > 2
		std::cerr << "round_probs_to_simplex(): number_of_bits = " << number_of_bits << ", scale = " << scale << ", total_fractional_part = " << total_fractional_part << ", upper = " << upper << ".\n" ;
		std::cerr << "round_probs_to_simplex(): p1 = " << *p << ".\n" ;
#endif
		std::multimap< double, double*, std::greater< double > >::const_iterator
			i = fractional_parts.begin(),
			end_i = fractional_parts.end() ;
		for( std::size_t count = 0; i != end_i; ++i, ++count ) {
			if( count < upper ) {
				*(i->second) = std::ceil( *(i->second) ) / scale ;
			} else {
				*(i->second) = std::floor( *(i->second) ) / scale ;
			}
		}
	}

	void write_probs( uint64_t** p, std::size_t* offset, uint64_t* const end_p, double const* probs, std::size_t const n, std::size_t number_of_bits ) {
		double const scale = uint64_t( 0xFFFFFFFFFFFFFFFF ) >> ( 64 - number_of_bits ) ; // 2^bits - 1.
		for( std::size_t i = 0; i < (n-1); ++i ) {
#if DEBUG > 1
			std::cerr << "write_probs(): p = " << std::hex << (*p)
				<< std::dec << ", offset = " << (*offset) << ", end_p = " << (end_p)
				<< ", number_of_bits = " << number_of_bits
				<< ", value = " << probs[i]
				<< ", unscaled value = " << uint64_t( probs[i] * scale )
				<< ".\n";
#endif
			uint64_t const storedValue = uint64_t( probs[i] * scale ) ;
			**p |= storedValue << *offset ;
			(*offset) += number_of_bits ;
			if( (*offset) >= 64 ) {
				// move to next 64-bit word
				(*offset) -= 64 ;
				++(*p) ;
				// Make sure and store unwritten remainder of value...
				(**p) |= ( storedValue >> ( number_of_bits - *offset ) ) ;
			}
		}
	}

	std::string construct_snp_probability_data_block_v12(
		uint32_t number_of_samples,
		uint16_t number_of_alleles,
		std::size_t const bits_per_probability,
		std::function< double ( std::size_t i, std::size_t g ) > get_probs,
		std::string const& type
	) {
		uint8_t buffer[10000] ;
		uint8_t* p = buffer ;
		uint8_t* end = buffer + 10000 ;
		std::ostringstream oStream ;
		// Write number of samples and ploidies.
		p = genfile::bgen::write_little_endian_integer( p, end, number_of_samples ) ;
		p = genfile::bgen::write_little_endian_integer( p, end, number_of_alleles ) ;
		uint8_t const ploidy = (number_of_samples > 0) ? 2 : 0 ;
		p = genfile::bgen::write_little_endian_integer( p, end, ploidy ) ;
		p = genfile::bgen::write_little_endian_integer( p, end, ploidy ) ;
		for( std::size_t i = 0; i < number_of_samples; ++i ) {
			p = genfile::bgen::write_little_endian_integer( p, end, ploidy ) ;
		}

		uint8_t const phased = ( type == "phased" || number_of_samples == 0 ) ? 1 : 0 ;
		p = genfile::bgen::write_little_endian_integer( p, end, phased ) ;
		p = genfile::bgen::write_little_endian_integer( p, end, uint8_t( bits_per_probability ) ) ;

		uint64_t const two_to_the_bits = ( uint64_t( 1 ) << bits_per_probability ) ;
		double scale = two_to_the_bits - 1 ;
		std::vector< char > probability_data( std::ceil( 2.0 * number_of_samples * bits_per_probability / 64.0 ) * 8, 0 ) ;
		uint64_t* const probBuffer = reinterpret_cast< uint64_t* >( &probability_data[0] ) ;
		uint64_t* const probEnd = reinterpret_cast< uint64_t* const >( &probability_data[0] + probability_data.size() ) ;
		uint64_t* probP = probBuffer ;
		std::size_t offset = 0 ;
		if( type == "unphased" ) {
			// Construct and write probability data.
			{
				double probs[3] ;
				for( std::size_t i = 0; i < number_of_samples; ++i ) {
					probs[0] = get_probs( i, 0 ) ;
					probs[1] = get_probs( i, 1 ) ;
					probs[2] = get_probs( i, 2 ) ;
					
					bool missing = ( probs[0] == 0.0 && probs[1] == 0.0 && probs[2] == 0.0 ) ;
					if( missing ) {
						*(buffer+8+i) |= 0x80 ;
					} else {
						round_probs_to_simplex( &probs[0], 3, bits_per_probability ) ;
					}
#if DEBUG
					std::cerr << format(
						"sample %d of %d, bits_per_probability = %d, two_to_the_bits=%d, scale = %f, AA=%f, AB=%f, sum = %f\n",
						i, number_of_samples, bits_per_probability, two_to_the_bits, scale, probs[0], probs[1], (probs[0]+probs[1])
					) ;
#endif
					write_probs( &probP, &offset, probEnd, probs, 3, bits_per_probability ) ;
				}
			}
		} else if( type == "phased" ) {
			// Construct and write probability data.
			double probs[2] ;
			for( std::size_t i = 0; i < number_of_samples; ++i ) {
				for( std::size_t hap = 0; hap < 2; ++hap ) {
					probs[0] = get_probs( i, 0+(2*hap) ) ;
					probs[1] = get_probs( i, 1+(2*hap) ) ;
					round_probs_to_simplex( &probs[0], 2, bits_per_probability ) ;
#if DEBUG
					std::cerr << format(
						"sample %d of %d, hap %d, bits_per_probability = %d, two_to_the_bits=%d, scale = %f, AA=%f, AB=%f, sum = %f\n",
						i, number_of_samples, hap, bits_per_probability, two_to_the_bits, scale, probs[0], probs[1], (probs[0]+probs[1])
					) ;
#endif
					write_probs( &probP, &offset, probEnd, &probs[0], 2, bits_per_probability ) ;
				}
			}
		} else {
			assert(0) ;
		}
		std::size_t const numBytes = ((probP-probBuffer)*8) + ((offset+7)/8) ;
		uint32_t const expected_size = (((number_of_samples*2*bits_per_probability)+7)/8) ;
		REQUIRE( numBytes == expected_size ) ;
		p = std::copy( &probability_data[0], &probability_data[0] + numBytes, p ) ;
		return std::string( buffer, p ) ;
	}

	std::string construct_variant_id_data_block_v12(
		uint32_t number_of_samples,
		std::string SNPID,
		std::string RSID,
		std::string chromosome,
		uint32_t SNP_position,
		std::vector< std::string > const alleles
	) {
		std::ostringstream oStream ;
		genfile::bgen::write_little_endian_integer( oStream, static_cast< uint16_t >( SNPID.size() )) ;
		oStream.write( SNPID.data(), SNPID.size() ) ;
		genfile::bgen::write_little_endian_integer( oStream, static_cast< uint16_t >( RSID.size() )) ;
		oStream.write( RSID.data(), RSID.size() ) ;
		genfile::bgen::write_little_endian_integer( oStream, static_cast< uint16_t >( chromosome.size() ) ) ;
		oStream.write( chromosome.data(), chromosome.size() ) ;
		genfile::bgen::write_little_endian_integer( oStream, SNP_position ) ;

		uint16_t const number_of_alleles = alleles.size() ;
		genfile::bgen::write_little_endian_integer( oStream, number_of_alleles ) ;
		for( std::size_t i = 0; i < alleles.size(); ++i ) {
			genfile::bgen::write_little_endian_integer( oStream, static_cast< uint32_t >( alleles[i].size() )) ;
			oStream.write( alleles[i].data(), alleles[i].size() ) ;
		}

		return oStream.str() ;
	}
	
	std::string construct_snp_block_v12(
		uint32_t number_of_samples,
		std::string SNPID,
		std::string RSID,
		std::string chromosome,
		uint32_t SNP_position,
		std::string a_allele,
		std::string b_allele,
		std::size_t const bits_per_probability,
		std::function< double ( std::size_t i, std::size_t g ) > get_probs,
		std::string const& type
	) {
		assert( bits_per_probability <= 64 ) ;
		std::ostringstream oStream ;

		std::vector< std::string > alleles(2) ;
		alleles[0] = a_allele ;
		alleles[1] = b_allele ;

		std::string const id_data = construct_variant_id_data_block_v12(
			number_of_samples, SNPID, RSID, chromosome, SNP_position, alleles
		) ;
		oStream.write( id_data.data(), id_data.size() ) ;

		// Total size of probability data is:
		uint32_t const stored_probability_size = (((number_of_samples*2*bits_per_probability)+7)/8) ;
		// To this we must add space for the header block
		uint32_t const buffer_size = 10 + number_of_samples + stored_probability_size ;
		// And four bytes indicating the total length.
		genfile::bgen::write_little_endian_integer( oStream, buffer_size ) ;

		std::string const prob_data = construct_snp_probability_data_block_v12(
			number_of_samples, 2,
			bits_per_probability,
			get_probs,
			type
		) ;
		oStream.write( prob_data.data(), prob_data.size() ) ;
		return oStream.str() ;
	}
	
	std::string construct_snp_block(
		std::string const& version,
		uint32_t number_of_samples,
		std::string SNPID,
		std::string RSID,
		std::string chromosome,
		uint32_t SNP_position,
		std::string a_allele,
		std::string b_allele,
		std::size_t const bits_per_probability,
		std::function< double ( std::size_t i, std::size_t g ) > get_probs,
		std::string const& type
	) {
#if DEBUG
		std::cerr << "construct_snp_block(): version=" << version << ".\n" ;
#endif
		if( version != "v12" ) {
			assert( bits_per_probability == 16 ) ;
		}

		if( version == "v10" ) {
			return construct_snp_block_v10(
				number_of_samples,
				std::max( SNPID.size(), RSID.size() ) + 1, SNPID, RSID,
				chromosome, SNP_position, a_allele, b_allele,
				get_probs
			) ;
		} else if( version == "v11" ) {
			return construct_snp_block_v11(
				number_of_samples,
				SNPID, RSID, chromosome, SNP_position, a_allele, b_allele,
				get_probs
			) ;
		} else if( version == "v12" ) {
			return construct_snp_block_v12(
				number_of_samples,
				SNPID, RSID, chromosome, SNP_position, a_allele, b_allele,
				bits_per_probability, get_probs, type
			) ;
		} else {
			assert(0) ;
		}
	}
}

namespace {
	enum { e_Layout0 = 0, e_Layout1 = 0x4, e_Layout2 = 0x8 } ;
	enum { e_CompressedSNPBlocks = 1 } ;
}

// The following section defines the needed objects for use with the bgen.hpp implementation.
template< typename T >
struct Setter
{
	Setter( T& field ): m_field( field ) {} ;
	template< typename T2 >
	void operator()( T2 const& value ) { m_field = T(value) ; }
private:
	T& m_field ;
} ;

template< typename T >
Setter< T > make_setter( T& field ) { return Setter<T>( field ) ; }


struct probabilities {
	double AA, AB, BB ;
} ;


double get_input_probability(
	std::size_t const number_of_samples,
	std::size_t i,
	std::size_t g,
	std::string const& type = "unphased"
) {
	double x = 0 ;
	if( type == "phased" ) {
		assert( g < 4 ) ;
		assert( number_of_samples > 0 ) ;
		// two haplotypes, each of whose probs sum to one
		if( g == 0 ) {
			x = double(i) / double(number_of_samples-1) ;
		}
		else if( g == 1 ) {
			x = 1.0 - (double(i) / double(number_of_samples-1)) ;
		} else if( g == 2 ) {
			x = double((i+1) % number_of_samples) / double(number_of_samples-1) ;
		} else if( g == 3 ) {
			x = 1.0 - double((i+1) % number_of_samples) / double(number_of_samples-1) ;
		}
	} else {
		assert( g < 3 ) ;
		if( i % 13 == 0 ) {
			// Every 13th sample is missing
			x = 0.0 ;
		}
		else {
			x = double(i) / double(number_of_samples-1) ;
			if( g == 1 ) {
				// values of second prob are 1 minus 1/3rd of first prob.
				x = 0.25 * ( 1.0 - x ) ;
			} else if( g == 2 ) {
				// values of second prob are 1 minus 2/3rd of first prob.
				x = 0.75 * ( 1.0 - x ) ;
			}
		}
	}
	return x ;
}

double get_expected_stored_probability(
	std::size_t const number_of_samples,
	std::size_t i,
	std::size_t g,
	std::string const& bgen_version,
	std::size_t bits_per_probability,
	std::string const& type = "unphased"
) {
	if( bgen_version == "v10" ) {
		return std::floor( 0.5 + get_input_probability( number_of_samples, i, g ) * 10000.0 ) / 10000.0 ;
	} else if( bgen_version == "v11" ) {
		return std::floor( 0.5 + get_input_probability( number_of_samples, i, g ) * 32768.0 ) / 32768.0 ;
	} else if( bgen_version == "v12" ){
		double v[4] ;
		if( type == "phased" ) {
			for( std::size_t l = 0; l < 4; ++l ) {
				v[l] = get_input_probability( number_of_samples, i, l, type ) ;
			}
			data::round_probs_to_simplex( &v[0]+(2*(g/2)), 2, bits_per_probability ) ;
#if DEBUG > 1
			std::cerr << format( "get_expected_stored_probability(): expected probs are: %f, %f\n", v[(2*(g/2))], v[1+(2*(g/2))] ) ;
#endif
			return *(v+g) ;
		} else {
			double v[3] ;
			for( std::size_t l = 0; l < 3; ++l ) {
				v[l] = get_input_probability( number_of_samples, i, l, type ) ;
			}
			data::round_probs_to_simplex( &v[0], 3, bits_per_probability ) ;
#if DEBUG > 1
			std::cerr << format( "get_expected_stored_probability(): expected probs are: %f, %f, %f\n", v[0], v[1], v[2] ) ;
#endif
			return *(v+g) ;
		}
	} else {
		assert(0) ;
	}
}

// The following section contains the main tests.
void do_snp_block_read_test( 
		std::string const& bgen_version,
		uint32_t number_of_individuals,
		std::string SNPID,
		std::string RSID,
		std::string chromosome,
		uint32_t SNP_position,
		std::string a,
		std::string b,
		std::size_t bits_per_probability = 16,
		std::string const& type = "unphased"
) {
	genfile::bgen::Context context ;
	context.number_of_samples = number_of_individuals ;
	std::function< double ( std::size_t i, std::size_t g ) > get_probs ;
	if( bgen_version == "v11" ) {
		context.flags = e_Layout1 ;
	} else if( bgen_version == "v12" ) {
		context.flags = e_Layout2 ;
	} else {
		context.flags = e_Layout0 ;
	}

	std::istringstream inStream ;
	inStream.str(
		data::construct_snp_block(
			bgen_version, number_of_individuals,
			SNPID, RSID, chromosome, SNP_position, a, b,
			bits_per_probability,
			std::bind( &get_input_probability, number_of_individuals, _1, _2, type ),
			type
		)
	) ;
#if DEBUG
	std::cerr << "do_snp_block_read_test(): bgen_version=" << bgen_version << ".\n" ;
	std::cerr << "do_snp_block_read_test():        data is: " << to_hex( inStream.str() ) << "\n" ;
#endif
	
	std::string SNPID2 ;
	std::string RSID2 ;
	std::string chromosome2 ;
	uint32_t SNP_position2 ;
	std::string a2 ;
	std::string b2 ;

	std::vector< genfile::byte_t > buffer1, buffer2 ;

	genfile::bgen::read_snp_identifying_data(
		inStream,
		context,
		&SNPID2,
		&RSID2,
		&chromosome2,
		&SNP_position2,
		&a2,
		&b2
	) ;

	REQUIRE( SNPID2 == SNPID ) ;
	REQUIRE( RSID2 == RSID ) ;
	bool const chromosomeCheck = ( chromosome2 == chromosome || chromosome2 == "0" + chromosome ) ;
	REQUIRE( chromosomeCheck ) ;
	REQUIRE( SNP_position2 == SNP_position ) ;
	REQUIRE( a2 == a ) ;
	REQUIRE( b2 == b ) ;

	ProbSetCheck setter(
		number_of_individuals,
		std::bind(
			&get_expected_stored_probability,
			number_of_individuals,
			_1,
			_2,
			bgen_version,
			bits_per_probability,
			type
		)
	) ;

	genfile::bgen::read_and_parse_genotype_data_block(
		inStream,
		context,
		setter,
		&buffer1,
		&buffer2
	) ;
}

void do_snp_block_write_test(
		std::string const& bgen_version,
		uint32_t number_of_individuals,
		std::string SNPID,
		std::string RSID,
		std::string chromosome,
		uint32_t SNP_position,
		std::string a,
		std::string b,
		std::size_t bits_per_probability = 16,
		std::string const& type = "unphased"
) {
	using namespace std::placeholders ;
	genfile::bgen::Context context ;
	context.number_of_samples = number_of_individuals ;
	if( bgen_version == "v11" ) {
		context.flags = e_Layout1 ;
	} else if( bgen_version == "v12" ) {
		context.flags = e_Layout2 ;
	} else {
		assert(0) ;	
	}
	
#if DEBUG
	std::cerr
		<< "do_snp_block_write_test(): bgen_version=" << bgen_version
		<< ", number_of_samples = " << number_of_individuals
		<< ", number_of_bits = " << bits_per_probability
		<< ".\n" ;
#endif
	std::ostringstream outStream ;
	
	std::vector< genfile::byte_t > buffer ;
	std::vector< genfile::byte_t > buffer2 ;
	genfile::byte_t* end = genfile::bgen::write_snp_identifying_data(
		&buffer,
		context,
		SNPID,
		RSID,
		chromosome,
		SNP_position,
		2, [&a,&b]( std::size_t i ) { if( i == 0 ) { return a ; } else if( i == 1 ) { return b ; } else { assert(0) ; }}
	) ;
	outStream.write( reinterpret_cast< char* >( &buffer[0] ), end - &buffer[0] ) ;

	genfile::bgen::GenotypeDataBlockWriter writer( &buffer, &buffer2, context, bits_per_probability ) ;
	writer.initialise( number_of_individuals, 2 ) ;
	for( std::size_t i = 0; i < number_of_individuals; ++i ) {
		writer.set_sample( i ) ;
		if( type == "unphased" ) {
			writer.set_number_of_entries( 2, 3, genfile::ePerUnorderedGenotype, genfile::eProbability ) ;
			writer.set_value( 0, get_input_probability( number_of_individuals, i, 0, type )) ;
			writer.set_value( 1, get_input_probability( number_of_individuals, i, 1, type )) ;
			writer.set_value( 2, get_input_probability( number_of_individuals, i, 2, type )) ;
		} else {
			writer.set_number_of_entries( 2, 4, genfile::ePerPhasedHaplotypePerAllele, genfile::eProbability ) ;
			writer.set_value( 0, get_input_probability( number_of_individuals, i, 0, type )) ;
			writer.set_value( 1, get_input_probability( number_of_individuals, i, 1, type )) ;
			writer.set_value( 2, get_input_probability( number_of_individuals, i, 2, type )) ;
			writer.set_value( 3, get_input_probability( number_of_individuals, i, 3, type )) ;
		}
	}
	writer.finalise() ;
	outStream.write( reinterpret_cast< char const* >( writer.repr().first ), writer.repr().second - writer.repr().first ) ;

	std::string expected = data::construct_snp_block(
		bgen_version,
		number_of_individuals,
		SNPID,
		RSID,
		chromosome,
		SNP_position,
		a,
		b,
		bits_per_probability,
		std::bind( &get_input_probability, number_of_individuals, _1, _2, type ),
		type
	) ;

#if DEBUG
	std::cerr << "          actual: " << to_hex( outStream.str() ) << "\n" ;
	std::cerr << "        expected: " << to_hex( expected ) << "\n" ;
#endif

	REQUIRE( to_hex( outStream.str() ) == to_hex( expected ) ) ;
	REQUIRE( outStream.str() == expected ) ;
}

TEST_CASE( "Test that compute_approximate_probabilities() works correctly", "[bgen]" ) {
	std::cout << "test_round_probs\n" ;
	double p[10] ;
	std::size_t anIndex[10] ;

	p[0] = 1; p[1] = 0; p[2] = 0 ;
	genfile::bgen::v12::impl::compute_approximate_probabilities( p, anIndex, 3, 1 ) ;
	REQUIRE( p[0] == 1 ) ;
	REQUIRE( p[1] == 0 ) ;
	REQUIRE( p[2] == 0 ) ;

	p[0] = 0; p[1] = 1; p[2] = 0 ;
	genfile::bgen::v12::impl::compute_approximate_probabilities( p, anIndex, 3, 1 ) ;
	REQUIRE( p[0] == 0 ) ;
	REQUIRE( p[1] == 1 ) ;
	REQUIRE( p[2] == 0 ) ;

	p[0] = 0; p[1] = 0; p[2] = 1 ;
	genfile::bgen::v12::impl::compute_approximate_probabilities( p, anIndex, 3, 1 ) ;
	REQUIRE( p[0] == 0 ) ;
	REQUIRE( p[1] == 0 ) ;
	REQUIRE( p[2] == 1 ) ;
	
	p[0] = 1; p[1] = 0; p[2] = 0 ;
	genfile::bgen::v12::impl::compute_approximate_probabilities( p, anIndex, 3, 8 ) ;
	REQUIRE( p[0] == 255 ) ;
	REQUIRE( p[1] == 0 ) ;
	REQUIRE( p[2] == 0 ) ;

	p[0] = 0; p[1] = 1; p[2] = 0 ;
	genfile::bgen::v12::impl::compute_approximate_probabilities( p, anIndex, 3, 8 ) ;
	REQUIRE( p[0] == 0 ) ;
	REQUIRE( p[1] == 255 ) ;
	REQUIRE( p[2] == 0 ) ;

	p[0] = 0; p[1] = 0; p[2] = 1 ;
	genfile::bgen::v12::impl::compute_approximate_probabilities( p, anIndex, 3, 8 ) ;
	REQUIRE( p[0] == 0 ) ;
	REQUIRE( p[1] == 0 ) ;
	REQUIRE( p[2] == 255 ) ;

	p[0] = 0.3; p[1] = 0.3; p[2] = 0.4 ;
	genfile::bgen::v12::impl::compute_approximate_probabilities( p, anIndex, 3, 1 ) ;
	REQUIRE( p[0] == 0 ) ;
	REQUIRE( p[1] == 0 ) ;
	REQUIRE( p[2] == 1 ) ;

	p[0] = 0.3; p[1] = 0.3; p[2] = 0.4 ;
	genfile::bgen::v12::impl::compute_approximate_probabilities( p, anIndex, 3, 2 ) ;
	REQUIRE( p[0] == 1 ) ;
	REQUIRE( p[1] == 1 ) ;
	REQUIRE( p[2] == 1 ) ;

	// up, down, exact.
	p[0] = 0.3; p[1] = 0.3; p[2] = 0.4 ;
	genfile::bgen::v12::impl::compute_approximate_probabilities( p, anIndex, 3, 4 ) ;
	REQUIRE( p[0] == 5 ) ;
	REQUIRE( p[1] == 4 ) ;
	REQUIRE( p[2] == 6 ) ;

	// up, down, exact.
	p[0] = 76.5/255.0; p[1] = 76.5/255.0; p[2] = 102.0/255.0 ;
	genfile::bgen::v12::impl::compute_approximate_probabilities( p, anIndex, 3, 8 ) ;
	REQUIRE( p[0] == 77 ) ;
	REQUIRE( p[1] == 76 ) ;
	REQUIRE( p[2] == 102 ) ;

	// down, up, down.
	p[0] = 76.2/255.0; p[1] = 76.9/255.0; p[2] = 101.9/255.0 ;
	genfile::bgen::v12::impl::compute_approximate_probabilities( p, anIndex, 3, 8 ) ;
	REQUIRE( p[0] == 76 ) ;
	REQUIRE( p[1] == 77 ) ;
	REQUIRE( p[2] == 102 ) ;
}

TEST_CASE( "Test that valid v1.0 variant data block can be read", "[bgen][v10]" ) {
	std::cout << "test_snp_block_input_v10\n" ;
//	do_snp_block_read_test( "v10", 0, "SNP01", "RS01", "1", 1000001, "A", "C" ) ;
	do_snp_block_read_test( "v10", 6, "SNP01", "RS01", "1", 1000001, "A", "C" ) ;
	do_snp_block_read_test( "v10", 15, "01234567890123456789012345678901234567890123456789", "01234567890123456789012345678901234567890123456789", "22", 4294967295u, "G", "T" ) ;
	do_snp_block_read_test( "v10", 100, "01234567890123456789012345678901234567890123456789", "01234567890123456789012345678901234567890123456789", "22", 4294967295u, "G", "T" ) ;
	do_snp_block_read_test( "v10", 1001, "01234567890123456789012345678901234567890123456789", "01234567890123456789012345678901234567890123456789", "22", 4294967295u, "G", "T" ) ;
}

TEST_CASE( "Test that valid v1.1 variant data block can be read", "[bgen][v11]" ) {
	std::cout << "test_snp_block_input_v11\n" ;
//	do_snp_block_read_test( "v11", 0, "SNP01", "RS01", "1", 1000001, "A", "C" ) ;
	do_snp_block_read_test( "v11", 6, "01234567890123456789012345678901234567890123456789", "01234567890123456789012345678901234567890123456789", "22", 4294967295u, "G", "T" ) ;
	do_snp_block_read_test( "v11", 15, "SNP01", "RS01", "1", 1000001, "A", "C" ) ;
	do_snp_block_read_test( "v11", 100, "01234567890123456789012345678901234567890123456789", "01234567890123456789012345678901234567890123456789", "22", 4294967295u, "G", "T" ) ;
	do_snp_block_read_test( "v11", 1001, "SNP01", "RS01", "1", 1000001, "A", "C" ) ;
}

TEST_CASE( "Test that valid v1.1 variant data block can be written", "[bgen][v11]" ) {
	std::cout << "test_snp_block_output_unphased_v11\n" ;
	do_snp_block_write_test( "v11", 6, "SNP01", "RS01", "1", 1000001, "A", "C" ) ;
	do_snp_block_write_test( "v11", 100, "01234567890123456789012345678901234567890123456789", "01234567890123456789012345678901234567890123456789", "22", 4294967295u, "G", "T" ) ;
}

TEST_CASE( "Test that valid v1.2 variant data block containing unphased data can be read", "[bgen][v12]" ) {
	std::cout << "test_snp_block_input_unphased_v12\n" ;
	for( std::size_t number_of_bits = 1; number_of_bits <= 32; ++number_of_bits ) {
		do_snp_block_read_test( "v12", 0, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits, "unphased" ) ;
		do_snp_block_read_test( "v12", 2, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits, "unphased" ) ;
		do_snp_block_read_test( "v12", 6, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits, "unphased" ) ;
		do_snp_block_read_test( "v12", 15, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits, "unphased" ) ;
		do_snp_block_read_test( "v12", 37, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits, "unphased" ) ;
		do_snp_block_read_test( "v12", 100, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits, "unphased" ) ;
		do_snp_block_read_test( "v12", 1001, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits, "unphased" ) ;
	}
}

TEST_CASE( "Test that valid v1.2 variant data block containing unphased data can be written", "[bgen][v12]" ) {
	std::cout << "test_snp_block_output_unphased_v12\n" ;
	do_snp_block_write_test( "v11", 6, "SNP01", "RS01", "1", 1000001, "A", "C" ) ;
	do_snp_block_write_test( "v11", 100, "01234567890123456789012345678901234567890123456789", "01234567890123456789012345678901234567890123456789", "22", 4294967295u, "G", "T" ) ;
	for( std::size_t number_of_bits = 1; number_of_bits <= 32; ++number_of_bits ) {
		do_snp_block_write_test( "v12", 0, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits ) ;
		do_snp_block_write_test( "v12", 2, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits ) ;
		do_snp_block_write_test( "v12", 6, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits ) ;
		do_snp_block_write_test( "v12", 15, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits ) ;
		do_snp_block_write_test( "v12", 37, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits ) ;
		do_snp_block_write_test( "v12", 100, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits ) ;
		do_snp_block_write_test( "v12", 1001, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits ) ;
	}
}

TEST_CASE( "Test that valid v1.2 variant data block containing phased data can be read", "[bgen][v12]" ) {
	std::cout << "test_snp_block_input_phased_v12\n" ;
	for( std::size_t number_of_bits = 1; number_of_bits <= 32; ++number_of_bits ) {
		do_snp_block_read_test( "v12", 0, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_read_test( "v12", 2, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_read_test( "v12", 6, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_read_test( "v12", 15, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_read_test( "v12", 37, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_read_test( "v12", 100, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_read_test( "v12", 1001, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits, "phased" ) ;
	}
}

TEST_CASE( "Test that valid variant data block containing phased data can be written", "[bgen][v12]" ) {
	std::cout << "test_snp_block_output_phased_v12\n" ;
	for( std::size_t number_of_bits = 1; number_of_bits <= 32; ++number_of_bits ) {
		do_snp_block_write_test( "v12", 0, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_write_test( "v12", 2, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_write_test( "v12", 6, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_write_test( "v12", 15, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_write_test( "v12", 37, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_write_test( "v12", 100, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_write_test( "v12", 1001, "SNP01", "RS01", "1", 1000001, "A", "C", number_of_bits, "phased" ) ;
	}
}
