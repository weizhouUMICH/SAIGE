
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

TEST_CASE( "Test probability bound", "[bgen][biallelic][unphased]" ) {
	std::cerr << "test_probability_bound\n" ;
	std::vector< genfile::byte_t > buffer( 1000 ) ;
	double epsilon = std::numeric_limits< double >::epsilon() ; // Something small and nonzero.

	std::vector< double > errors ;
	errors.push_back( 1.0/2048.0 ) ;
	errors.push_back( 1.0/1024.0 ) ;
	errors.push_back( 0.0005 ) ;
	errors.push_back( 0.005 ) ;
	errors.push_back( 0.02 ) ;
	errors.push_back( 0.05 ) ;
	for( std::size_t i = 0; i < errors.size(); ++i ) {
		double const& error = errors[i] ;
		for( std::size_t g = 0; g < 3; ++g ) {
			{
				genfile::bgen::v12::ProbabilityDataWriter writer( 16, error ) ;
				writer.initialise( 1, 2, &buffer[0], &buffer[0] + buffer.size() ) ;
				writer.set_sample( 0 ) ;
				writer.set_number_of_entries( 3, 3, genfile::ePerUnorderedGenotype, genfile::eProbability ) ;
				for( std::size_t k = 0; k < 3; ++k ) {
					if( k == g ) {
						REQUIRE_NOTHROW( writer.set_value( k, 1.0 + error ) ) ;
					} else {
						REQUIRE_NOTHROW( writer.set_value( k, 0 ) ) ;
					}
				}
				REQUIRE_NOTHROW( writer.finalise() ) ;
			}

			{
				genfile::bgen::v12::ProbabilityDataWriter writer( 16, error ) ;
				writer.initialise( 1, 2, &buffer[0], &buffer[0] + buffer.size() ) ;
				writer.set_sample( 0 ) ;
				writer.set_number_of_entries( 3, 3, genfile::ePerUnorderedGenotype, genfile::eProbability ) ;
				for( std::size_t k = 0; k < 3; ++k ) {
					if( k == g ) {
						REQUIRE_THROWS_AS( writer.set_value( k, 1.0 + error + 2 * epsilon ), genfile::bgen::BGenError ) ;
						break ;
					} else {
						REQUIRE_NOTHROW( writer.set_value( k, 0 ) ) ;
					}
				}
			}

			{
				genfile::bgen::v12::ProbabilityDataWriter writer( 16, error ) ;
				writer.initialise( 1, 2, &buffer[0], &buffer[0] + buffer.size() ) ;
				writer.set_sample( 0 ) ;
				writer.set_number_of_entries( 3, 3, genfile::ePerUnorderedGenotype, genfile::eProbability ) ;
				for( std::size_t k = 0; k < 3; ++k ) {
					double v = ((k==g) ? 0.5 : 0.25) + error ;
					REQUIRE_NOTHROW( writer.set_value( k, v ) ) ;
				}
				REQUIRE_NOTHROW( writer.finalise() ) ;
			}

			{
				genfile::bgen::v12::ProbabilityDataWriter writer( 16, error ) ;
				writer.initialise( 1, 2, &buffer[0], &buffer[0] + buffer.size() ) ;
				writer.set_sample( 0 ) ;
				writer.set_number_of_entries( 3, 3, genfile::ePerUnorderedGenotype, genfile::eProbability ) ;
				for( std::size_t k = 0; k < 3; ++k ) {
					double v = ((k==g) ? 0.5 : 0.25) - error ;
					REQUIRE_NOTHROW( writer.set_value( k, v ) ) ;
				}
				REQUIRE_NOTHROW( writer.finalise() ) ;
			}

			{
				genfile::bgen::v12::ProbabilityDataWriter writer( 16, error ) ;
				writer.initialise( 1, 2, &buffer[0], &buffer[0] + buffer.size() ) ;
				writer.set_sample( 0 ) ;
				writer.set_number_of_entries( 3, 3, genfile::ePerUnorderedGenotype, genfile::eProbability ) ;
				for( std::size_t k = 0; k < 3; ++k ) {
					double v = ((k==g) ? (0.5+6*epsilon) : 0.25) + error;
					if( k < 2 ) {
						REQUIRE_NOTHROW( writer.set_value( k, v ) ) ;
					} else {
						// last value will throw.
						REQUIRE_THROWS_AS(  writer.set_value( k, v ), genfile::bgen::BGenError ) ;
					}
				}
			}

			{
				genfile::bgen::v12::ProbabilityDataWriter writer( 16, error ) ;
				writer.initialise( 1, 2, &buffer[0], &buffer[0] + buffer.size() ) ;
				writer.set_sample( 0 ) ;
				writer.set_number_of_entries( 3, 3, genfile::ePerUnorderedGenotype, genfile::eProbability ) ;
				for( std::size_t k = 0; k < 3; ++k ) {
					double v = ((k==g) ? 0.5-6*epsilon : 0.25) - error ;
					if( k < 2 ) {
						REQUIRE_NOTHROW( writer.set_value( k, v ) ) ;
					} else {
						// last value will throw
						REQUIRE_THROWS_AS( writer.set_value( k, v ), genfile::bgen::BGenError ) ;
					}
				}
			}
		}
	}
}

TEST_CASE( "Single sample (biallelic unphased)", "[bgen][biallelic][unphased]" ) {
	std::cout << "test_single_sample_biallelic_unphased\n" ;

	std::vector< genfile::byte_t > data( 1000 ) ;

	genfile::bgen::Context context ;
	context.flags = genfile::bgen::e_Layout2 ;
	context.number_of_samples = 1 ;

	enum { ePloidy = 8, eNumberOfBits = 10, eData = 11 } ;

	std::vector< genfile::byte_t > expected{
		0x1,0x0,0x0,0x0,		// sample count
		0x2,0x0,				// allele count
		0x2,0x2,				// min/max ploidy
		0x82, 					// ploidy
		0x0, 0x0,				// phased, # bits to be filled in below
		// For data we need a maximum of 63 (max ploidy) x 2 (number of alleles) x 4 bytes = 504 bytes.
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0
	} ;
	
	std::pair< genfile::byte_t const*, genfile::byte_t const* > repr ;
	std::size_t const number_of_alleles = 2 ;
	for( uint8_t ploidy = 0; ploidy < 63; ++ploidy ) {
		expected[6] = expected[7] = ploidy ;
		expected[8] = ploidy ;
		for( unsigned char number_of_bits = 1; number_of_bits <= 32; ++number_of_bits ) {
			expected[ eNumberOfBits ] = number_of_bits ;
			// Always treat 0- or haploid as phased.
			uint8_t const phased = ( ploidy > 1 ) ? 0 : 1 ;
			expected[9] = phased ;
			// For a biallelic (K=2) variant the number of probability values is ploidy+1 if unphased
			// or ploidy*2 if phased.
			std::size_t const numberOfValues = (ploidy+1) ;
			// We also put a 1.0 prob in every possible position
			std::size_t const numberOfDistinctGenotypes = numberOfValues ;
			std::size_t const numberOfInferredValues = 1 ;
			// Expected data block size for one sample:
			std::size_t expected_size = 11 + (( number_of_bits * (numberOfValues - numberOfInferredValues)) + 7 ) / 8 ;

			for( std::size_t g = 0; g <= numberOfDistinctGenotypes; ++g ) { // g=numberOfValues => missing.
#if DEBUG
				std::cerr << "test_single_sample(): ploidy=" << int( ploidy ) << ", phased=" << int( phased ) << ", bits=" << int( number_of_bits ) << ", g=" << g << ", numberOfValues=" << numberOfValues << ".\n" ;
#endif				

				// Now we need to set the expected genotype.
				// We do it with a trick as follows.
				expected[8] = ploidy ;
				std::fill( &expected[0] + eData, &expected[0] + expected.size(), uint8_t(0) ) ;
				if( g < (numberOfValues-1) ) {
					std::size_t bit = g * number_of_bits ;
					std::size_t byte = bit / 8 ;
					bit = bit % 8 ;
					// Create a bitmask with the right number of bits.
					uint64_t mask = ( uint64_t( 0xFFFFFFFF ) >> ( 32 - number_of_bits )) << bit ;

#if DEBUG
					std::cerr << "test_single_sample(): setting byte " << eData + byte << " to "  << std::hex << mask << std::dec << ".\n" ;
#endif
					for( std::size_t i = 0; i < number_of_bits*ploidy; i += 8 ) {
						expected[ eData + byte + (i/8) ] = mask & 0xFF ;
						mask >>= 8 ;
					}
				} else if( g == numberOfDistinctGenotypes ) {
					expected[8] |= 0x80 ;
				}
				{
					genfile::bgen::v12::ProbabilityDataWriter writer( number_of_bits ) ;
					writer.initialise( 1, number_of_alleles, &data[0], &data[0] + data.size() ) ;
					writer.set_sample( 0 ) ;
					writer.set_number_of_entries(
						ploidy, numberOfValues,
						( phased ? genfile::ePerPhasedHaplotypePerAllele : genfile::ePerUnorderedGenotype ),
						genfile::eProbability
					) ;
					for( std::size_t i = 0; i < numberOfValues; ++i ) {
						if( g == numberOfDistinctGenotypes ) {
							writer.set_value( i, genfile::MissingValue() ) ;
						} else if( i == g ) {
							writer.set_value( i, 1.0 ) ;
						} else {
							writer.set_value( i, 0.0 ) ;
						}
					}
					writer.finalise() ;
					repr = writer.repr() ;
				}
	#if DEBUG
				std::cerr << "test_single_sample(): (bits=" << int( number_of_bits ) << ", ploidy=" << int(ploidy) <<", g=" << g << "):   result is:" << to_hex( repr.first, repr.second ) << ".\n" ;
				std::cerr << "test_single_sample(): (bits=" << int( number_of_bits ) << ", ploidy=" << int(ploidy) <<", g=" << g << "): expected is:" << to_hex( &expected[0], &expected[0] + expected_size ) << ".\n" ;
	#endif
				REQUIRE( ( repr.second - repr.first ) == expected_size ) ;
				REQUIRE( std::memcmp( repr.first, expected.data(), expected_size ) == 0 ) ;
			}
		}
	}
}

TEST_CASE( "Single sample (biallelic phased)", "[bgen][biallelic][phased]" ) {
	std::cout << "test_single_sample_biallelic_phased\n" ;

	std::vector< genfile::byte_t > data( 1000 ) ;

	genfile::bgen::Context context ;
	context.flags = genfile::bgen::e_Layout2 ;
	context.number_of_samples = 1 ;

	enum { ePloidy = 8, eNumberOfBits = 10, eData = 11 } ;

	std::vector< genfile::byte_t > expected{
		0x1,0x0,0x0,0x0,		// sample count
		0x2,0x0,				// allele count
		0x2,0x2,				// min/max ploidy
		0x82, 					// ploidy
		0x0, 0x0,				// phased, # bits to be filled in below
		// For data we need a maximum of 63 (max ploidy) x 2 (number of alleles) x 4 bytes = 504 bytes.
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0
	} ;
	
	std::pair< genfile::byte_t const*, genfile::byte_t const* > repr ;
	std::size_t const number_of_alleles = 2 ;
	for( uint8_t ploidy = 0; ploidy < 63; ++ploidy ) {
		expected[6] = expected[7] = ploidy ;
		expected[8] = ploidy ;
		for( unsigned char number_of_bits = 1; number_of_bits <= 32; ++number_of_bits ) {
			expected[ eNumberOfBits ] = number_of_bits ;
			uint8_t const phased = 1 ;
			expected[9] = phased ;
			// For a biallelic (K=2) variant the number of probability values is ploidy+1 if unphased
			// or ploidy*2 if phased.
			std::size_t const numberOfValues = (2*ploidy) ;
			// For low ploidy we test with a 1.0 prob in every possible position of each haplotype.
			// For higher ploidy this is too many values so we just put a 1 in each possible
			// position of the combined genotype.
			bool useSimpleGenotype = (ploidy == 0 || ploidy > 8) ;
			std::size_t const numberOfDistinctGenotypes = useSimpleGenotype ? numberOfValues : ( 1 << ploidy) ;
			// Expected block size for one sample
			std::size_t const numberOfInferredValues = ploidy ;
			std::size_t expected_size = 11 + (( number_of_bits * (numberOfValues - numberOfInferredValues) ) + 7 ) / 8 ;

			for( uint32_t g = 0; g <= numberOfDistinctGenotypes; ++g ) { // g=numberOfValues => missing.
#if DEBUG
				std::cerr << "test_single_sample(): ploidy=" << int( ploidy ) << ", phased=" << int( phased ) << ", bits=" << int( number_of_bits ) << ", g=" << g << ", numberOfValues=" << numberOfValues << ".\n" ;
#endif				

				// Now we need to set the expected genotype.
				// We do it with a trick as follows.
				expected[8] = ploidy ;
				std::fill( &expected[0] + eData, &expected[0] + expected.size(), uint8_t(0) ) ;
				if( g < numberOfDistinctGenotypes ) {
					for( std::size_t elt = 0; elt < ploidy; ++elt ) {
						// Use the bit-pattern of g to determine genotype.
						// Note: bgen stores the first probability (i.e. prob of alleleA per ploid copy)
						// so we store nonzero value iff the corresponding bit is 0.
						bool setBit = ( useSimpleGenotype ? (g == 2*elt) : ((( g >> elt ) & 0x1) == 0 ) ) ;
						if( setBit ) { //
							std::size_t bit = elt * number_of_bits ;
							std::size_t byte = bit / 8 ;
							bit = bit % 8 ;
							// Create a bitmask with the right number of bits.
							uint64_t mask = ( uint64_t( 0xFFFFFFFF ) >> ( 32 - number_of_bits )) << bit ;

	#if DEBUG
							std::cerr << "test_single_sample(): setting byte " << eData + byte << " to "  << std::hex << mask << std::dec << ".\n" ;
	#endif
							for( std::size_t i = 0; i < number_of_bits*ploidy; i += 8 ) {
								expected[ eData + byte + (i/8) ] |= mask & 0xFF ;
								mask >>= 8 ;
							}
						}
					}
				} else if( ploidy > 0 && g == numberOfDistinctGenotypes ) {
					expected[8] |= 0x80 ;
				}

				// Now use bgen writer to set expected genotype
				{
					genfile::bgen::v12::ProbabilityDataWriter writer( number_of_bits ) ;
					writer.initialise( 1, number_of_alleles, &data[0], &data[0] + data.size() ) ;
					writer.set_sample( 0 ) ;
					writer.set_number_of_entries(
						ploidy, numberOfValues,
						genfile::ePerPhasedHaplotypePerAllele,
						genfile::eProbability
					) ;
					for( std::size_t i = 0; i < ploidy; ++i ) {
						if( g == numberOfDistinctGenotypes ) {
							writer.set_value( 2*i+0, genfile::MissingValue() ) ;
							writer.set_value( 2*i+1, genfile::MissingValue() ) ;
						} else {
							bool setFirstAllele = ( useSimpleGenotype ? (g == (2*i)) : ((( g >> i ) & 0x1) == 0 ) ) ;
							writer.set_value( 2*i+0, setFirstAllele ? 1.0 : 0.0 ) ;
							writer.set_value( 2*i+1, setFirstAllele ? 0.0 : 1.0 ) ;
						}
					}
					writer.finalise() ;
					repr = writer.repr() ;
				}
	#if DEBUG
				std::cerr << "test_single_sample(): (bits=" << int( number_of_bits ) << ", ploidy=" << int(ploidy) <<", g=" << g << "):   result is:" << to_hex( repr.first, repr.second ) << ".\n" ;
				std::cerr << "test_single_sample(): (bits=" << int( number_of_bits ) << ", ploidy=" << int(ploidy) <<", g=" << g << "): expected is:" << to_hex( &expected[0], &expected[0] + expected_size ) << ".\n" ;
	#endif
				REQUIRE( ( repr.second - repr.first ) == expected_size ) ;
				REQUIRE( std::memcmp( repr.first, expected.data(), expected_size ) == 0 ) ;
			}
		}
	}
}

TEST_CASE( "Test two samples, unphased", "[bgen][small]" ) {
	std::cout << "test_two_samples\n" ;

	std::vector< genfile::byte_t > data( 1000 ) ;
	std::vector< double > values{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;
	std::function< double( std::size_t ) > get_AA = [&]( std::size_t i ) { return values[i*3] ; } ;
	std::function< double( std::size_t ) > get_AB = [&]( std::size_t i ) { return values[i*3+1] ; } ;
	std::function< double( std::size_t ) > get_BB = [&]( std::size_t i ) { return values[i*3+2] ; } ;

	genfile::bgen::Context context ;
	context.flags = genfile::bgen::e_Layout2 ;
	context.number_of_samples = 2 ;

	std::vector< genfile::byte_t > expected{
		0x2,0x0,0x0,0x0,		// sample count
		0x2,0x0,				// allele count
		0x2,0x2,				// min/max ploidy
		0x82, 0x82,				// ploidy
		0x0, 0x0,				// phased, # bits to be filled in below
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,		// data (padded to maximum size).
		0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0		// ditto.
	} ;
	enum { ePloidy = 8, eNumberOfBits = 11, eData = 12 } ;

	for( std::size_t g = 0; g < 4; ++g ) {
		values[0] = values[1] = values[2] = 0 ;
		values[g] = ( g == 3 ? 0.0 : 1.0 ) ;
		expected[ ePloidy ] = (g == 3) ? 0x82 : 0x2 ;
		expected[ eData ] = expected[ eData + 1 ] = 0x0 ;
		if( g < 3 ) {
			expected[ eData + g ] = 0xFF ;
		}
		for( std::size_t g2 = 0; g2 < 3; ++g2 ) {
			values[3] = values[4] = values[5] = 0 ;
			values[g2+3] = ( g2 == 3 ? 0.0 : 1.0 ) ;
			expected[ ePloidy+1 ] = (g2 == 3) ? 0x82 : 0x2 ;
			expected[ eData + 2 ] = expected[ eData + 3 ] = 0x0 ;
			if( g2 < 2 ) {
				expected[ eData + 2 + g2 ] = 0xFF ;
			}
			unsigned char number_of_bits = 8 ;
			{
				expected[ eNumberOfBits ] = number_of_bits ;
				std::size_t expected_size = 12 + ((( number_of_bits * 4 ) + 7 ) / 8 ) ;

				genfile::bgen::v12::ProbabilityDataWriter writer( number_of_bits ) ;
				writer.initialise( 2, 2, &data[0], &data[0] + data.size() ) ;
				for( std::size_t i = 0; i < 2; ++i ) {
					writer.set_sample( i ) ;
					writer.set_number_of_entries( 2, 3, genfile::ePerUnorderedGenotype, genfile::eProbability ) ;
					writer.set_value( 0, get_AA(i) ) ;
					writer.set_value( 1, get_AB(i) ) ;
					writer.set_value( 2, get_BB(i) ) ;
				}
				writer.finalise() ;

				std::pair< genfile::byte_t const*, genfile::byte_t const* > repr = writer.repr() ;
	#if DEBUG
				std::cerr << "test_two_samples(): (bits=" << int( number_of_bits ) << "):   result is:" << to_hex( repr.first, repr.second ) << ".\n" ;
				std::cerr << "test_two_samples(): (bits=" << int( number_of_bits ) << "): expected is:" << to_hex( &expected[0], &expected[0] + expected_size ) << ".\n" ;
	#endif
				REQUIRE( ( repr.second - repr.first ) == expected_size ) ;
				REQUIRE( std::memcmp( repr.first, expected.data(), expected_size ) == 0 ) ;
			}
		}
	}
}

TEST_CASE( "Test truncated data", "[bgen]" ) {
	std::cout << "test_truncated\n" ;

	std::vector< char > data( 1000 ) ;
	std::vector< double > values{ 0.0, 0.0, 1.0, 1.0, 0.0, 0.0 } ;
	std::function< double( std::size_t ) > get_AA = [&]( std::size_t i ) { return values[i*3] ; } ;
	std::function< double( std::size_t ) > get_AB = [&]( std::size_t i ) { return values[i*3+1] ; } ;
	std::function< double( std::size_t ) > get_BB = [&]( std::size_t i ) { return values[i*3+2] ; } ;

	genfile::bgen::Context context ;
	context.number_of_samples = 2 ;
	context.flags = genfile::bgen::e_Layout2 ;

	enum { ePloidy = 8, eNumberOfBits = 10, eData = 11 } ;

	std::vector< unsigned char > expected{
		0x2,0x0,0x0,0x0,		// sample count
		0x2,0x0,				// allele count
		0x2,0x2,				// min/max ploidy
		0x02, 0x02,				// ploidy
		0x0, 0x8,				// phased, # bits to be filled in below
		0x0, 0x0, 0xFF, 0x0		// data (padded to maximum size).
	} ;

	// Check that truncated data fails...
	for( std::size_t i = 1; i < expected.size(); ++i ) {
#if DEBUG
		std::cerr << "test_malformed(): i = " << i << ": " << to_hex( &expected[0], &expected[0] + expected.size() - i ) << "\n" ;
#endif
		ProbSetCheck setter(
			2,
			[&values]( std::size_t i, std::size_t g ) { return values[ i*3 + g ] ; }
		) ;
		setter.expect_fail() ;
		REQUIRE_THROWS(
			genfile::bgen::parse_probability_data(
				&expected[0],
				&expected[0] + expected.size() - i,
				context,
				setter
			)
		) ;
	}
}

