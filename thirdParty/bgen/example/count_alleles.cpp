
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <stdexcept>
#include <memory>
#include <numeric>
#include "genfile/bgen/bgen.hpp"
#include "genfile/bgen/View.hpp"
#include "genfile/bgen/IndexQuery.hpp"

// AlleleCounter is a callback object appropriate
// for passing to bgen::read_genotype_data_block() or the synonymous method of genfile::bgen::View.

// AlleleCounter is a callback object appropriate for passing to bgen::read_genotype_data_block() or
// the synonymous method of genfile::bgen::View. See the comments below, comments in bgen.hpp,
// or the bgen wiki for a description of the API.
// The purpose of AlleleCounter is to accumulate the probability mass on each
// allele to compute expected allele counts at each variant.
struct AlleleCounter {
	AlleleCounter() {}
	
	// Called once per variant allowing us to set storage.
	void initialise( std::size_t number_of_samples, std::size_t number_of_alleles ) {
		m_number_of_alleles = number_of_alleles ;
		m_expected_allele_counts.assign( number_of_alleles, 0.0 ) ;
	}

	// If present with this signature, called once after initialise()
	// to set the minimum and maximum ploidy and number of probabilities among samples in the data.
	// This enables us to set up storage for the data ahead of time.
	void set_min_max_ploidy( uint32_t min_ploidy, uint32_t max_ploidy, uint32_t min_entries, uint32_t max_entries ) {
		// Make sure we've enough space to store probs
		m_data.reserve( max_entries ) ;
	}

	// Called once per sample to determine whether we want data for this sample
	bool set_sample( std::size_t i ) {
		// Yes, here we want info for all samples.
		return true ;
	}

	// Called once per sample to set the number of probabilities
	// that are present for this sample, as well as whether the data is phased
	// or unphased, etc.
	void set_number_of_entries(
		std::size_t ploidy,
		std::size_t number_of_entries,
		genfile::OrderType order_type,
		genfile::ValueType value_type
	) {
		assert( value_type == genfile::eProbability ) ;
		m_data.resize( number_of_entries ) ;
		m_ploidy = ploidy ;
		m_order_type = order_type ;
		m_missing = false ;
	}

	// Called once for each genotype (or haplotype) probability per sample.
	void set_value( uint32_t entry_i, double value ) {
		m_data[ entry_i ] = value ;
		if( entry_i == m_data.size() - 1 ) {
			// We have read all the data now, so let's compute the allele counts
			compute_expected_allele_counts( m_data, &m_expected_allele_counts ) ;
		}
	}

	// Ditto, but called if data is missing for this sample.
	void set_value( uint32_t entry_i, genfile::MissingValue value ) {
		m_data[ entry_i ] = -1 ;
		m_missing = true ;
		if( entry_i == m_data.size() - 1 ) {
			// We have read all the data now, so let's compute the allele counts
			compute_expected_allele_counts( m_data, &m_expected_allele_counts ) ;
		}
	}

	// If present with this signature, called once after all samples have been processed.
	void finalise() {
		// Nothing to do here.
	}
	
	// Report the results
	std::vector< double > const& expected_allele_counts() const {
		return m_expected_allele_counts ;
	}

private:
	std::vector< double > m_data ;
	std::size_t m_number_of_alleles ;

	// Used to keep track of what we're doing.
	std::size_t m_ploidy ;
	genfile::OrderType m_order_type ;
	bool m_missing ;
	
	// These fields are used to enumerate genotypes for the GT field.
	std::vector< uint16_t > m_genotype_allele_limits ;
	std::vector< uint16_t > m_genotype ;
	
	std::vector< double > m_expected_allele_counts ;

private:
	// Compute expected allele counts given genotype or haplotype probabilities
	// This implementation handles arbitrary ploidy and numbers of alleles.
	// This could be replaced with a simpler implementation if all samples are diploid.
	void compute_expected_allele_counts(
		std::vector< double > const& probs,
		std::vector< double >* result
	) {
		if( !m_missing ) {
			 if( m_order_type == genfile::ePerPhasedHaplotypePerAllele ) {
				compute_expected_allele_counts_phased( probs, result ) ;
			} else {
				compute_expected_allele_counts_unphased( probs, result ) ;
			}
		}
	}

	// compute expected allele counts given haplotype probabilities
	void compute_expected_allele_counts_phased(
		std::vector< double > const& probs,
		std::vector< double >* result
	) {
		assert( result->size() == m_number_of_alleles ) ;
		for( uint32_t i = 0; i < m_ploidy; ++i ) {
			uint32_t j = 0 ;
			for( ; j < m_number_of_alleles; ++j ) {
				(*result)[j] += probs[i*m_number_of_alleles+j] ;
			}
		}
	}

	// compute expected allele counts given haplotype probabilities
	void compute_expected_allele_counts_unphased(
		std::vector< double > const& probs,
		std::vector< double >* result
	) {
		assert( result->size() == m_number_of_alleles ) ;
		// Genotype probabilities are stored in colex order of the allele count representation.
		// Specifically, suppose we have m_ploidy = n chromosomes in total.
		// Genotypes are all ways to put n_alleles = k alleles into those chromosomes.
		// We represent these as k-vectors that sum to n (i.e. v=(v_i) where v_i is the count of allele i).
		// Colex order is lexicographical order of these vectors, reading them right-to-left.
		// E.g. for ploidy = 3 and 3 alleles, the order is
		// 3,0,0 = AAA
		// 2,1,0 = AAB
		// 1,2,0 = ABB
		// 0,3,0 = BBB
		// 2,0,1 = AAC
		// 1,1,1 = ABC
		// 0,2,1 = BBC
		// 1,0,2 = ACC
		// 0,1,2 = BCC
		// 0,0,3 = CCC
		// Here we enumerate these and accumulate, multiplying the probabilities by the allele counts.
		m_genotype_allele_limits.assign( (m_number_of_alleles-1), m_ploidy ) ;
		m_genotype.assign( m_number_of_alleles, 0 ) ;
		// Set up first genotype - all ref allele
		m_genotype[0] = m_ploidy ;

		// Iterate through all genotypes.
		for( std::size_t index = 0; true; ++index ) {
			// Accumulate probs.
			for( std::size_t k = 0; k < m_number_of_alleles; ++k ) {
				(*result)[k] += probs[index] * m_genotype[k] ;
			}
			
			// Move to next possible genotype
			std::size_t j = 0 ;
			for( ; j < (m_number_of_alleles-1); ++j ) {
				uint16_t value = m_genotype[j+1] ;
				if( value < m_genotype_allele_limits[ j ] ) {
					++m_genotype[j+1] ;
					--m_genotype[0] ;
					for( std::size_t k = 0; k < j; ++k ) {
						--m_genotype_allele_limits[k] ;
					}
					break ;
				} else {
					// allele count has reached its limit.
					// Reset it to zero.
					// Note that to get here all lower-order counts must be zero.
					m_genotype[j+1] = 0 ;
					m_genotype[0] += value ;
					for( std::size_t k = 0; k < j; ++k ) {
						m_genotype_allele_limits[k] += value ;
					}
				}
			}
			if( j == (m_number_of_alleles-1) ) {
				break ;
			}
		}
	}
} ;

void output_allele_counts(
	std::string const& rsid,
	std::vector< std::string > const& alleles,
	std::vector< double > const& expected_allele_counts
) {
	assert( expected_allele_counts.size() == alleles.size() ) ;

	std::cout << std::setprecision(2) << std::fixed ;
	std::cout << rsid << ":  " ;
	for( std::size_t i = 0; i < alleles.size(); ++i ) {
		std::cout << ((i>0) ? " " : "" ) << std::setw( 6 ) << std::right << alleles[i] ;
	}
	std::cout << " " << std::setw(6) << "total" ;

	double const total = std::accumulate( expected_allele_counts.begin(), expected_allele_counts.end(), 0.0 ) ;

	std::cout << "\n" ;
	std::cout << std::string( rsid.size() + 3, ' ' ) ;
	for( std::size_t i = 0; i < expected_allele_counts.size(); ++i ) {
		std::cout << ((i>0) ? " " : "" ) << std::setw( 6 ) << std::right << expected_allele_counts[i] ;
	}
	std::cout << " " << std::setw( 6 ) << std::right << total ;
	std::cout << "\n" ;
}

// This example program reads data from a bgen file
// and computes the expected count of each allele at each variant.
// Optionally, if a list of IDs are given, the index file will be used to
// restrict results to the IDs listed.
int main( int argc, char** argv ) {
	if( argc < 2 ) {
		std::cerr << "Usage: count_alleles <name of bgen file> [id1...]\n" ;
		exit(-1) ;
	}
	std::string const filename = argv[1] ;
	try {
		using namespace genfile ;
		bgen::View::UniquePtr bgenView = bgen::View::create( filename ) ;

		// If further arguments are given, use them as a query on the file.
		if( argc > 2 ) {
			genfile::bgen::IndexQuery::UniquePtr query = bgen::IndexQuery::create( filename + ".bgi" ) ;
			query
				->include_rsids( std::vector< std::string >( &argv[0] + 2, &argv[0] + argc ) )
				.initialise() ;
			bgenView->set_query( query ) ;
		}
		
		// Now iterate through variants.
		std::string chromosome ;
		uint32_t position ;
		std::string SNPID, rsid ;
		std::vector< std::string > alleles ;
		std::vector< std::vector< double > > probs ;

		AlleleCounter allele_counter ;
		while( bgenView->read_variant( &SNPID, &rsid, &chromosome, &position, &alleles ) ) {
			bgenView->read_genotype_data_block( allele_counter ) ;
			output_allele_counts( rsid, alleles, allele_counter.expected_allele_counts() ) ;
		}
		return 0 ;
	}
	catch( std::invalid_argument const& e ) {
		std::cerr << "!! Error: " << e.what() << ".\n" ;
		return -1 ;
	}
	catch( genfile::bgen::BGenError const& e ) {
		std::cerr << "!! Uh-oh, error parsing bgen file.\n" ;
		return -1 ;
	}
}

