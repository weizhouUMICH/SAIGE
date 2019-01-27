
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <memory>
#include "genfile/bgen/bgen.hpp"
#include "genfile/bgen/View.hpp"

// DosageSetter is a callback object appropriate for passing to bgen::read_genotype_data_block() or
// the synonymous method of genfile::bgen::View.
// Its job is to compute genotype dosage (dosage of the 2nd allele).
// This implementation assumes samples are diploid, data unphased, and variants are biallelic
// and will assert() if not.
struct DosageSetter {
	typedef std::vector< double > Data ;
	DosageSetter():
		m_sample_i(0)
	{}
		
	// Called once allowing us to set storage.
	void initialise( std::size_t number_of_samples, std::size_t number_of_alleles ) {
		m_result = 0 ;
		m_n = 0 ;
	}
	
	// If present with this signature, called once after initialise()
	// to set the minimum and maximum ploidy and numbers of probabilities among samples in the data.
	// This enables us to set up storage for the data ahead of time.
	void set_min_max_ploidy( uint32_t min_ploidy, uint32_t max_ploidy, uint32_t min_entries, uint32_t max_entries ) {
		assert( min_ploidy == 2 ) ;
		assert( max_ploidy == 2 ) ;
	}
	
	// Called once per sample to determine whether we want data for this sample
	bool set_sample( std::size_t i ) {
		m_sample_i = i ;
		// Yes, here we want info for all samples.
		return true ;
	}
	
	// Called once per sample to set the number of probabilities that are present.
	void set_number_of_entries(
		std::size_t ploidy,
		std::size_t number_of_entries,
		genfile::OrderType order_type,
		genfile::ValueType value_type
	) {
		assert( value_type == genfile::eProbability ) ;
		assert( order_type == genfile::ePerUnorderedGenotype ) ;
		assert( ploidy == 2 ) ;
		assert( number_of_entries == 3 ) ;

		++m_n ;
	}

	// Called once for each genotype (or haplotype) probability per sample.
	void set_value( uint32_t entry_i, double value ) {
		m_result += entry_i * value ;
	}

	// Ditto, but called if data is missing for this sample.
	void set_value( uint32_t, genfile::MissingValue value ) {
	}

	// If present with this signature, called once after all data has been set.
	void finalise() {
		// nothing to do in this implementation.
		m_result /= m_n ;
	}

	double const& result() const { return m_result ; }
private:
	double m_result ;
	std::size_t m_n ;
	std::size_t m_sample_i ;
	std::size_t m_entry_i ;
} ;

void process_data_using_dosage_setter( genfile::bgen::View& view ) ;
void process_data_using_lookup_table( genfile::bgen::View& view ) ;
std::vector< double > compute_lookup_table() ;

// This example program reads data from a bgen file specified as the first argument
// and outputs it as a VCF file.
int main( int argc, char** argv ) {
	std::ios_base::sync_with_stdio( false ) ;
	if( argc < 2 ) {
		std::cerr << "You must supply an argument, the name of the bgen file to process.\n" ;
		std::cerr << "Usage: compute_expected_dosage <bgen filename> [<method>]\n" ;
		std::cerr << "Where <method> can be \"default\" or \"lookup-table\".\n" ;
		exit(-1) ;
	}
	if( argc > 3 ) {
		std::cerr << "Usage: compute_expected_dosage <bgen filename> [<method>]\n" ;
		std::cerr << "Where <method> can be \"default\" or \"lookup-table\".\n" ;
		std::cerr << "The lookup-table method is only implemented for 8-bit BGEN files.\n" ;
		exit(-1) ;
	}
	std::string const filename = argv[1] ;
	std::string mode = "default" ;
	if( argc == 3 ) {
		mode = argv[2] ;
	}
	try {
		genfile::bgen::View view( filename ) ;
		if( mode == "default" ) {
			process_data_using_dosage_setter( view ) ;
		} else if( mode == "lookup-table" ) {
			process_data_using_lookup_table( view ) ;
		} else {
			std::cerr << "method \"" + mode + "\" is not supported.  Use \"default\" or \"lookup-table\".\n" ;
			exit(-1) ;
		}
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

void process_data_using_dosage_setter( genfile::bgen::View& view ) {
	std::string SNPID, rsid, chromosome ;
	uint32_t position ;
	std::vector< std::string > alleles ;

	std::vector< double > m_data ;
	for( std::size_t i = 0; i < view.number_of_variants(); ++i ) {
		bool success = view.read_variant(
			&SNPID, &rsid, &chromosome, &position, &alleles
		) ;
		assert( success ) ;
		
		DosageSetter dc ;
		view.read_genotype_data_block( dc ) ;
		
		std::cout << rsid << ": " << dc.result() << "\n" ;	
	}
}


void process_data_using_lookup_table( genfile::bgen::View& view ) {
	std::vector< double > const lookup_table = compute_lookup_table() ;

	std::string SNPID, rsid, chromosome ;
	uint32_t position ;
	std::vector< std::string > alleles ;

	genfile::bgen::v12::GenotypeDataBlock pack ;
	for( std::size_t i = 0; i < view.number_of_variants(); ++i ) {
		bool success = view.read_variant(
			&SNPID, &rsid, &chromosome, &position, &alleles
		) ;
		assert( success ) ;
		
		view.read_and_unpack_v12_genotype_data_block( &pack ) ;
		assert( pack.numberOfAlleles == 2 ) ;
		assert( pack.ploidyExtent[0] == 2 ) ;
		assert( pack.ploidyExtent[1] == 2 ) ;
		assert( pack.phased == false ) ;
		assert( pack.bits == 8 ) ;
	
		std::size_t n = 0 ;
		double totalDosage = 0.0 ;
		genfile::byte_t const* buffer = pack.buffer ;
		genfile::byte_t const* ploidy = pack.ploidy ;
		for( int sample = 0; buffer < pack.end; buffer += 2, ++ploidy, ++sample ) {
			if( !(*ploidy & 0x80 )) {
				totalDosage += lookup_table[ *reinterpret_cast< uint16_t const* >( buffer ) ] ;
				++n ;
			}
		}
		std::cout << rsid << ": " << (totalDosage / n) << "\n" ;
	}
}

std::vector< double > compute_lookup_table() {
	std::vector< double > result( 256*256 ) ;
	uint32_t p0, p1 ;
	for( p0 = 0; p0 < 256; ++p0 ) {
		for( p1 = 0; p1 < (256-p0); ++p1 ) {
			uint32_t p2 = 255 - p0 - p1 ;
			double dosage = (2 * double(p2) + double(p1)) / 255.0 ;
			result[ (p1 << 8) | p0 ] = dosage ;
		}
	}
	return result ;
}
