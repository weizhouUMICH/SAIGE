#include <sstream>
#include <map>
#include "genfile/bgen/View.hpp"
#include "genfile/bgen/IndexQuery.hpp"
#include <Rcpp.h>

// #define DEBUG 1

namespace {
	template< typename Integer >
	std::string atoi( Integer const value ) {
		std::ostringstream stream ;
		stream << value ;
		return stream.str() ;
	}
	
	struct DataSetter {
		typedef Rcpp::IntegerVector IntegerVector ;
		typedef Rcpp::NumericVector NumericVector ;
		typedef Rcpp::LogicalVector LogicalVector ;
		typedef Rcpp::Dimension Dimension ;

		DataSetter(
			IntegerVector* ploidy,
			Dimension const& ploidy_dimension,
			NumericVector* data,
			Dimension const& data_dimension,
			LogicalVector* phased,
			std::size_t variant_i,
			std::map< std::size_t, std::size_t > const& requested_samples
		):
			m_ploidy( ploidy ),
			m_ploidy_dimension( ploidy_dimension ),
			m_data( data ),
			m_data_dimension( data_dimension ),
			m_phased( phased ),
			m_variant_i( variant_i ),
			m_requested_samples( requested_samples ),
			m_requested_sample_i( m_requested_samples.begin() ),
			m_storage_i( 0 ),
			m_order_type( genfile::eUnknownOrderType )
		{
			assert( m_data_dimension[0] == m_ploidy_dimension[0] ) ;
			assert( m_data_dimension[1] == m_ploidy_dimension[1] ) ;
			assert( m_data_dimension[1] == m_requested_samples.size() ) ;
			assert( m_variant_i < m_data_dimension[0] ) ;
			assert( m_data_dimension[2] >= 3 ) ;
			assert( m_phased->size() == m_data_dimension[0] ) ;
		}
	
		// Called once allowing us to set storage.
		void initialise( std::size_t number_of_samples, std::size_t number_of_alleles ) {
		}

		// If present with this signature, called once after initialise()
		// to set the minimum and maximum ploidy and numbers of probabilities among samples in the data.
		// This enables us to set up storage for the data ahead of time.
		void set_min_max_ploidy(
			genfile::bgen::uint32_t min_ploidy, genfile::bgen::uint32_t max_ploidy,
			genfile::bgen::uint32_t min_entries, genfile::bgen::uint32_t max_entries
		) {
			if( max_entries > m_data_dimension[2] ) {
				throw std::invalid_argument( "max_entries=" + atoi( max_entries )
					+ " (expected at most " + atoi( m_data_dimension[2] ) + ")" ) ;
			}
		}

		// Called once per sample to determine whether we want data for this sample
		bool set_sample( std::size_t i ) {
			if( m_requested_sample_i->first == i ) {
				m_storage_i = m_requested_sample_i->second ;
				++m_requested_sample_i ;
#if DEBUG
				std::cerr << "DataSetter::set_sample(): sample " << i << " has storage index " << m_storage_i << ".\n" ;
#endif
				return true ;
			} else {
				// Don't want this sample
				return false ;
			}
		}

		// Called once per sample to set the number of probabilities that are present.
		void set_number_of_entries(
			std::size_t ploidy,
			std::size_t number_of_entries,
			genfile::OrderType order_type,
			genfile::ValueType value_type
		) {
			if( value_type != genfile::eProbability ) {
				throw std::invalid_argument(
					"value_type ("
					+ atoi( value_type ) + ", expected "
					+ atoi( genfile::eProbability ) + "=genfile::eProbability)"
				) ;
			}
			if( m_order_type == genfile::eUnknownOrderType ) {
				m_order_type = order_type ;
				(*m_phased)( m_variant_i ) = ( m_order_type == genfile::ePerPhasedHaplotypePerAllele ) ;
			} else {
				assert( order_type == m_order_type ) ;
			}
			int const index = m_variant_i + m_storage_i * m_ploidy_dimension[0] ;
			(*m_ploidy)[ index ] = ploidy ;
		}

		void set_value( genfile::bgen::uint32_t entry_i, double value ) {
			int const index = m_variant_i + m_storage_i * m_data_dimension[0] + entry_i * m_data_dimension[0] * m_data_dimension[1] ;
#if DEBUG
			std::cerr << "Setting data for index " << m_variant_i << ", " << m_storage_i << ", " << entry_i << ": index " << index << "...\n" << std::flush ;
#endif
			(*m_data)[ index ] = value ;
		}

		void set_value( genfile::bgen::uint32_t entry_i, genfile::MissingValue value ) {
			int const index = m_variant_i + m_storage_i * m_data_dimension[0] + entry_i * m_data_dimension[0] * m_data_dimension[1] ;
#if DEBUG
			std::cerr << "Setting data for index " << m_variant_i << ", " << m_storage_i << ", " << entry_i << ": index " << index << "...\n" << std::flush ;
#endif
			(*m_data)[ index ] = NA_REAL ;
		}

		void finalise() {
		// nothing to do
		}

	private:
		Rcpp::IntegerVector* m_ploidy ;
		Rcpp::Dimension const m_ploidy_dimension ;
		Rcpp::NumericVector* m_data ;
		Rcpp::Dimension const m_data_dimension ;
		Rcpp::LogicalVector* m_phased ;

		std::size_t const m_variant_i ;

		std::map< std::size_t, std::size_t > const& m_requested_samples ;
		std::map< std::size_t, std::size_t >::const_iterator m_requested_sample_i ;
		std::size_t m_storage_i ;

		genfile::OrderType m_order_type ;
	} ;

	struct set_sample_names {
		typedef std::map< std::size_t, std::size_t > SampleIndexMap ;
		
		set_sample_names( std::vector< std::string >* result, SampleIndexMap* sample_indices ):
			m_result( result ),
			m_sample_indices( sample_indices ),
			m_index(0)
		{
			assert( result != 0 ) ;
			assert( sample_indices != 0 ) ;
			assert( sample_indices->size() == result->size() ) ;
		}

		void operator()( std::string const& value ) {
			m_sample_indices->insert( std::make_pair( m_index, m_index ) ) ;
			(*m_result)[m_index++] = value ;
		}
	private:
		std::vector< std::string >* m_result ;
		SampleIndexMap* m_sample_indices ;
		std::size_t m_index ;
	} ;
	
	struct set_requested_sample_names {
		typedef std::map< std::string, std::size_t > RequestedSamples ;
		typedef std::map< std::size_t, std::size_t > SampleIndexMap ;
		
		set_requested_sample_names(
			std::vector< std::string >* result,
			SampleIndexMap* sample_indices,
			RequestedSamples const& requested_samples
		):
			m_result( result ),
			m_sample_indices( sample_indices ),
			m_requested_samples( requested_samples ),
			m_index(0),
			m_value_index(0)
		{
			assert( result != 0 ) ;
			assert( sample_indices != 0 ) ;
			assert( sample_indices.size() == 0 ) ;
			assert( result->size() == requested_samples.size() ) ;
		}
		
		void operator()( std::string const& value ) {
			RequestedSamples::const_iterator where = m_requested_samples.find( value ) ;
			if( where != m_requested_samples.end() ) {
				(*m_result)[ where->second ] = value ;
				m_sample_indices->insert( std::make_pair( m_value_index, where->second ) ) ;
			}
			++m_value_index ;
		}
	private:
		std::vector< std::string >* m_result ;
		SampleIndexMap* m_sample_indices ;
		RequestedSamples const& m_requested_samples ;
		std::size_t m_index ;
		std::size_t m_value_index ;
	} ;
	
	genfile::bgen::View::UniquePtr construct_view(
		std::string const& filename,
		std::string const& index_filename,
		Rcpp::DataFrame const& ranges,
		std::vector< std::string > const& rsids = std::vector< std::string >()
	) {
		using namespace genfile::bgen ;
		using namespace Rcpp ;

		View::UniquePtr view = View::create( filename ) ;
		{
			StringVector const& chromosome = ranges["chromosome"] ;
			IntegerVector const& start = ranges["start"] ;
			IntegerVector const& end = ranges["end"] ;
			
			IndexQuery::UniquePtr query = IndexQuery::create( index_filename ) ;
			for( int i = 0; i < ranges.nrows(); ++i ) {
				if( end[i] < start[i] ) {
					throw std::invalid_argument( "Range (" + chromosome[i] + ":" + atoi( start[i] ) + "-" + atoi( end[i] ) + ") is malformed." ) ;
				}
				query->include_range( IndexQuery::GenomicRange( std::string( chromosome[i] ), start[i], end[i] )) ;
			}
			query->include_rsids( rsids ) ;
			query->initialise() ;
			view->set_query( query ) ;
		}
		return view ;
	}
}

void get_all_samples(
	genfile::bgen::View const& view,
	std::size_t* number_of_samples,
	std::vector< std::string >* sampleNames,
	std::map< std::size_t, std::size_t >* requestedSamplesByIndexInDataIndex
) {
	*number_of_samples = view.number_of_samples() ;
	sampleNames->resize( *number_of_samples ) ;
	view.get_sample_ids( set_sample_names( sampleNames, requestedSamplesByIndexInDataIndex ) ) ;
}

void get_requested_samples(
	genfile::bgen::View const& view,
	Rcpp::StringVector const& requestedSamples,
	std::size_t* number_of_samples,
	std::vector< std::string >* sampleNames,
	std::map< std::size_t, std::size_t >* requestedSamplesByIndexInDataIndex
) {
	// convert requested sample IDs to a map of requested indices.
	std::map< std::string, std::size_t > requestedSamplesByName ;
	for( std::size_t i = 0; i < requestedSamples.size(); ++i ) {
		requestedSamplesByName.insert( std::map< std::string, std::size_t >::value_type( std::string( requestedSamples[i] ), i )) ;
	}
	if( requestedSamplesByName.size() != requestedSamples.size() ) {
		throw std::invalid_argument(
			"load_unsafe(): requiredSamples: expected a list of unique samples with no repeats."
		) ;
	}

	*number_of_samples = requestedSamples.size() ;
	sampleNames->resize( requestedSamples.size() ) ;
	view.get_sample_ids( set_requested_sample_names( sampleNames, requestedSamplesByIndexInDataIndex, requestedSamplesByName ) ) ;

	// Check each requested sample has been asked for exactly once
	// We count distinct samples, among those requested, that we've found in the data
	// And we also count the min and max index of those samples.
	// If min = 0 and max = (#requested samples-1) and each sample was unique, we're ok.
	std::set< std::size_t > checkSamples ;
	std::size_t minIndex = std::numeric_limits< std::size_t >::max() ;
	std::size_t maxIndex = 0 ;
	
	for(
		std::map< std::size_t, std::size_t >::const_iterator p = requestedSamplesByIndexInDataIndex->begin();
		p != requestedSamplesByIndexInDataIndex->end();
		++p
	 ) {
		checkSamples.insert( p->second ) ;
		minIndex = std::min( minIndex, p->second ) ;
		maxIndex = std::max( maxIndex, p->second ) ;
	}
	if( checkSamples.size() != requestedSamples.size() || minIndex != 0 || maxIndex != (requestedSamples.size()-1) ) {
		// Huh.  To be most useful, let's print diagnostics
		std::cerr << "!! Uh-oh: requested sample indices (data, request) are:\n" ;
		for(
			std::map< std::size_t, std::size_t >::const_iterator p = requestedSamplesByIndexInDataIndex->begin();
			p != requestedSamplesByIndexInDataIndex->end();
			++p
		 ) {
			std::cerr << p->first << ", " << p->second << ".\n" ;
		}
		
		throw std::invalid_argument(
			"load_unsafe(): requiredSamples contains a sample not present in the data, or data contains a repeated sample ID."
		) ;
	}
}

Rcpp::List load_unsafe(
	std::string const& filename,
	std::string const& index_filename,
	Rcpp::DataFrame const& ranges,
	std::vector< std::string > const& requested_rsids,
	std::size_t max_entries_per_sample,
	Rcpp::StringVector const* const requestedSamples
) {
	using namespace genfile::bgen ;
	using namespace Rcpp ;

	View::UniquePtr view = construct_view( filename, index_filename, ranges, requested_rsids ) ;

	std::size_t const number_of_variants = view->number_of_variants() ;
	std::size_t number_of_samples = 0 ;

	// Build list of sample names as a std::vector.
	// Note: Rcpp will convert it to StringVector on return
	std::vector< std::string > sampleNames ;
	std::map< std::size_t, std::size_t > requestedSamplesByIndexInDataIndex ;
	
	if( requestedSamples ) {
		get_requested_samples( *view, *requestedSamples, &number_of_samples, &sampleNames, &requestedSamplesByIndexInDataIndex ) ;
	} else {
		get_all_samples( *view, &number_of_samples, &sampleNames, &requestedSamplesByIndexInDataIndex ) ;
	}

	// Declare storage for all the things we need
	StringVector chromosomes( number_of_variants ) ;
	IntegerVector positions( number_of_variants ) ;
	StringVector rsids( number_of_variants ) ;
	IntegerVector number_of_alleles( number_of_variants ) ;
	StringVector allele0s( number_of_variants ) ;
	StringVector allele1s( number_of_variants ) ;

	Dimension data_dimension = Dimension( number_of_variants, number_of_samples, max_entries_per_sample ) ;
	Dimension ploidy_dimension = Dimension( number_of_variants, number_of_samples ) ;

	NumericVector data = NumericVector( data_dimension, NA_REAL ) ;
	IntegerVector ploidy = IntegerVector( ploidy_dimension, NA_INTEGER ) ;
	LogicalVector phased = LogicalVector( number_of_variants, NA_LOGICAL ) ;

	std::string SNPID, rsid, chromosome ;
	genfile::bgen::uint32_t position ;
	std::vector< std::string > alleles ;

	// Iterate through variants
	for( std::size_t variant_i = 0; variant_i < number_of_variants; ++variant_i ) {
		view->read_variant( &SNPID, &rsid, &chromosome, &position, &alleles ) ;
		chromosomes[variant_i] = chromosome ;
		positions[variant_i] = position ;
		rsids[variant_i] = rsid ;
		number_of_alleles[variant_i] = alleles.size() ;
		allele0s[variant_i] = alleles[0] ;
		allele1s[variant_i] = alleles[1] ;

		DataSetter setter(
			&ploidy, ploidy_dimension,
			&data, data_dimension,
			&phased,
			variant_i,
			requestedSamplesByIndexInDataIndex
		) ;
		
		view->read_genotype_data_block( setter ) ; // will be fixed later
	}

	DataFrame variants = DataFrame::create(
		Named("chromosome") = chromosomes,
		Named("position") = positions,
		Named("rsid") = rsids,
		Named("number_of_alleles") = number_of_alleles,
		Named("allele0") = allele0s,
		Named("allele1") = allele1s
	) ;
	variants.attr( "row.names" ) = rsids ;

	StringVector genotypeNames(max_entries_per_sample) ;
	for( std::size_t i = 0; i < max_entries_per_sample; ++i ) {
		genotypeNames[i] = "g=" + atoi(i) ;
	}

	List dataNames = List(3) ;
	dataNames[0] = rsids ;
	dataNames[1] = sampleNames ;
	dataNames[2] = genotypeNames ;

	List ploidyNames = List(2) ;
	ploidyNames[0] = rsids ;
	ploidyNames[1] = sampleNames ;

	data.attr( "dimnames" ) = dataNames ;
	ploidy.attr( "dimnames" ) = ploidyNames ;

	List result ;
	result[ "variants" ] = variants ;
	result[ "samples" ] = sampleNames ;
	result[ "ploidy" ] = ploidy ;
	result[ "phased" ] = phased ;
	result[ "data" ] = data ;

	return( result ) ;
}

// [[Rcpp::export]]
Rcpp::List load(
	std::string const& filename,
	std::string const& index_filename,
	Rcpp::DataFrame const& ranges,
	std::vector< std::string > const& rsids,
	std::size_t max_entries_per_sample
) {
	try {
		return load_unsafe( filename, index_filename, ranges, rsids, max_entries_per_sample, 0 ) ;
	}
	catch( std::exception const& e ) {
		forward_exception_to_r( e ) ;
	}
	catch( ... ) {
		::Rf_error("A C++ exception occurred (unknown reason)") ;
	}
	return Rcpp::List() ;
}

// [[Rcpp::export]]
Rcpp::List load(
	std::string const& filename,
	std::string const& index_filename,
	Rcpp::DataFrame const& ranges,
	std::vector< std::string > const& rsids,
	std::size_t max_entries_per_sample,
	Rcpp::StringVector const& requestedSamples
) {
	try {
		return load_unsafe( filename, index_filename, ranges, rsids, max_entries_per_sample, &requestedSamples ) ;
	}
	catch( std::exception const& e ) {
		forward_exception_to_r( e ) ;
	}
	catch( ... ) {
		::Rf_error("A C++ exception occurred (unknown reason)") ;
	}
	return Rcpp::List() ;
}

