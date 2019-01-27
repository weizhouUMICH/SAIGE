
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#include <set>
#include <utility>
#include <functional>
#include "genfile/types.hpp"
#include "genfile/MissingValue.hpp"

/*
* ProbSetCheck
* This class checks sequence and values of probabilities set using
* the bgen API
*/
struct ProbSetCheck {
	typedef std::function< double( std::size_t i, std::size_t g ) > GetExpectedProbs ;
	enum State { eNone, eSetNumberOfSamples, eSetSample, eSetNumberOfEntries, eSetValue, eFinalised } ;

	ProbSetCheck(
		std::size_t n,
		GetExpectedProbs get_expected_probs
	) ;

	~ProbSetCheck() throw() ;

	void expect_fail() ;

	void initialise( std::size_t nSamples, std::size_t nAlleles ) ;
	void set_min_max_ploidy( uint32_t min_ploidy, uint32_t max_ploidy, uint32_t min_entries, uint32_t max_entries ) ;
	bool set_sample( std::size_t i ) ;
	void set_number_of_entries( std::size_t ploidy, std::size_t n, genfile::OrderType const order_type, genfile::ValueType const value_type ) ;
	void set_value( uint32_t value_i, genfile::MissingValue const value ) ;
	void set_value( uint32_t value_i, double const value ) ;
	void finalise() ;
	
private:
	std::size_t m_number_of_samples ;
	GetExpectedProbs m_get_expected_probs ;
	std::size_t m_sample_i ;
	std::size_t m_number_of_entries ;
	std::size_t m_entry_i ;
	genfile::OrderType m_order_type ;
	std::set< std::pair< std::size_t, std::size_t > > m_set_values ;

	State m_state ;
	bool m_expect_fail ;
} ;
