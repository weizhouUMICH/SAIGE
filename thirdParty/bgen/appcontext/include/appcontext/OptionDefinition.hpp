
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_APPCONTEXT_OPTION_DEFINITION_HPP
#define QCTOOL_APPCONTEXT_OPTION_DEFINITION_HPP

#include <vector>
#include <limits>
#include <set>
#include <string>
#include <sstream>
#include <cassert>

namespace appcontext {
	struct OptionDefinition {
		enum { eUntilNextOption = 2147483647 } ;
		public:
			OptionDefinition() ;
			OptionDefinition( OptionDefinition const& other ) ;
			OptionDefinition& operator=( OptionDefinition const& other ) ;

			typedef void (*value_checker_t)( std::string const&, std::vector< std::string > const& ) ;
			typedef std::vector< std::string > (*value_preprocessor_t)( std::string const&, std::vector< std::string > const& ) ;

			std::string const& group() const { return m_group ; }
			std::string description() const { return m_description ; } 
			bool is_required() const { return m_minimum_multiplicity > 0 ; }
			unsigned int number_of_values_per_use() const { return m_number_of_values_per_use ; }
		
			bool takes_values() const { return m_number_of_values_per_use > 0  ; }
			std::vector< value_checker_t > value_checkers() const { return m_value_checkers ; } 
			bool has_default_value() const { return m_default_values.size() > 0 ; }
			std::string default_value() const { assert( m_default_values.size() == 1 ) ; return m_default_values[0] ; }
			std::vector< std::string > default_values() const { assert( m_default_values.size() > 0 ) ; return m_default_values ; }
			bool takes_value_by_position() const { return m_position > 0 ; }
			int position() const { return m_position ; }
			bool hidden() const { return m_hidden ; }
		
			OptionDefinition& set_group( std::string const& group ) { m_group = group ; return *this ; }
			OptionDefinition& set_description( char const* desc ) { m_description = desc ; return *this ; }
			OptionDefinition& set_description( std::string const& desc ) { m_description = desc ; return *this ; }
			OptionDefinition& set_is_required() {
				m_minimum_multiplicity = std::max( m_minimum_multiplicity, 1u ) ;
				m_maximum_multiplicity = std::max( m_maximum_multiplicity, 1u ) ;
				return *this ;
			}
			OptionDefinition& set_minimum_multiplicity( unsigned int multiplicity ) {
				m_minimum_multiplicity = multiplicity ;
				m_maximum_multiplicity = std::max( m_maximum_multiplicity, multiplicity ) ;
				return *this ;
			}
			OptionDefinition& set_maximum_multiplicity( unsigned int multiplicity ) {
				assert( multiplicity >= m_minimum_multiplicity ) ;
				m_maximum_multiplicity = multiplicity ;
				return *this ;
			}
			OptionDefinition& set_takes_values( unsigned int n ) {
				m_number_of_values_per_use = n ;
				return *this ;
			}
			OptionDefinition& set_takes_single_value() {
				m_number_of_values_per_use = 1u ;	
				m_maximum_multiplicity = 1u ;
				m_minimum_multiplicity = std::min( m_minimum_multiplicity, m_maximum_multiplicity ) ;
				return *this ;
			}
			OptionDefinition& set_takes_values_until_next_option() {
				m_number_of_values_per_use = eUntilNextOption ;	
				return *this ;
			}
			OptionDefinition& add_value_checker( value_checker_t value_checker ) {
				assert( value_checker != 0 ) ;
				m_value_checkers.push_back( value_checker ) ; 
				return *this ; 
			}
			template< typename T > OptionDefinition& set_default_value( T const& value ) {
				std::ostringstream aStream ;
				aStream << value ;
				m_default_values.push_back( aStream.str() ) ;
				return *this ;
			}
			OptionDefinition& set_takes_value_by_position( int position ) {
				m_position = position ;
				return *this ;
			}
			OptionDefinition& set_hidden() {
				m_hidden = true ;
				return *this ;
			}

			void check_option_values( std::string const& option_name, std::vector< std::string > const& option_values ) const ;

	    private:

			std::string m_group ;
			std::string m_description ;
			unsigned int m_number_of_values_per_use ;
			unsigned int m_minimum_multiplicity, m_maximum_multiplicity ;
			std::vector< value_checker_t > m_value_checkers ;
			std::vector< std::string > m_default_values ;
			int m_position ;
			bool m_hidden ;
	} ;
}

#endif

