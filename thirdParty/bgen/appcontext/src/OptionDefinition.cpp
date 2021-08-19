
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <sstream>
#include <cassert>
#include "appcontext/OptionDefinition.hpp"
#include "appcontext/OptionProcessor.hpp"

namespace appcontext {
	OptionDefinition::OptionDefinition():
		m_number_of_values_per_use( 0 ),
		m_minimum_multiplicity( 0 ),
		m_maximum_multiplicity( 1 ),
		m_default_values(),
		m_position( -1 ),
		m_hidden( false )
	{}

	OptionDefinition::OptionDefinition( OptionDefinition const& other ):
		m_group( other.m_group ),
		m_description( other.m_description ),
		m_number_of_values_per_use( other.m_number_of_values_per_use ),
		m_minimum_multiplicity( other.m_minimum_multiplicity ),
		m_maximum_multiplicity( other.m_maximum_multiplicity ),
		m_default_values( other.m_default_values ),
		m_position( other.m_position ),
		m_hidden( false )
	{}

	OptionDefinition& OptionDefinition::operator=( OptionDefinition const& other ) {
		m_group = other.m_group ;
		m_description = other.m_description ;
		m_number_of_values_per_use = other.m_number_of_values_per_use ;
		m_minimum_multiplicity = other.m_minimum_multiplicity ;
		m_maximum_multiplicity = other.m_maximum_multiplicity ;
		m_default_values = other.m_default_values ;
		m_position = other.m_position ;
		m_hidden = other.m_hidden ;
		return *this ;
	}

	void OptionDefinition::check_option_values( std::string const& option_name, std::vector< std::string > const& option_values ) const {
		// Check if number or values is too few...
		if( number_of_values_per_use() != eUntilNextOption ) {
			if( option_values.size() < m_minimum_multiplicity * number_of_values_per_use() ) {
				std::ostringstream ostr ;
				ostr << "Option \"" << option_name << "\" requires at least " << m_minimum_multiplicity * number_of_values_per_use() << " value" ;
				if( m_minimum_multiplicity * number_of_values_per_use() > 1 ) ostr << "s" ;
		
				throw OptionValueInvalidException( option_name, option_values, ostr.str() ) ;
			}

			// ...or too many.
			if( option_values.size() > m_maximum_multiplicity * number_of_values_per_use() ) {
				std::ostringstream ostr ;
				ostr << "Option \"" << option_name << "\" takes at most " << m_maximum_multiplicity * number_of_values_per_use() << " value" ;
				if( m_maximum_multiplicity * number_of_values_per_use() > 1 ) ostr << "s" ;

				throw OptionValueInvalidException( option_name, option_values, ostr.str() ) ;
			}
		}

		std::vector< value_checker_t >::const_iterator
			i = m_value_checkers.begin(),
			end_i = m_value_checkers.end() ;
		for( ; i != end_i; ++i ) {
			(*i)( option_name, option_values ) ;
		}
	}
}
