
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef APPCONTEXT_OSTREAM_TEE_HPP
#define APPCONTEXT_OSTREAM_TEE_HPP

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <memory>
#include <cassert>

namespace appcontext {
	class OstreamTee {
	public:
		OstreamTee() {} ;
		~OstreamTee() ;

		void add_stream( std::string const& name, std::ostream& stream ) ; 
		void add_stream( std::string const& name, std::unique_ptr< std::ostream > stream ) ; 
		std::ostream& operator[]( std::string const& name ) ;
		std::ostream& operator[]( char const* name ) ;

		typedef std::ostream& (*Manipulator)( std::ostream& ) ;

		template< typename T >
		friend OstreamTee& operator<<( OstreamTee & ostream_tee, T const& t ) ;

		friend OstreamTee& operator<<( OstreamTee & ostream_tee, Manipulator t ) ;

	private:
		std::map< std::string, std::ostream* > m_streams ;
		std::vector< std::ostream* > m_managed_streams ;
		
		OstreamTee( OstreamTee const& ) ;
		OstreamTee& operator=( OstreamTee const& other ) ;
	} ;

	template< typename T >
	OstreamTee& operator<<( OstreamTee& ostream_tee, T const& t ) {
		for(
			std::map< std::string, std::ostream* >::const_iterator i = ostream_tee.m_streams.begin() ;
			i != ostream_tee.m_streams.end() ;
			++i
		) {
			(*(i->second)) << t ;
		}

		return ostream_tee ;
	}
}
#endif
