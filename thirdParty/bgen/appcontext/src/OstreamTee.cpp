
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <map>
#include <string>
#include "appcontext/OstreamTee.hpp"

namespace appcontext {
	OstreamTee::~OstreamTee() {
		for( std::size_t i = 0; i < m_managed_streams.size(); ++i ) {
			delete m_managed_streams[i] ;
		}
	}

	void OstreamTee::add_stream( std::string const& name, std::ostream& stream ) {
		m_streams[ name ] = &stream ;
	}

	void OstreamTee::add_stream( std::string const& name, std::unique_ptr< std::ostream > stream ) {
		m_managed_streams.push_back( stream.release() ) ;
		m_streams[ name ] = m_managed_streams.back() ;
	}

	std::ostream& OstreamTee::operator[]( std::string const& name ) {
		std::map< std::string, std::ostream* >::iterator
			where = m_streams.find( name ) ;
		assert( where != m_streams.end() ) ;
		return *(where->second) ;
	}
    
	std::ostream& OstreamTee::operator[]( char const* name ) {
        return operator[]( std::string( name )) ;
    }
    

	OstreamTee& operator<<( OstreamTee& ostream_tee, OstreamTee::Manipulator t ) {
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

