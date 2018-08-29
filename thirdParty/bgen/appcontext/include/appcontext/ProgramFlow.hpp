
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef APPCONTEXT_PROGRAM_FLOW_HPP
#define APPCONTEXT_PROGRAM_FLOW_HPP

#include <exception>

namespace appcontext {
	struct HaltProgramWithReturnCode: public std::exception {
		HaltProgramWithReturnCode( int return_code = -1 ) : m_return_code( return_code ) {}
		~HaltProgramWithReturnCode() throw() {}
		char const* what() const throw() { return "HaltProgramWithReturnCode" ; }
		int return_code() const { return m_return_code ; }
	private:
		int m_return_code ;
	} ;
}

#endif
