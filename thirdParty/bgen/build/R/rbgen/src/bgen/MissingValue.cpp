
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include "genfile/MissingValue.hpp"

namespace genfile {
	std::ostream& operator<<( std::ostream& o, MissingValue const& ) {
		return o << "NA" ;
	}
	
	bool MissingValue::operator==( MissingValue const& ) const {
		return true ;
	}

	bool MissingValue::operator<( MissingValue const& ) const {
		return false ;
	}

	bool MissingValue::operator<=( MissingValue const& ) const {
		return true ;
	}
}
