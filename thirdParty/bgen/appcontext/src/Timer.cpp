
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <cassert>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <boost/timer/timer.hpp>
#include "appcontext/Timer.hpp"

namespace appcontext {
	double Timer::elapsed() const {
		return double( m_boost_timer.elapsed().wall ) / 1000000000.0 ;
	}
	
	void Timer::restart() {
		m_boost_timer.start() ;
	}
	
	std::string Timer::display() const {
		std::ostringstream os ;
		os << std::fixed << std::setprecision(1) << elapsed() << "s" ;
		return os.str() ;
	}
}
