
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef TIMER_HPP
#define TIMER_HPP

#include <string>
#include <cassert>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <boost/timer/timer.hpp>

namespace appcontext {
	struct Timer
	{
		double elapsed() const ;
		void restart() ;
		std::string display() const ;
	private:
		boost::timer::cpu_timer m_boost_timer ;
	} ;
}

#endif
