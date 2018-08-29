
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include "appcontext/progress_bar.hpp"

namespace appcontext {
	// Given a width and a number from 0 to 1 representing the amount of
	// progress through a task, generate a progress bar representing that amount.
	std::string get_progress_bar( std::size_t width, double progress ) {
		progress = std::min( std::max( progress, 0.0 ), 1.0 ) ;
		if( progress == progress && progress > 0.0 ) {
			std::size_t visible_progress = progress * width ;
			return
				"["
				+ std::string( std::size_t( visible_progress ), '*' )
				+ std::string( std::size_t( width - visible_progress ), ' ' )
				+ "]" ;
		} else {
			return
				"["
				+ std::string( std::size_t( width ), ' ' )
				+ "]" ;
		}
	}
}

