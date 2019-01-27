
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <exception>
#include <string>

#include "appcontext/ProgramFlow.hpp"
#include "appcontext/CmdLineOptionProcessor.hpp"

namespace appcontext {
	CmdLineOptionProcessor::~CmdLineOptionProcessor() {
	}

	void CmdLineOptionProcessor::process( int argc, char** argv ) {
		declare_options( *this ) ;
		OptionProcessor::process( argc, argv ) ;
	}

	void CmdLineOptionProcessor::declare_options( OptionProcessor& options ) {
	}
}
