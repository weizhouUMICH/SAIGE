
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef APPCONTEXT_CMD_LINE_OPTION_PROCESSOR_HPP
#define APPCONTEXT_CMD_LINE_OPTION_PROCESSOR_HPP

#include <string>
#include <map>
#include <set>

#include "appcontext/ProgramFlow.hpp"
#include "appcontext/OptionProcessor.hpp"

namespace appcontext {
	struct CmdLineOptionProcessor: public OptionProcessor {
		virtual ~CmdLineOptionProcessor() ;

		void process( int argc, char** argv ) ;
		virtual std::string get_program_name() const = 0 ;
		virtual void declare_options( OptionProcessor& options ) ;
	} ;
}

#endif
