
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef APPLICATION_CONTEXT_HPP
#define APPLICATION_CONTEXT_HPP

#include <memory>
#include <string>
#include <boost/function.hpp>
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/Timer.hpp"
#include "appcontext/OstreamTee.hpp"
#include "appcontext/CmdLineUIContext.hpp"
#include "appcontext/OptionProcessor.hpp"

namespace appcontext {
	struct ApplicationContext
	{
	public:
		typedef appcontext::UIContext UIContext ;
		typedef boost::function< void( appcontext::OptionProcessor const& ) > OptionChecker ;
	public:
		ApplicationContext(
			std::string const& application_name,
			std::unique_ptr< OptionProcessor > options,
			int argc,
			char** argv,
			std::string const& log_filename,
			OptionChecker optionChecker = OptionChecker()
		) ;

		ApplicationContext(
			std::string const& application_name,
			std::string const& application_revision,
			std::unique_ptr< OptionProcessor > options,
			int argc,
			char** argv,
			std::string const& log_filename,
			OptionChecker optionChecker = OptionChecker()
		) ;

		virtual ~ApplicationContext() ;
		OptionProcessor& options() const ;
		virtual UIContext& ui() const ;

		std::string const& application_name() const ;

	private:
		
		void process_options( int argc, char** argv, std::string const& log_option, OptionChecker checker ) ;
		void construct_logger( std::string const& log_option ) ;
		void write_start_banner() ;
		void write_end_banner() ;
	
		std::string const m_application_name ;
		std::string const m_application_version ;
		std::unique_ptr< OptionProcessor > m_options ;
		std::unique_ptr< UIContext > m_ui ;
	} ;
}

#endif
