
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <memory>
#include <string>
#include <fstream>
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/OstreamTee.hpp"
#include "appcontext/progress_bar.hpp"
#include "appcontext/ApplicationContext.hpp"
#include "appcontext/null_ostream.hpp"
#include "appcontext/ProgramFlow.hpp"

namespace appcontext {
	ApplicationContext::ApplicationContext(
		std::string const& application_name,
		std::unique_ptr< OptionProcessor > options,
		int argc,
		char** argv,
		std::string const& log_option,
		OptionChecker checker
	 ):
		m_application_name( application_name ),
		m_application_version( "" ),
		m_options( std::move( options ) ),
		m_ui( new appcontext::CmdLineUIContext )
	{
		process_options( argc, argv, log_option, checker ) ;
		write_start_banner() ;
	}

	ApplicationContext::ApplicationContext(
		std::string const& application_name,
		std::string const& application_version,
		std::unique_ptr< OptionProcessor > options,
		int argc,
		char** argv,
		std::string const& log_option,
		OptionChecker checker
	 ):
		m_application_name( application_name ),
		m_application_version( application_version ),
		m_options( std::move( options ) ),
		m_ui( new appcontext::CmdLineUIContext )
	{
		process_options( argc, argv, log_option, checker ) ;
		write_start_banner() ;
	}
	
	void ApplicationContext::process_options( int argc, char** argv, std::string const& log_option, OptionChecker checker ) {
		try {
			ui().logger().add_stream( "screen", std::cerr ) ;
			m_options->process( argc, argv ) ;
			if( checker ) {
				checker( *m_options ) ;
			}
			construct_logger( log_option ) ;
		}
		catch( OptionProcessorMutuallyExclusiveOptionsSuppliedException const& e ) {
			ui().logger() << "!! Error (" << e.what() << "):\n" ;
			ui().logger() << "Options \"" << e.first_option()
			<< "\" and \"" << e.second_option()
			<< "\" cannot be supplied at the same time.\n"
			<< "Please consult the documentation, or use \""
			<< m_application_name << " " << m_options->get_help_option_name()
			<< "\" for more information.\n" ;
			throw HaltProgramWithReturnCode( -1 ) ;
		}
		catch( OptionProcessorImpliedOptionNotSuppliedException const& e ) {
			ui().logger() << "!! Error (" << e.what() << "):\n" ;
			ui().logger() << "When using option \"" << e.first_option()
			<< "\", option \"" << e.second_option()
			<< "\" must also be supplied.\n"
			<< "Please consult the documentation, or use \""
			<< m_application_name << " " << m_options->get_help_option_name()
			<< "\" for more information.\n" ;
			throw HaltProgramWithReturnCode( -1 ) ;
		}
		catch( OptionProcessingException const& e ) {
			ui().logger() << "!! Error (" << e.what() << "): " << e.message() << ".\n" ;
			ui().logger() << "Please consult the documentation, or use \""
			<< m_application_name << " " << m_options->get_help_option_name()
			<< "\" for more information.\n" ;
			throw HaltProgramWithReturnCode( 0 );
		}
		catch( OptionProcessorHelpRequestedException const& ) {
			write_start_banner() ;
			std::cout << "Usage: "
				<< m_application_name << " <options>\n"
				<< "\nOPTIONS:\n"
				<< *(m_options)
				<< "\n" ;
			throw HaltProgramWithReturnCode( 0 );
		}
		catch( std::runtime_error const& e ) {
			ui().logger() << "!! Error: " << e.what() << ".\n" ;
			throw HaltProgramWithReturnCode( -1 ) ;
		}
		catch( std::exception const& e ) {
			ui().logger() << "!! Error (" << e.what() << "): \n";
			ui().logger() << "Please use \""
			<< m_application_name << " " << m_options->get_help_option_name()
			<< "\" for more information.\n" ;
			throw HaltProgramWithReturnCode( -1 ) ;
		}
		
	}

	ApplicationContext::~ApplicationContext() {
		write_end_banner() ;
	}

	namespace {
		std::auto_ptr< std::ostream > open_file_for_output( std::string const& filename ) {
			std::ios_base::openmode open_mode = std::ios_base::out ;

			std::auto_ptr< std::ostream > result( new std::ofstream( filename, open_mode )) ;
			if( !(*result) ) {
				throw std::runtime_error( "open_file_for_output(): could not open file \"" + filename + "\"." ) ;
			}

			return std::auto_ptr< std::ostream >( result.release() ) ;
		}
	}

	void ApplicationContext::construct_logger( std::string const& log_option ) {
		std::string filename ;
		if( options().check_if_option_is_defined( log_option ) && options().check_if_option_has_value( log_option )) {
			std::string const& filename = options().get_value< std::string > ( log_option ) ;
			if( filename != "" ) {
				ui().logger().add_stream( "log", open_file_for_output( filename )) ;
			}
			else {
				ui().logger().add_stream( "log", std::unique_ptr< std::ostream >( new null_ostream )) ;
			}
		}
		else {
			ui().logger().add_stream( "log", std::unique_ptr< std::ostream >( new null_ostream )) ;
		}
	}

	OptionProcessor& ApplicationContext::options() const { return *m_options ; }

	ApplicationContext::UIContext& ApplicationContext::ui() const {
		return *m_ui ;
	}

	void ApplicationContext::write_start_banner() {
		m_ui->logger() << "\nWelcome to " << m_application_name << "\n" ;
		if( m_application_version != "" ) {
			m_ui->logger() << "(version: " << m_application_version << ")\n" ;
		}
		m_ui->logger() << "\n(C) 2009-2017 University of Oxford\n\n";
	}

	void ApplicationContext::write_end_banner() {
		m_ui->logger() << "\n"
			<< "Thank you for using " << m_application_name << ".\n" ;
	}
	
	std::string const& ApplicationContext::application_name() const {
		return m_application_name ;
	}
	
}
