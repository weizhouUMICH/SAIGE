
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef UICONTEXT_CMDLINEUICONTEXT_HPP
#define UICONTEXT_CMDLINEUICONTEXT_HPP

#include <memory>
#include <string>
#include <vector>
#include "appcontext/Timer.hpp"
#include "appcontext/OstreamTee.hpp"
#include "appcontext/UIContext.hpp"

namespace appcontext {
	struct CmdLineUIContext: public UIContext
	{
	public:
		CmdLineUIContext() ;
		~CmdLineUIContext() ;
	
		OstreamTee& logger() const ;
		ProgressContext get_progress_context( std::string const& name = "", std::string const& type = "bar" ) ;
	private:
		void remove_progress_context_impl( std::string const& name ) ;
	private:
		mutable std::map< std::string, ProgressContextImpl* > m_progress_contexts ;
		std::auto_ptr< OstreamTee > m_logger ;
	} ;


	struct ProgressBarProgressContext: public ProgressContextImpl
	{
	public:
		ProgressBarProgressContext( CmdLineUIContext const& ui_context, std::string const& name ) ;

		void notify_progress(
			std::size_t const count,
			boost::optional< std::size_t > const total_count
		) const ;

		void notify_progress() const ;

		void finish() const ;

		std::string name() const { return m_name ; }

		void restart_timer() const {
			m_timer.restart() ;
		}

	private:
		void print_progress(
			std::size_t const count,
			boost::optional< std::size_t > const total_count,
			std::string const& msg,
			std::size_t max_msg_length
		) const ;

		CmdLineUIContext const& m_ui_context ;
		std::string const m_name ;
		mutable Timer m_timer ;
		mutable double m_last_time ;
		
		mutable std::size_t m_last_count ;
		mutable boost::optional< std::size_t > m_last_total_count ;
	} ;
}

#endif
