
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <memory>
#include <string>
#include "appcontext/OstreamTee.hpp"
#include "appcontext/progress_bar.hpp"
#include "appcontext/CmdLineUIContext.hpp"

namespace appcontext {
	CmdLineUIContext::CmdLineUIContext()
		: m_logger( new OstreamTee() )
	{
	}

	CmdLineUIContext::~CmdLineUIContext()
	{
		for(
			std::map< std::string, ProgressContextImpl* >::iterator i = m_progress_contexts.begin();
			i != m_progress_contexts.end();
			++i
		) {
			delete i->second ;
		}
	}

	OstreamTee& CmdLineUIContext::logger() const {
		return *m_logger ;
	}

	ProgressContextProxy CmdLineUIContext::get_progress_context( std::string const& name, std::string const& type ) {
		assert( m_progress_contexts.size() == 0 ) ;
		std::map< std::string, ProgressContextImpl* >::iterator where = m_progress_contexts.find( name ) ;
		assert( where == m_progress_contexts.end() ) ;

		if( type == "bar" ) {
			m_progress_contexts[ name ] = new ProgressBarProgressContext( *this, name ) ;
			where = m_progress_contexts.find( name ) ;
		}
		else {
			assert(0) ;
		}
		return ProgressContextProxy( *this, *(where->second) ) ;
	}

	void CmdLineUIContext::remove_progress_context_impl( std::string const& name ) {
		std::map< std::string, ProgressContextImpl* >::iterator where = m_progress_contexts.find( name ) ;
		assert( where != m_progress_contexts.end() ) ;
		delete where->second ;
		m_progress_contexts.erase( where ) ;
	}

	ProgressBarProgressContext::ProgressBarProgressContext( CmdLineUIContext const& ui_context, std::string const& name ):
		m_ui_context( ui_context ),
		m_name( name ),
		m_last_time( m_timer.elapsed() ),
		m_last_count( 0 )
	{}

	void ProgressBarProgressContext::notify_progress(
		std::size_t const count,
		boost::optional< std::size_t > const total_count
	) const {
		double time_now = m_timer.elapsed() ;
		if((count == 0) || (total_count && count == *total_count) || (time_now - m_last_time) > 1 ) {
			print_progress( count, total_count, m_name, 60 ) ;
			m_last_time = time_now ;
		}
		m_last_count = count ;
		m_last_total_count = total_count ;
	}

	void ProgressBarProgressContext::notify_progress() const {
		notify_progress( m_last_count + 1, m_last_total_count ) ;
	}

	void ProgressBarProgressContext::finish() const {
		print_progress( m_last_count, m_last_total_count, m_name, 60 ) ;
		m_ui_context.logger() << "\n" ;
	}

	void ProgressBarProgressContext::print_progress(
		std::size_t const count,
		boost::optional< std::size_t > const total_count,
		std::string const& msg,
		std::size_t const max_msg_length
	) const {

		m_ui_context.logger()["screen"] << "\r" ;
	
		if( msg != "" ) {
			if( max_msg_length >= 3 && msg.size() > max_msg_length ) {
				m_ui_context.logger() << std::setw( max_msg_length - 3 ) << std::left << msg.substr( 0, max_msg_length - 3 ) << "...: " ;
			}
			else {
				m_ui_context.logger() << std::setw( max_msg_length ) << std::left << msg.substr( 0, max_msg_length ) << ": " ;
			}
		}
	
		if( total_count ) {	
			double progress = (total_count == 0) ? 1.0 : (static_cast< double >( count ) / *total_count ) ;
			
			if( count == *total_count ) {
				m_ui_context.logger()
					<< get_progress_bar( 30, progress )
					<< " (" << count << "/" << *total_count
					<< "," << m_timer.display()
					<< "," << std::fixed << std::setprecision(1) << (static_cast< double >( count ) / m_timer.elapsed()) << "/s"
					<< ")" << std::flush ;
			}
			else {
				m_ui_context.logger()["screen"]
					<< get_progress_bar( 30, progress )
					<< " (" << count << "/" << *total_count
					<< "," << m_timer.display() ;
				if( m_timer.elapsed() > 0.0 ) {
					m_ui_context.logger()["screen"]
						<< "," << std::fixed << std::setprecision(1) << (static_cast< double >( count ) / m_timer.elapsed()) << "/s" ;
				}
				m_ui_context.logger()["screen"] << ")" << std::flush ;
			}
		}
		else {
			m_ui_context.logger()["screen"]
				<< " (" << count << "/" << "?" << "," << m_timer.display() ;
			if( m_timer.elapsed() > 0.0 ) {
				m_ui_context.logger()["screen"]
					<< "," << std::fixed << std::setprecision(1) << (static_cast< double >( count ) / m_timer.elapsed()) << "/s" ;
			}
			m_ui_context.logger()["screen"] << ")" << std::flush ;
		}
	}
}
