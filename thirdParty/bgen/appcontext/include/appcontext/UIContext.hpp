
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef UICONTEXT_UICONTEXT_HPP
#define UICONTEXT_UICONTEXT_HPP

#include <string>
#include <boost/optional.hpp>
#include "appcontext/OstreamTee.hpp"

namespace appcontext {
	struct ProgressContextImpl
	{
		virtual void notify_progress(
			std::size_t const count,
			boost::optional< std::size_t > const total_count
		) const = 0 ;
		virtual void notify_progress() const = 0 ;
		virtual void finish() const = 0 ;
		virtual std::string name() const = 0 ;
		virtual void restart_timer() const = 0 ;
		virtual ~ProgressContextImpl() {} ;
	protected:
		friend struct ProgressContextProxy ;
	} ;

	struct UIContext ;

	struct ProgressContextProxy
	{
		ProgressContextProxy( UIContext&, ProgressContextImpl const& progress_context ) ;
		ProgressContextProxy( ProgressContextProxy const& other ) ;
		ProgressContextProxy& operator=( ProgressContextProxy const& other ) ;
		~ProgressContextProxy() ;

		void operator()(
			std::size_t const count,
			boost::optional< std::size_t > const total_count
		) const {
			notify_progress( count, total_count ) ;
		}

		void notify_progress(
			std::size_t const count,
			boost::optional< std::size_t > const total_count = boost::optional< std::size_t >()
		) const ;

		void notify_progress() const ;

		void finish() const ;

		void restart_timer() const ;
		
	private:
		UIContext* m_ui_context ;
		mutable ProgressContextImpl const* m_progress_context ;
	} ;

	struct UIContext
	{
		typedef ProgressContextProxy ProgressContext ;
	
		virtual ~UIContext() {}
		virtual OstreamTee& logger() const = 0 ;
		virtual ProgressContext get_progress_context( std::string const& name = "", std::string const& type = "bar" ) = 0 ;
	private:
		friend struct ProgressContextProxy ;
		void remove_progress_context( std::string const& name ) {
			remove_progress_context_impl( name ) ;
		}
		virtual void remove_progress_context_impl( std::string const& name ) = 0 ;
	} ;
}
#endif
