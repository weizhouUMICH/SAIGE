
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef DB_SQL_ERROR_HPP
#define DB_SQL_ERROR_HPP

namespace db {
	struct Error: public std::exception
	{
		Error( std::string const& caller, std::string const& db_spec, int error, std::string const& sql = "(unknown)" ) ;
		~Error() throw() ;
		char const* what() const throw() { return "db::Error" ; }

		int const& error_code() const { return m_error ; }
		std::string const& spec() const { return m_spec ; }
		std::string description() const ;
		std::string sql() const { return m_sql ; }
	private:
		std::string const m_caller ;
		std::string const m_spec ;
		int const m_error ;
		std::string const m_sql ;
	} ;
	
	struct ConnectionError: public Error
	{
		ConnectionError( std::string const& caller, std::string const& db_spec, int error_code ): Error( caller, db_spec, error_code ) {}
		char const* what() const throw() { return "db::ConnectionError" ; }
	} ;

	struct TransactionError: public Error
	{
		TransactionError( std::string const& caller, std::string const& db_spec, int error_code, std::string const& sql = "(unknown)" ): Error( caller, db_spec, error_code, sql ) {}
		char const* what() const throw() { return "db::TransactionError" ; }
	} ;

	struct StatementPreparationError: public Error
	{
		StatementPreparationError( std::string const& caller, std::string const& db_spec, int error_code, std::string const& sql = "(unknown)" ): Error( caller, db_spec, error_code, sql ) {}
		~StatementPreparationError() throw() {}
		char const* what() const throw() { return "db::StatementPreparationError" ; }
	} ;

	struct StatementStepError: public Error
	{
		StatementStepError( std::string const& caller, std::string const& db_spec, int error_code, std::string const& sql = "(unknown)" ): Error( caller, db_spec, error_code, sql ) {}
		char const* what() const throw() { return "db::StatementStepError" ; }
	} ;

	struct ValueBindError: public Error
	{
		~ValueBindError() throw() {}
		ValueBindError(
			std::string const& caller,
			std::string const& db_spec,
			int error_code,
			std::string const& slot
		): Error( caller, db_spec, error_code ), m_slot( slot ) {}
		char const* what() const throw() { return "db::StatementStepError" ; }
		std::string const& slot() const { return m_slot ; }
	private:
		std::string const m_slot ;
	} ;
	
}

#endif
