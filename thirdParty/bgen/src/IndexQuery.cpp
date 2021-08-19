
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <memory>
#include <boost/format.hpp>
#include <boost/optional.hpp>
#include "db/sqlite3.hpp"
#include "genfile/bgen/IndexQuery.hpp"

// #define DEBUG 1

namespace genfile {
	namespace bgen {
		IndexQuery::UniquePtr IndexQuery::create( std::string const& filename, std::string const& table_name ) {
			return IndexQuery::UniquePtr( new SqliteIndexQuery( filename, table_name )) ;
		}

		SqliteIndexQuery::SqliteIndexQuery( std::string const& filename, std::string const& table_name ):
			m_connection( open_connection( filename ) ),
			m_metadata( load_metadata( *m_connection ) ),
			m_index_table_name( table_name ),
			m_initialised( false )
		{
		}
	
		boost::optional< SqliteIndexQuery::FileMetadata > const&
		SqliteIndexQuery::file_metadata() const {
			return m_metadata ;
		}

		void SqliteIndexQuery::initialise( ProgressCallback callback ) {
			db::Connection::StatementPtr stmt = build_query() ;
			if( callback ) {
				callback( 0, boost::optional< std::size_t >() ) ;
			}
			m_positions.reserve( 1000000 ) ;
			std::size_t batch_i = 0 ;
			for( stmt->step() ; !stmt->empty(); stmt->step(), ++batch_i ) {
				int64_t const pos = stmt->get< int64_t >( 0 ) ;
				int64_t const size = stmt->get< int64_t >( 1 ) ;
				assert( pos >= 0 ) ;
				assert( size >= 0 ) ;
				m_positions.push_back( std::make_pair( int64_t( pos ), int64_t( size ))) ;
				if( callback ) {
					callback( m_positions.size(), boost::optional< std::size_t >() ) ;
				}
			}
	#if DEBUG
			std::cerr << "SqliteIndexQuery::initialise(): read positions for " << m_positions.size() << " variants.\n" ;
	#endif
	
			m_initialised = true ;
		}

		std::size_t SqliteIndexQuery::number_of_variants() const {
			assert( m_initialised ) ;
			return m_positions.size() ;
		}

		SqliteIndexQuery::FileRange SqliteIndexQuery::locate_variant( std::size_t index ) const {
	#if DEBUG
			std::cerr << "SqliteIndexQuery::locate_variant(" << index << ")...\n" ;
	#endif
			assert( m_initialised ) ;
			assert( index < m_positions.size() ) ;
			return m_positions[index] ;
		}

		SqliteIndexQuery& SqliteIndexQuery::include_range( GenomicRange const& range ) {
			m_query_parts.inclusion += ((m_query_parts.inclusion.size() > 0) ? " OR " : "" ) + (
				boost::format( "( chromosome == '%s' AND position BETWEEN %d AND %d )" ) % range.chromosome() % range.start() % range.end()
			).str() ;
			m_initialised = false ;
			return *this ;
		}

		SqliteIndexQuery& SqliteIndexQuery::exclude_range( GenomicRange const& range ) {
			m_query_parts.exclusion += ( m_query_parts.exclusion.size() > 0 ? " AND" : "" ) + (
				boost::format( " NOT ( chromosome == '%s' AND position BETWEEN %d AND %d )" )
					% range.chromosome() % range.start() % range.end()
			).str() ;
			m_initialised = false ;
			return *this ;
		}

		SqliteIndexQuery& SqliteIndexQuery::include_rsids( std::vector< std::string > const& ids ) {
			m_connection->run_statement( "CREATE TEMP TABLE IF NOT EXISTS tmpIncludedId( identifier TEXT NOT NULL PRIMARY KEY ) WITHOUT ROWID" ) ;
			db::Connection::StatementPtr insert_stmt = m_connection->get_statement( "INSERT INTO tmpIncludedId( identifier ) VALUES( ? )" ) ;
			for( std::size_t i = 0; i < ids.size(); ++i ) {
				insert_stmt->bind( 1, ids[i] ).step() ;
				insert_stmt->reset() ;
			}
			if( m_query_parts.join.find( "tmpIncludedId" ) == std::string::npos ) {
				m_query_parts.join += " LEFT OUTER JOIN tmpIncludedId TI ON TI.identifier == V.rsid" ;
				m_query_parts.inclusion += ( m_query_parts.inclusion.size() > 0 ? " OR" : "" ) + std::string( " TI.identifier IS NOT NULL" ) ;
			}
			m_initialised = false ;
			return *this ;
		}

		SqliteIndexQuery& SqliteIndexQuery::exclude_rsids( std::vector< std::string > const& ids ) {
			m_connection->run_statement( "CREATE TEMP TABLE IF NOT EXISTS tmpExcludedId( identifier TEXT NOT NULL PRIMARY KEY ) WITHOUT ROWID" ) ;
			db::Connection::StatementPtr insert_stmt = m_connection->get_statement( "INSERT INTO tmpExcludedId( identifier ) VALUES( ? )" ) ;
			for( std::size_t i = 0; i < ids.size(); ++i ) {
				insert_stmt->bind( 1, ids[i] ).step() ;
				insert_stmt->reset() ;
			}
			if( m_query_parts.join.find( "tmpExcludedId" ) == std::string::npos ) {
				m_query_parts.join += " LEFT OUTER JOIN tmpExcludedId TE ON TE.identifier == V.rsid" ;
				m_query_parts.exclusion += ( m_query_parts.exclusion.size() > 0 ? " AND" : "" ) + std::string( " TE.identifier IS NULL" ) ;
			}
			m_initialised = false ;
			return *this ;
		}

		db::Connection::UniquePtr SqliteIndexQuery::open_connection( std::string const& filename ) const {
			db::Connection::UniquePtr result ;
			try {
				result = db::Connection::create( "file:" + filename + "?nolock=1", "r" ) ;
				//result = db::Connection::create( filename, "r" ) ;
			} catch( db::ConnectionError const& e ) {
				throw std::invalid_argument( "Could not open the index file \"" + filename + "\"" ) ;
			}
			return result ;
		}

		SqliteIndexQuery::OptionalFileMetadata
		SqliteIndexQuery::load_metadata( db::Connection& connection ) const {
			OptionalFileMetadata result ;
			db::Connection::StatementPtr stmt = connection.get_statement( "SELECT * FROM sqlite_master WHERE name == 'Metadata' AND type == 'table'" ) ;
			stmt->step() ;
			if( !stmt->empty() ) {
				db::Connection::StatementPtr mdStmt = connection.get_statement( "SELECT filename, file_size, last_write_time, first_1000_bytes FROM Metadata" ) ;
				mdStmt->step() ;

				if( mdStmt->empty() ) {
					throw std::invalid_argument( "Index file appears malformed (empty \"Metadata\" table)" ) ;
				}
				FileMetadata metadata ;
				// Get metadata fields for comparison
				metadata.filename = mdStmt->get< std::string >( 0 ) ;
				metadata.size = mdStmt->get< int64_t >( 1 ) ;
				metadata.last_write_time = mdStmt->get< int64_t >( 2 ) ;
				metadata.first_bytes = mdStmt->get< std::vector< uint8_t > >( 3 ) ;
				
				result = metadata ;
			}
			return result ;
		}

		db::Connection::StatementPtr SqliteIndexQuery::build_query() const {
			std::string const select = "SELECT file_start_position, size_in_bytes FROM `"
				+ m_index_table_name + "` V" ;
			std::string const inclusion = ( m_query_parts.inclusion.size() > 0 ) ? ("(" + m_query_parts.inclusion + ")") : "" ;
			std::string const exclusion = ((m_query_parts.inclusion.size() > 0 && m_query_parts.exclusion.size() > 0 ) ? "AND " : "" )
				+ (( m_query_parts.exclusion.size() > 0 ) ? ("(" + m_query_parts.exclusion + ")") : "" ) ;
			std::string const where = (inclusion.size() > 0 || exclusion.size() > 0) ? ("WHERE " + inclusion + exclusion) : "" ;
			std::string const orderBy = "ORDER BY chromosome, position, rsid, allele1, allele2" ;
			std::string const select_sql = select + " " + m_query_parts.join + " " + where + " " + orderBy ;
	#if DEBUG
			std::cerr << "BgenIndex::build_query(): SQL is: \"" << select_sql << "\"...\n" ;
	#endif
	
			return m_connection->get_statement( select_sql ) ;
		}
	}
}

