
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef BGEN_INDEX_QUERY_HPP
#define BGEN_INDEX_QUERY_HPP

#include <boost/function.hpp>
#include <boost/optional.hpp>
#include <stdint.h>
#include <vector>
#include <string>
#include <ctime>
#include "db/sqlite3.hpp"

namespace genfile {
	namespace bgen {
		// Base class representing a query against a BGEN file index
		struct IndexQuery {
		public:
			// We use std::auto_ptr to avoid using C++11 features here.
			typedef std::auto_ptr< IndexQuery > UniquePtr ;
			typedef uint8_t byte_t ;

			static UniquePtr create( std::string const& filename, std::string const& table_name = "Variant" ) ;

		public:
			struct FileMetadata ;
			typedef boost::optional< FileMetadata > OptionalFileMetadata ;
			struct GenomicRange {
				GenomicRange(): m_start(0), m_end(0) {}
				GenomicRange( GenomicRange const& other ):
					m_chromosome( other.m_chromosome ),
					m_start( other.m_start ),
					m_end( other.m_end )
				{}
					
				GenomicRange& operator=( GenomicRange const& other ) {
					m_chromosome = other.m_chromosome ;
					m_start = other.m_start ;
					m_end = other.m_end ;
					return *this ;
				}
				GenomicRange(
					std::string const& chromosome,
					uint32_t start,
					uint32_t end
				):
					m_chromosome( chromosome ),
					m_start( start ),
					m_end( end )
				{
					if( m_end < m_start ) {
						throw std::invalid_argument( "end" ) ;
					}
				}

				std::string const& chromosome() const { return m_chromosome ; }
				uint32_t const& start() const { return m_start ; }
				uint32_t const& end() const { return m_end ; }
				
			private:
				std::string m_chromosome ;
				uint32_t m_start ;
				uint32_t m_end ;
			} ;
			//typedef boost::tuple< std::string, uint32_t, uint32_t > GenomicRange ;
			typedef std::pair< int64_t, int64_t> FileRange ;
			typedef boost::function< void ( std::size_t n, boost::optional< std::size_t > total ) > ProgressCallback ;

		public:
			virtual ~IndexQuery() {} ;
			virtual OptionalFileMetadata const& file_metadata() const = 0 ;

			// Methods for building queries
			// Each method returns this object, allowing methods to be chained
			virtual IndexQuery& include_range( GenomicRange const& range ) = 0 ;
			virtual IndexQuery& exclude_range( GenomicRange const& range ) = 0 ;
			virtual IndexQuery& include_rsids( std::vector< std::string > const& ids ) = 0 ;
			virtual IndexQuery& exclude_rsids( std::vector< std::string > const& ids ) = 0 ;

			// Initialise must be called before calling number_of_variants() or locate_variant().
			virtual void initialise( ProgressCallback callback = ProgressCallback() ) = 0 ;
			// Report the number of variants in this query.
			virtual std::size_t number_of_variants() const = 0 ;
			// Report the number of variants in this query.
			virtual FileRange locate_variant( std::size_t index ) const = 0 ;

			struct FileMetadata {
				FileMetadata():
					size(-1)
				{}

				FileMetadata( FileMetadata const& other ):
					filename( other.filename ),
					size( other.size ),
					last_write_time( other.last_write_time ),
					first_bytes( other.first_bytes )
				{}

				FileMetadata& operator=( FileMetadata const& other ) {
					filename = other.filename ;
					size = other.size ;
					last_write_time = other.last_write_time ;
					first_bytes = other.first_bytes ;
					return *this ;
				}

				std::string filename ;
				int64_t size ;
				std::time_t last_write_time ;
				std::vector< byte_t > first_bytes ;
			} ;
		} ;
		
		// Class for index queries implemented using a sqlite file, a la bgenix.
		struct SqliteIndexQuery: public IndexQuery {
		public:
			// We use auto_ptr to avoid using C++11 features here.
			typedef std::auto_ptr< SqliteIndexQuery > UniquePtr ;

		public:
			// Construct given an index file and an index table name
			SqliteIndexQuery( std::string const& filename, std::string const& table_name = "Variant" ) ;

			// Methods for building queries
			// Each method returns this object, allowing methods to be chained

			// Include variants in a range
			SqliteIndexQuery& include_range( GenomicRange const& range ) ;
			// Exclude variants in a range
			SqliteIndexQuery& exclude_range( GenomicRange const& range ) ;
			// Include variants with one of the given rsids.  The list provided must be unique.
			SqliteIndexQuery& include_rsids( std::vector< std::string > const& ids ) ;
			// Exclude variants with one of the given rsids.  The list provided must be unique.
			SqliteIndexQuery& exclude_rsids( std::vector< std::string > const& ids ) ;

		public:
			// IndexQuery methods
			void initialise( ProgressCallback callback = ProgressCallback() ) ;
			OptionalFileMetadata const& file_metadata() const ;
			std::size_t number_of_variants() const ;
			FileRange locate_variant( std::size_t index ) const ;

		private:
			db::Connection::UniquePtr open_connection( std::string const& filename ) const ;
			OptionalFileMetadata load_metadata( db::Connection& connection ) const ;
			db::Connection::StatementPtr build_query() const ;
	
		private:
			db::Connection::UniquePtr m_connection ;
			OptionalFileMetadata const m_metadata ;
			std::string const m_index_table_name ;
			struct QueryParts {
				std::string join ;
				std::string inclusion ;
				std::string exclusion ;
			} ;
			QueryParts m_query_parts ;
			bool m_initialised ;
			std::vector< std::pair< int64_t, int64_t> > m_positions ;
		} ;
		
	}
}

#endif
