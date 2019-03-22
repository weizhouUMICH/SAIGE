
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <memory>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <boost/tuple/tuple.hpp>
#include <boost/format.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include "genfile/bgen/bgen.hpp"
#include "genfile/bgen/IndexQuery.hpp"
#include "genfile/bgen/View.hpp"

namespace genfile {
	namespace bgen {
		View::UniquePtr View::create( std::string const& filename ) {
			return View::UniquePtr( new View( filename )) ;
		}

		/* View implementation */
		View::View( std::string const& filename ):
			m_filename( filename ),
			m_variant_i(0),
			m_have_sample_ids( false ),
			m_state( e_NotOpen )
		{
			setup( m_filename ) ;
			m_file_position = m_stream->tellg() ;
		}

		std::size_t View::number_of_samples() const {
			return m_context.number_of_samples ;
		}

		void View::set_query( IndexQuery::UniquePtr query ) {
			m_index_query = query ;
			if( m_index_query->number_of_variants() > 0 ) {
				m_stream->seekg( m_index_query->locate_variant(0).first ) ;
			}
			m_file_position = m_stream->tellg() ;
		}

		View::FileMetadata const& View::file_metadata() const {
			return m_file_metadata ;
		}

		genfile::bgen::Context const& View::context() const {
			return m_context ;
		}

		std::streampos View::current_file_position() const {
			return m_file_position ;
		}

		uint32_t View::number_of_variants() const {
			if( m_index_query.get() ) {
				return m_index_query->number_of_variants() ;
			} else {
				return m_context.number_of_variants ;
			}
		}

		std::ostream& View::summarise( std::ostream& o ) const {
			o << "View: bgen file ("
				<< ( m_context.flags & genfile::bgen::e_Layout2 ? "layout = 2" : "layout = 1" )
				<< ", " ;
			std::string compression = "no" ;
			if( (m_context.flags & genfile::bgen::e_CompressedSNPBlocks) == genfile::bgen::e_ZlibCompression ) {
				compression = "zlib" ;
			} else if( (m_context.flags & genfile::bgen::e_CompressedSNPBlocks) == genfile::bgen::e_ZstdCompression ) {
				compression = "zstd" ;
			}
			o << compression << " compression)" ;
			o << " with " 
				<< m_context.number_of_samples << " " << ( m_have_sample_ids ? "named" : "anonymous" ) << " samples and "
				<< m_context.number_of_variants << " variants.\n" ;
			if( m_index_query.get() ) {
				o << "IndexQuery: query will return " << m_index_query->number_of_variants() << " variants.\n" ;
			}
			return o ;
		}

		namespace {
			// callbacks to avoid using C++11 lambda functions
			void resize_vector( std::vector< std::string >* target, std::size_t n ) {
				target->resize(n) ;
			}
			
			void set_vector_element( std::vector< std::string >* target, std::size_t i, std::string const& value ) {
				target->at(i) = value ;
			}

			void push_back_vector( std::vector< std::string >* target, std::string const& value ) {
				target->push_back( value ) ;
			}
		}

		bool View::read_variant(
			std::string* SNPID,
			std::string* rsid,
			std::string* chromosome,
			uint32_t* position,
			std::vector< std::string >* alleles
		) {
			assert( m_state == e_ReadyForVariant ) ;

			if( m_index_query.get() ) {
				if( m_variant_i == m_index_query->number_of_variants() ) {
					return false ;
				}
				IndexQuery::FileRange const range = m_index_query->locate_variant( m_variant_i ) ;
				m_stream->seekg( range.first ) ;
			}

			if(
				genfile::bgen::read_snp_identifying_data(
					*m_stream, m_context,
					SNPID, rsid, chromosome, position,
//					[&alleles]( std::size_t n ) { alleles->resize( n ) ; },
					boost::bind( &resize_vector, alleles, _1 ),
//					[&alleles]( std::size_t i, std::string const& allele ) { alleles->at(i) = allele ; }
					boost::bind( &set_vector_element, alleles, _1, _2 )
				)
			) {
				m_state = e_ReadyForProbs ;
				return true ;
			} else {
				return false ;
			}
		}

		// Read and uncompress genotype probability data, and unpack
		// it into constituent parts, but don't do a full parse.
		// This can lead to more efficient code paths than a full parse for some operations.
		//
		// Currently this function works for 'layout=2' files, e.g. v1.2 and above only.
		// Data is returned in the fields of the supplied 'pack' object, which
		// is defined in bgen.hpp.
		void View::read_and_unpack_v12_genotype_data_block(
			genfile::bgen::v12::GenotypeDataBlock* pack
		) {
			assert( (m_context.flags & genfile::bgen::e_Layout) == genfile::bgen::e_Layout2 ) ;
			std::vector< byte_t > const& buffer = read_and_uncompress_genotype_data_block() ;
			pack->initialise( m_context, &buffer[0], &buffer[0] + buffer.size() ) ;
			++m_variant_i ;
		}

		// Ignore genotype probability data for the SNP just read using read_variant()
		// After calling this method it should be safe to call read_variant()
		// to fetch the next variant from the file.
		void View::ignore_genotype_data_block() {
			assert( m_state == e_ReadyForProbs ) ;
			genfile::bgen::ignore_genotype_data_block( *m_stream, m_context ) ;
			m_file_position = m_stream->tellg() ;
			m_state = e_ReadyForVariant ;
			++m_variant_i ;
		}

		// Open the bgen file, read header data and gather metadata.
		void View::setup( std::string const& filename ) {
			m_file_metadata.filename = filename ;
			m_file_metadata.last_write_time = boost::filesystem::last_write_time( filename ) ;

			// Open the stream
			m_stream.reset(
				new std::ifstream( filename.c_str(), std::ifstream::binary )
			) ;
			if( !*m_stream ) {
				throw std::invalid_argument( filename ) ;
			}

			// get file size
			{
				std::ios::streampos origin = m_stream->tellg() ;
				m_stream->seekg( 0, std::ios::end ) ;
				m_file_metadata.size = m_stream->tellg() - origin ;
				m_stream->seekg( 0, std::ios::beg ) ;
			}
			// read first (up to) 1000 bytes.
			{
				m_file_metadata.first_bytes.resize( 1000, 0 ) ;
				m_stream->read( reinterpret_cast< char* >( &m_file_metadata.first_bytes[0] ), 1000 ) ;
				m_file_metadata.first_bytes.resize( m_stream->gcount() ) ;
				m_stream->clear() ;
			}

			m_state = e_Open ;

			// Read the offset, header, and sample IDs if present.
			m_stream->seekg( 0, std::ios::beg ) ;
			genfile::bgen::read_offset( *m_stream, &m_offset ) ;
			genfile::bgen::read_header_block( *m_stream, &m_context ) ;

			if( m_context.flags & genfile::bgen::e_SampleIdentifiers ) {
				genfile::bgen::read_sample_identifier_block(
					*m_stream, m_context,
					boost::bind( &push_back_vector, &m_sample_ids, _1 )
					//[this]( std::string id ) { m_sample_ids.push_back( id ) ; }
				) ;
				m_have_sample_ids = true ;
			}
	
			// read data up to first data block.
			m_postheader_data.resize( m_offset+4 - m_stream->tellg() ) ;
			m_stream->read( reinterpret_cast< char* >( &m_postheader_data[0] ), m_postheader_data.size() ) ;
			if( m_stream->gcount() != m_postheader_data.size() ) {
				throw std::invalid_argument(
					(
						boost::format(
							"BGEN file (\"%s\") appears malformed - offset specifies more bytes (%d) than are in the file."
						) % filename % m_offset
					).str()
				) ;
			}

			// Jump to the first variant data block.
			// m_stream->seekg( m_offset + 4 ) ;

			// We keep track of state (though it's not really needed for this implementation.)
			m_state = e_ReadyForVariant ;
		}

		// Utility function to read and uncompress variant genotype probability data
		// without further processing.
		std::vector< byte_t > const& View::read_and_uncompress_genotype_data_block() {
			assert( m_state == e_ReadyForProbs ) ;
			genfile::bgen::read_genotype_data_block( *m_stream, m_context, &m_buffer1 ) ;
			m_file_position = m_stream->tellg() ;
			m_state = e_ReadyForVariant ;
			genfile::bgen::uncompress_probability_data( m_context, m_buffer1, &m_buffer2 ) ;
			++m_variant_i;
			return m_buffer2 ;
		}
	}
}
