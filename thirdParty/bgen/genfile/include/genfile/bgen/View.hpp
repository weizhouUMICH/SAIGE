
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_BGEN_READVIEW_HPP
#define GENFILE_BGEN_READVIEW_HPP

#include <memory>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include "genfile/bgen/bgen.hpp"
#include "genfile/bgen/IndexQuery.hpp"

namespace {
	std::string to_string( std::size_t i ) {
		std::stringstream s ;
		s << i ;
		return s.str() ;
	}
}

namespace genfile {
	namespace bgen {
		struct View {
		public:
			typedef std::auto_ptr< View > UniquePtr ;
			typedef genfile::bgen::IndexQuery IndexQuery ;
			typedef genfile::bgen::IndexQuery::FileMetadata FileMetadata ;

			static UniquePtr create( std::string const& filename ) ;

		public:
			View( std::string const& filename ) ;

			// Restrict this reader to a set of variants specified by the given query
			void set_query( IndexQuery::UniquePtr query ) ;

			// Report high-level information about the file
			uint32_t number_of_variants() const ;
			std::size_t number_of_samples() const ;
			std::ostream& summarise( std::ostream& o ) const ;

			// Report the sample IDs in the file using the given setter object
			// (If there are no sample IDs in the file, report a dummy identifier).
			// Setter object must be callable as setter( index of sample, sample identifier ).
			template< typename Setter >
			void get_sample_ids( Setter setter ) const {
				if( m_have_sample_ids ) {
					for( std::size_t i = 0; i < m_context.number_of_samples; ++i ) {
						setter( m_sample_ids[i] ) ;
					}
				} else {
					for( std::size_t i = 0; i < m_context.number_of_samples; ++i ) {
						setter( "(anonymous_sample_" + to_string( i+1 ) + ")" ) ;
					}
				}
			}

			// Report low-level information about the file.
			genfile::bgen::Context const& context() const ;
			FileMetadata const& file_metadata() const ;
			std::streampos current_file_position() const ;

			// Attempt to read identifying information about the next available variant from the
			// returning data in the given fields.
			// If this method returns true, data was successfully read, and it should be safe to call
			// read_genotype_data_block(), ignore_genotype_data_block(), or read_and_unpack_v12_genotype_data_block().
			// If this method returns false, data could not be read indicating end of the file.
			bool read_variant(
				std::string* SNPID,
				std::string* rsid,
				std::string* chromosome,
				uint32_t* position,
				std::vector< std::string >* alleles
			) ;

			// Read, uncompress, and parse genotype probability data for the variant just read by read_variant().
			// Data is returned via a setter object, using the parse_genotype_data API documented on the wiki.
			// An example using this API is found in the bgen_to_vcf.cpp example program.
			template< typename ProbSetter >
			void read_genotype_data_block( ProbSetter& setter ) {
				assert( m_state == e_ReadyForProbs ) ;
				genfile::bgen::read_and_parse_genotype_data_block< ProbSetter >(
					*m_stream,
					m_context,
					setter,
					&m_buffer1,
					&m_buffer2
				) ;
				m_file_position = m_stream->tellg() ;
				m_state = e_ReadyForVariant ;
				++m_variant_i ;
			}
	
			// Read, uncompress, and unpack data for the variant just read by read_variant()
			// without doing a full parse of all the probability data.
			// This method is useful for applications that want to manipulate the bgen-formatted
			// data directly, e.g. for efficiency reasons.
			// Currently this function works for 'layout=2' files, e.g. v1.2 and above only.
			// The function will assert() if the data is not in this format.
			// Data is returned in the fields of the supplied 'pack' object.  See bgen.hpp for the
			// declaration of this object.
			void read_and_unpack_v12_genotype_data_block(
				genfile::bgen::v12::GenotypeDataBlock* pack
			) ;

			// Skip over (i.e. ignore) genotype probability data for the current variant.
			void ignore_genotype_data_block() ;

		private:
			// Open the bgen file, read header data and gather metadata.
			void setup( std::string const& filename ) ;
		public:
			// Utility function to read and uncompress variant genotype probability data
			// without further processing.
			std::vector< byte_t > const& read_and_uncompress_genotype_data_block() ;

		private:
			std::string const m_filename ;
			std::auto_ptr< std::istream > m_stream ;
			std::size_t m_variant_i ;
			IndexQuery::UniquePtr m_index_query ;

			// meta data used to avoid stale index files.
			FileMetadata m_file_metadata ;

			// offset byte from top of bgen file.
			uint32_t m_offset ;

			// bgen::Context object holds information from the header block,
			// including bgen flags
			genfile::bgen::Context m_context ;
	
			bool m_have_sample_ids ;
			std::vector< std::string > m_sample_ids ;

			// All data following header up to the first variant data block.
			std::vector< byte_t > m_postheader_data ;

			// We keep track of our state in the file.
			// This is not strictly necessary for this implentation but makes it clear that
			// the sequence of calls must be read_variant() followed by
			// ignore_genotype_data_block() or ignore_genotype_data_block() repeatedly.
			enum State { e_NotOpen = 0, e_Open = 1, e_ReadyForVariant = 2, e_ReadyForProbs = 3, eComplete = 4 } ;
			State m_state ;
	
			// To avoid issues with tellg() and failbit, we store the stream position at suitable points
			std::streampos m_file_position ;
	
			// Two buffers for processing
			std::vector< byte_t > m_buffer1 ;
			std::vector< byte_t > m_buffer2 ;
		} ;
	}
}

#endif
