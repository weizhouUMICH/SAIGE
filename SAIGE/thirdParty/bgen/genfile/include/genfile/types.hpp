//          Copyright Gavin Band 2008 - 2012.
//          Distributed under the Boost Software License, Version 1.0.
//          (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)
//
#ifndef GENFILE_TYPES_HPP
#define GENFILE_TYPES_HPP

#include <stdint.h>

namespace genfile {
	enum OrderType {
		eUnknownOrderType = 0,
		eUnorderedList = 1, 				// a list, treated as unordered
		eOrderedList = 2,					// a list, treated as ordered
		ePerUnorderedGenotype = 3, 			// E.g. genotype probabilities; GEN, BGEN v1.1, etc.
		ePerOrderedHaplotype = 4,			// E.g. Phased GT, SHAPEIT or IMPUTE haplotypes
		ePerUnorderedHaplotype = 5,			// E.g. Unphased GT, binary PED file.
		ePerPhasedHaplotypePerAllele = 6,	// E.g. BGEN v1.2-style haplotype probabilities.
		ePerAllele = 7,						// E.g. assay intensities.
		ePerSample = 8					// 'B' allele dosage, one value per sample.
	} ;
	enum ValueType {
		eUnknownValueType = 0,
		eProbability = 1,
		eAlleleIndex = 2,
		eDosage = 3
	} ;

	typedef uint8_t byte_t ;
}

#endif

