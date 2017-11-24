/*
 *  Copyright (C) 2012  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef __VCF_GENOTYPE_FORMAT_H__
#define __VCF_GENOTYPE_FORMAT_H__

#include "VcfGenotypeField.h"

/// This header file provides interface to read/write VCF files.
class VcfGenotypeFormat : public VcfGenotypeField
{
public:
    static const int GENOTYPE_INDEX_NA = -1;

    /// Default Constructor, initializes the variables.
    VcfGenotypeFormat();

    /// Destructor
    virtual ~VcfGenotypeFormat();
    
    /// Read the genotype format from the file up until the next \t,\n, or EOF.
    /// \param filePtr IFILE to read from.
    /// \return true if a tab ended the field, false if it was \n or EOF.
    bool read(IFILE filePtr);

    /// Get the index of the specified key.
    /// \param key to find the index for.
    /// \return index of the specified key or GENOTYPE_INDEX_NA if the key
    /// is not found.
    int getIndex(const std::string& key);

    /// Get the GT index, returns GENOTYPE_INDEX_NA if it is not found..
    inline int getGTIndex() { return(myGTIndex); }

    /// Return true if the specified subField should be read/stored.
    bool storeIndex(unsigned int index);

    /// Get Original number of fields - this is different than the number
    /// of stored fields.
    int getOrigNumFields() { return(myStoreIndices.size()); }

protected:
    /// reset the sample for a new entry.
    virtual void internal_reset();

private:
    VcfGenotypeFormat(const VcfGenotypeFormat& vcfGenotypeFormat);
    VcfGenotypeFormat& operator=(const VcfGenotypeFormat& vcfGenotypeFormat);

    int myGTIndex;

    // Set when reading the format by checking the readFields in
    // VcfRecordGenotype and queried when reading the samples fields.
    std::vector<bool> myStoreIndices;
};

#endif
