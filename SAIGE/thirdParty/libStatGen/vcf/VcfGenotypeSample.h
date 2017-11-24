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


#ifndef __VCF_GENOTYPE_SAMPLE_H__
#define __VCF_GENOTYPE_SAMPLE_H__

#include "VcfGenotypeFormat.h"

/// This header file provides interface to read/write VCF files.
class VcfGenotypeSample : public VcfGenotypeField
{
public:
    static const int INVALID_GT;
    static const int MISSING_GT;
    static const std::string MISSING_FIELD;

    /// Default Constructor, initializes the variables.
    VcfGenotypeSample();

    /// Destructor
    virtual ~VcfGenotypeSample();
    
    /// Read this sample from the file up until the next \t,\n, or EOF.
    /// \param filePtr IFILE to read from.
    /// \param format the VCF Genotype Format field description.
    /// \return true if a tab ended the field, false if it was \n or EOF.
    bool read(IFILE filePtr, VcfGenotypeFormat& format);

    virtual bool write(IFILE filePtr);

    /// Get a pointer to the string containing the value associated with the
    /// specified key(the pointer will be invalid if the field is
    /// changed/reset).  
    /// \param key to find the falue for.
    /// \return const pointer to the string value for this key, NULL if
    /// the sample or the key was not found.
    const std::string* getString(const std::string& key);

    /// Set the string associated with the specified key, returns true if set,
    /// false if not.
    bool setString(const std::string& key, const std::string& value);

    inline bool isPhased()              { return(myPhased); }
    inline bool isUnphased()            { return(myUnphased); }
    inline bool hasAllGenotypeAlleles() { return(myHasAllGenotypeAlleles); }

    /// Return the integer allele at the specified index of the GT field.
    /// For example, a GT of 0/3, getGT(1) returns 3 and getGT(0) returns 0.
    /// Returns INVALID_GT if there is no GT at the specified index or GT was
    /// not specified and returns MISSING_GT if it is '.'.
    int getGT(unsigned int index);

    /// Set the integer allele at the specified index of the GT field.
    /// Requires the GT index to already exist.
    void setGT(unsigned int index, int newGt);

    /// Return the number of GT fields for this sample.
    int getNumGTs();

protected:
    /// reset the sample for a new entry.
    virtual void internal_reset();

    bool parseGT();

private:
    VcfGenotypeSample(const VcfGenotypeSample& vcfGenotypeSample);
    VcfGenotypeSample& operator=(const VcfGenotypeSample& vcfGenotypeSample);

    void updateGTString();

    VcfGenotypeFormat* myFormatPtr;

    bool myPhased;
    bool myUnphased;
    bool myHasAllGenotypeAlleles;
    bool myNewGT;

    std::vector<int> myGTs;
};

#endif
