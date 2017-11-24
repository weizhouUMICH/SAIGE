/*
 *  Copyright (C) 2011-2012  Regents of the University of Michigan,
 *                      Hyun Min Kang, Matthew Flickenger, Matthew Snyder,
 *                      and Goncalo Abecasis
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


#ifndef __VCF_RECORD_GENOTYPE_H__
#define __VCF_RECORD_GENOTYPE_H__

#include <list>
#include "VcfRecordField.h"
#include "ReusableVector.h"
#include "VcfSubsetSamples.h"
#include "VcfGenotypeFormat.h"
#include "VcfGenotypeSample.h"

/// This header file provides interface to read/write VCF files.
class VcfRecordGenotype : public VcfRecordField
{
public:
    /// When reading, store all fields.
    static void storeAllFields();

    /// When reading, store the specified field in addition to any others
    /// that have previously been specified.  By default, all fields are stored.
    /// Use this method if you only want to store certain fields.
    static void addStoreField(const char* field);

    /// Return true if the specified field has been set to be stored.
    static bool storeField(std::string& field);

    /// Default Constructor, initializes the variables.
    VcfRecordGenotype();

    /// Destructor
    virtual ~VcfRecordGenotype();
    
    /// Read this genotype field from the file up until the next \t,\n, or EOF.
    /// \param filePtr IFILE to read from.
    /// \return true if a tab ended the field, false if it was \n or EOF (always
    /// returns false since this is the last field on the line).
    bool read(IFILE filePtr);

    /// Read this genotype field from the file up until the next \t,\n, or EOF.
    /// \param filePtr IFILE to read from.
    /// \param subsetInfo pointer to optional subsetting information.
    /// \return true if a tab ended the field, false if it was \n or EOF (always
    /// returns false since this is the last field on the line).
    bool read(IFILE filePtr, VcfSubsetSamples* subsetInfo);

    /// Write the genotype field to the file, without printing the
    // starting/trailing '\t'.
    /// \param filePtr IFILE to write to.
    /// \return true if the field was successfully written to the specified
    ///  filePtr, false if not.
    bool write(IFILE filePtr);

    /// reset the field for a new entry.
    void reset();
    /// reset the field for a new entry.
    void clear() {reset();}

    /// Get a pointer to the string containing the value associated with the
    /// specified key for the specified sample
    /// (the pointer will be invalid if the field is changed/reset).  
    /// \param key to find the value for.
    /// \param sampleNum which sample to get the value for (starts at 0).
    /// \return const pointer to the string value for this key, NULL if
    /// the sample or the key was not found.
    const std::string* getString(const std::string& key, 
                                 int sampleNum);

    /// Set the string associated with the specified key for the specified
    /// sample, returns true if set, false if not.
    bool setString(const std::string& key, int sampleNum, 
                   const std::string& value);

    int getGT(int sampleNum, unsigned int gtIndex);
    void setGT(int sampleNum, unsigned int gtIndex, int newGt);

    int getNumGTs(int sampleNum);

    /// Determine if all samples are phased.  Returns true if all are phased
    /// and false if any are not phased or if any are unphased.
    bool allPhased();

    /// Determine if all samples are unphased.  Returns true if all are unphased
    /// and false if any are not unphased or if any are phased.
    bool allUnphased();

    /// Determine if all samples have all the genotype alleles specified. 
    /// Returns true if all genotype alleles are specified and false if
    /// any are missing ('.') or if GT is not specified.
    bool hasAllGenotypeAlleles();

    /// Determine if the specified sample number is phased, returns true
    /// if it is phased and false if it is unphased or the sample number is
    /// invalid.
    bool isPhased(int sampleNum);

    /// Determine if the specified sample number is unphased, returns true
    /// if it is unphased and false if it is phased or the sample number is
    /// invalid.
    bool isUnphased(int sampleNum);

    /// Determine if the specified sample number has all of its genotype
    /// alleles specified.  Returns true if all genotype alleles are specified
    /// and false if any are missing ('.') or if GT is not specified.
    bool hasAllGenotypeAlleles(int sampleNum);

    /// Get the number of samples.
    inline int getNumSamples() const { return(mySamples.size()); }

protected:

private:
    VcfRecordGenotype(const VcfRecordGenotype& gt);
    VcfRecordGenotype& operator=(const VcfRecordGenotype& gt);

    // Fields that should be stored when reading for all records.
    static std::set<std::string> ourStoreFields;

    VcfGenotypeFormat myFormat;
    ReusableVector<VcfGenotypeSample> mySamples;
};

#endif
