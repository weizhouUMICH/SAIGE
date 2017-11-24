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

#ifndef __VCF_SUBSET_SAMPLES_H__
#define __VCF_SUBSET_SAMPLES_H__

#include <vector>
#include <set>
#include <string>
#include "VcfHeader.h"

class VcfSubsetSamples
{
public:
    VcfSubsetSamples()
        : mySampleSubsetIndicator(),
          mySampleNames()
    {}

    ~VcfSubsetSamples()
    {
        mySampleSubsetIndicator.clear();
    }

    void reset();

    /// Read the samples from the header initiallizing them all to be
    /// included/excluded based on the include paramater.  The header is
    /// not stored or updated based on any include/excludes.  The mapping
    /// between sampleNames & indexes is stored to be used for addIncludeSample
    /// and addExcludeSample.
    void init(const VcfHeader& header, bool include);

    /// Include the specified sample.  initSample must first be called to 
    /// specify the header mapping between index and sample.  
    /// \return true if the sample could be included, false if the sample
    /// was not found in the sample list so cannot be included.
    bool addIncludeSample(const char* sampleName);

    /// Exclude the specified sample.  initSample must first be called to 
    /// specify the header mapping between index and sample.  
    /// \return true if the sample was found in teh sample list and could be
    /// excluded, false if the sample was not found in the sample list.
    bool addExcludeSample(const char* sampleName);

    /// Initialize this object based on the sample names found in sampleFileName
    /// delimited by any of the characters in delims or '\n' and update the
    /// header to only include the specified samples.
    /// This also initializes this class to identify which samples should
    /// be kept/removed when reading records.
    bool init(VcfHeader& header, const char* sampleFileName, 
              const char* excludeSample, const char* excludeFileName,
              const char* delims = "\n");

    /// Return if the specified original sample index should be kept.
    /// This is only applicable after calling init.
    /// If the index is out of range, it will return false (do not keep).
    /// \param sampleIndex index into the original sample set to check if
    /// it should be kept.
    /// \return true if the sample index should be kept, false if not or if
    /// the index is out of range.
    bool keep(unsigned int sampleIndex);

private:
    VcfSubsetSamples(const VcfSubsetSamples& vcfSubsetSamples);
    VcfSubsetSamples& operator=(const VcfSubsetSamples& vcfSubsetSamples);

    // Read a list of samples from the specified file delimited by any of the
    // characters in delims or '\n' and store them in the specified container.
    bool readSamplesFromFile(const char* fileName, 
                             std::set<std::string>& sampleList,
                             const char* delims="\n");


    std::vector<bool> mySampleSubsetIndicator;

    // Used for initSample & addIncludeSample & addExcludeSample for
    // mapping between original sample names and indexes in
    // mySampleSubsetIndicator.
    std::vector<std::string>mySampleNames;
};

#endif
