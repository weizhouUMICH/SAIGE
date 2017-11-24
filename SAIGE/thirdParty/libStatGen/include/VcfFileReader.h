/*
 *  Copyright (C) 2010-2012  Regents of the University of Michigan,
 *                           Hyun Min Kang, Matthew Flickenger, Matthew Snyder,
 *                           and Goncalo Abecasis
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

#ifndef __VCF_FILE_READER_H__
#define __VCF_FILE_READER_H__

#include "VcfFile.h"
#include "VcfRecord.h"
#include "VcfRecordDiscardRules.h"
#include "VcfSubsetSamples.h"
#include "Tabix.h"

#ifdef __GXX_EXPERIMENTAL_CXX0X__
#include <unordered_set>
#else
#include <set>
#endif


/// This header file provides interface to read/write VCF files.
class VcfFileReader : public VcfFile
{
public:
    static const uint64_t DISCARD_NON_PHASED = 0x1;
    static const uint64_t DISCARD_MISSING_GT = 0x2;
    static const uint64_t DISCARD_FILTERED = 0x4;
    static const uint64_t DISCARD_MULTIPLE_ALTS = 0x8;

    /// Default Constructor, initializes the variables, but does not open
    /// any files.
    VcfFileReader();
    /// Destructor
    virtual ~VcfFileReader();

    /// Open the vcf file with the specified filename for reading.
    /// This method does no sample subsetting.
    /// \param  filename the vcf file to open for reading.
    /// \param header to be read from the file
    /// \return true = success; false = failure.
    virtual bool open(const char* filename, VcfHeader& header);

    /// Open the vcf file with the specified filename for reading
    /// subsetting the samples in the file to just the samples specified
    /// in the sample file.
    /// \param filename the vcf file to open for reading.
    /// \param header to be read from the file
    /// \param includeFileName file containing samples to keep
    ///                        or NULL if all samples except those specified as
    ///                        excluded should be included.
    /// \param excludeSample sample to be excluded or NULL
    /// \param excludeFileName file containing samples to remove or NULL
    ///                        or NULL if there is no file containing samples
    ///                        to exclude.                        
    /// \param delims deliminators separating the samples in the files ('\n'
    /// is always considered a delimiter even if it isn't specified).  When
    /// any of the specified delimiter characters is found in the file it
    /// indicates the end of a sample name.
    /// \return true = success; false = failure.
    virtual bool open(const char* filename, VcfHeader& header,
                      const char* includeFileName, const char* excludeSample,
                      const char* excludeFileName, const char* delims = "\n");

    /// Read the specified vcf index file.  It must be read prior to setting a
    /// read section, for seeking and reading portions of a vcf file.
    /// \param filename the name of the vcf index file to be read.
    /// \return true = success; false = failure.
    bool readVcfIndex(const char * filename);

    /// Read the bam index file using the VCF filename as a base. 
    /// It must be read prior to setting a read section, for seeking
    /// and reading portions of a vcf file.
    /// Must be read after opening the VCF file since it uses the
    /// VCF filename as a base name for the index file.
    /// First it tries filename.vcf.tbi. If that fails, it tries
    /// it without the .vcf extension, filename.tbi.
    /// \return true = success; false = failure.
    bool readVcfIndex();

    /// Get the VCF index (it must have already been read).
    /// \return a const pointer to the tabix object, or NULL if the
    /// index has not been read.
    const Tabix* getVcfIndex();

    /// Read the next Vcf record from the file until a line passes all
    /// discard rules (if any) or until the end of the file is found..
    /// \param record record to populate with the next record.
    /// \param subset ptr to subset of samples to keep.  This overrides
    ///         mySampleSubset that may have been set at open.
    /// \return true if successful, false if not.
    bool readRecord(VcfRecord& record, VcfSubsetSamples* subset = NULL);

    /// Only read the specified chromosome when readRecord is called.
    /// If an index is not used, the read section can only be set prior to 
    /// reading any records.
    bool setReadSection(const char* chromName);

    /// Only read the specified chromosome/positions when readRecord is called.
    /// If an index is not used, the read section can only be set prior to 
    /// reading any records.
    /// \param chromName chromosome name to read.
    /// \param start inclusive 1-based start positions of records that should be
    /// read for this chromosome.
    /// \param end exclusive 1-based end positions of records that should be
    /// read for this chromosome (this position is not read).
    /// \param overlap bool indicating whether or not records overlapping
    /// the specified region should be included even if they do not start
    /// in the region.  False (DEFAULT) means only read records that start
    /// in the region.  True means to read record's whose deletions extend 
    /// into the region.
    bool set1BasedReadSection(const char* chromName, 
                              int32_t start, int32_t end, 
                              bool overlap = false);

    /// Returns whether or not the end of the file has been reached.
    /// \return true = EOF; false = not eof.
    /// If the file is not open, true is returned.
    bool isEOF();

    /// Get the number of VCF records that have been processed (read/written)
    /// so far excluding any discarded records.
    int getNumKeptRecords() {return(myNumKeptRecords);}

    /// Get the number of VCF records that were read, even those outside
    /// a specified region and those discarded.
    int getTotalReadRecords() {return(myTotalRead);}

    /////////////////////////////
    /// @name  Discard Rules Methods
    /// Methods for setting up the automatic discard rules when reading the file
    //@{

    /// When reading records, skip all variants with the ids specified
    /// in the passed in filename.
    /// Returns false, if the file could not be read.
    bool setExcludeIDs(const char* filename);

    /// When reading records, keep only variants with the ids specified
    /// in the passed in filename.
    /// Returns false, if the file could not be read.
    bool setIncludeIDs(const char* filename);

    /// Set which rules should be applied for discarding records.
    /// OR in all discard rules to be applied.
    /// For example:: reader.setDiscards(DISCARD_NON_PHASED | 
    ///                                  DISCARD_MISSING_GT);
    /// NOTE: Discard rules are NOT reset when a file is reset, closed, or a new
    /// one opened, but are reset when this is called.
    void setDiscardRules(uint64_t discards) { myDiscardRules = discards; }
    
    /// Add additional discard rule,  OR in all additional discard rules to
    /// be applied.
    /// For example:: reader.setDiscards(DISCARD_NON_PHASED | 
    ///                                  DISCARD_MISSING_GT);
    /// NOTE: Discard rules are NOT reset when a file is reset, closed, or a new
    /// one opened, and are NOT reset when this is called.
    void addDiscardRules(uint64_t discards) { myDiscardRules |= discards; }
    
    /// Add a discard rule based on the minimum allele count of alternate
    /// alleles in the specified samples.  If the sum of all alternate allele
    /// counts in the specified samples (or in all samples if NULL is passed)
    ///is greater than or equal to the amount specified, the record is kept.
    /// If not, then the record is discarded.
    /// \param minAltAlleleCount minimum count of all alternate alleles for
    /// a record that should be kept, any with fewer will be discarded.
    /// \param subset only count alternate alleles in samples within the
    /// specified subset or NULL if all samples should be counted.  The
    /// pointer is stored in this object, but is not cleaned up by this object.
    void addDiscardMinAltAlleleCount(int32_t minAltAlleleCount,
                                     VcfSubsetSamples* subset);

    /// Remove the discard rule for minimum alternate allele count.
    void rmDiscardMinAltAlleleCount();
    
    /// Add a discard rule based on the minimum allele count of alternate
    /// alleles in the specified samples.  If the sum of all alternate allele
    /// counts in the specified samples (or in all samples if NULL is passed)
    ///is greater than or equal to the amount specified, the record is kept.
    /// If not, then the record is discarded.
    /// \param minAltAlleleCount minimum count of all alternate alleles for
    /// a record that should be kept, any with fewer will be discarded.
    /// \param subset only count alternate alleles in samples within the
    /// specified subset or NULL if all samples should be counted.  The
    /// pointer is stored in this object, but is not cleaned up by this object.
    void addDiscardMinMinorAlleleCount(int32_t minMinorAlleleCount,
                                       VcfSubsetSamples* subset);

    /// Remove the discard rule for minimum alternate allele count.
    void rmDiscardMinMinorAlleleCount();

    //@}

protected: 
    virtual void resetFile();

private:
    VcfFileReader(const VcfFileReader& vcfFileReader);
    VcfFileReader& operator=(const VcfFileReader& vcfFileReader);

    static const int32_t UNSET_MIN_MINOR_ALLELE_COUNT = -1;
    static const int32_t UNSET_MIN_ALT_ALLELE_COUNT = -1;

    // Set1BasedReadSection was called so process the section prior to reading.
    bool processNewSection();

    // New section information.
    Tabix* myVcfIndex;
    bool myNewSection;
    std::string mySectionChrom;
    int32_t mySection1BasedStartPos;
    int32_t mySection1BasedEndPos;
    bool mySectionOverlap;

    VcfRecordDiscardRules myRecordDiscardRules;

    VcfSubsetSamples mySampleSubset;
    bool myUseSubset;

    int32_t myMinAltAlleleCount;
    VcfSubsetSamples* myAltAlleleCountSubset;

    int32_t myMinMinorAlleleCount;
    VcfSubsetSamples* myMinorAlleleCountSubset;

    uint64_t myDiscardRules;

    // Number of records read/written so far that were not discarded.
    int myNumKeptRecords;
    int myTotalRead;
};

#endif
