/*
 *  Copyright (C) 2011-2012  Regents of the University of Michigan,
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

#ifndef __VCF_RECORD_H__
#define __VCF_RECORD_H__

#include <vector>
#include <stdlib.h>
#include "VcfRecordFilter.h"
#include "VcfRecordInfo.h"
#include "VcfRecordGenotype.h"
#include "StatGenStatus.h"
#include "VcfRecordDiscardRules.h"

/// This header file provides interface to read/write VCF files.
class VcfRecord
{
public:
    /// Default Constructor, initializes the variables, but does not open
    /// any files.
    VcfRecord();
    /// Destructor
    virtual ~VcfRecord();
    
    /// Read the next Vcf data line from the file.
    /// \param filePtr IFILE to read from.
    /// \param siteOnly only store the first 8 columns
    /// \param sampleSubset pointer to sample subset information, 
    ///        but NULL if sample subsetting is not to be done.
    /// \param discardRules pointer to record discard information,
    ///        but NULL if record discarding is not to be done.
    /// \return true if a line was successfully read from the specified filePtr,
    /// false if not.
    bool read(IFILE filePtr, bool siteOnly,
              VcfRecordDiscardRules& discardRules,
              VcfSubsetSamples* sampleSubset = NULL);

    /// Write this data line to the file (including the newline).
    /// \param filePtr IFILE to write to.
    /// \param siteOnly only write the first 8 columns
    /// \return true if a line was successfully written to the specified filePtr,
    /// false if not.
    bool write(IFILE filePtr, bool siteOnly);

    /// Reset this header, preparing for a new one.
    void reset();

    /// Returns the status associated with the last method that sets the status.
    /// \return StatGenStatus of the last command that sets status.
    const StatGenStatus& getStatus();

//     bool isValid();

    ///////////////////////
    /// @name  Get Vcf Fields
    /// Get methods for record fields (do not set status).
    //@{
    const char* getChromStr() {return(myChrom.c_str());}
    int get1BasedPosition() {return(my1BasedPosNum);}
    const char* getIDStr() {return(myID.c_str());}
    const char* getRefStr() {return(myRef.c_str());}
    int getNumRefBases() {return(myRef.size());}
    const char* getAltStr() {return(myAlt.c_str());}
    /// Return a pointer to the alleles at the specified index with index 0
    /// being the reference string for this position and index 1 starting
    /// the alternate alleles, throwing an exception if the index is out of
    /// range.
    /// \param index allele index (0 for reference, 1 for first alt, 
    ///  2 for second, etc)
    /// \return string of the alleles at the specified index 
    const char* getAlleles(unsigned int index);
    /// Return the int value of the first allele in the string at the 
    /// specified index with index 0 being the reference string for
    /// this position and index 1 starting the alternate alleles, 
    /// throwing an exception if the index is out of range.
    /// \param index allele index (0 for reference, 1 for first alt, 
    ///  2 for second, etc)
    /// \return int allele at the specified index, 1=A, 2=C, 3=G, 4=T 
    int getIntAllele(unsigned int index);
    /// Return the number of alternates listed in the Alts string.
    unsigned int getNumAlts();

    float getQual() {return(myQualNum);}
    const char* getQualStr() {return(myQual.c_str());}

    /// Return a reference to the filter information.
    VcfRecordFilter& getFilter(){return(myFilter);}
    /// Return whether or not all filters were passed.
    int passedAllFilters() { return(myFilter.passedAllFilters()); }

    /// Get a reference to the information field.
    VcfRecordInfo& getInfo() {return myInfo;}

    /// Get a reference to the genotype fields.
    VcfRecordGenotype& getGenotypeInfo() {return myGenotype;}

    inline int getNumSamples() { return(myGenotype.getNumSamples()); }

    inline int getNumGTs(int index) { return(myGenotype.getNumGTs(index)); }

    inline int getGT(int sampleNum, unsigned int gtIndex)
    { return(myGenotype.getGT(sampleNum, gtIndex)); }

    inline void setGT(int sampleNum, unsigned int gtIndex, int newGt)
    { myGenotype.setGT(sampleNum, gtIndex, newGt); }

    /// Return true if all of the samples are phased and none are unphased,
    /// false if any are unphased or not phased.
    inline bool allPhased() { return(myGenotype.allPhased()); }

    /// Return true if all of the samples are unphased and none are phased,
    /// false if any are phased or not unphased.
    inline bool allUnphased() { return(myGenotype.allUnphased()); }

    /// Return true if all samples of all records have all the genotype alleles
    /// specified, false if not or if any GT field is missing.
    bool hasAllGenotypeAlleles() { return(myGenotype.hasAllGenotypeAlleles()); }

    /// Return the number of occurances of the specified allele index in the
    /// genotypes for this record.  Index 0 for the reference.  The alternate
    /// alleles start with index 1.  An exception is thrown if the index is
    /// out of range.  Optionally, the specified subset of samples can be
    /// skipped when determining allele counts.  (If the record is read with
    /// just a subset of samples, those are automatically excluded here
    /// regardless of the passed in sampleSubset.
    /// \param index allele index (0 for reference, 1 for first alt, 
    ///  2 for second, etc)
    /// \param sampleSubset pointer to sample subset information, 
    ///        but NULL if additional sample subsetting is not to be done.
    /// \return int allele count for the specified ref/alt.
    int getAlleleCount(unsigned int index,
                       VcfSubsetSamples* sampleSubset = NULL);

    //@}


    ///////////////////////
    /// @name  Set Vcf Fields
    /// Set methods for record fields (do not set status).
    //@{
    void setChrom(const char* chrom) {myChrom = chrom;}
    void set1BasedPosition(int pos) {my1BasedPosNum = pos;}
    void setID(const char* id) {myID = id;}
    void setRef(const char* ref) {myRef = ref;}
    void setAlt(const char* alt) {myAlt = alt; myAltArray.clear();}
    //    void setQual(float qual) {myQualNum = qual; }
    void setQual(const char* qual)
    {
        myQual = qual; 
        if(myQual != ".")
        {
            myQualNum = atof(qual);
        }
        else
        {
            myQualNum = -1;
        }
    }

protected: 

    /// Read the specified file until a tab, '\n', or EOF is found
    /// and appending the read characters to stringRef (except for the 
    /// stopping character).
    /// \param filePtr open file to be read.
    /// \param stringRef reference to a string that should be appended to
    /// with the characters read from the file until a '\t', '\n', or EOF
    /// is found.
    /// \return true if a '\t' stopped the reading, false if '\n' or
    /// EOF stopped the reading.
    bool readTilTab(IFILE filePtr, std::string& stringRef);

private:
    VcfRecord(const VcfRecord& vcfRecord);
    VcfRecord& operator=(const VcfRecord& vcfRecord);

    static const char ALT_DELIM = ',';

    std::string myChrom;
    int my1BasedPosNum;
    //std::string my1BasedPos;
    vcfIDtype myID;
    std::string myRef;
    std::string myAlt;
    float myQualNum;
    std::string myQual;
    VcfRecordFilter myFilter;
    VcfRecordInfo myInfo;
    VcfRecordGenotype myGenotype;
    ReusableVector<std::string> myAltArray;
    std::vector<int> myAlleleCount;


    // The status of the last failed command.
    StatGenStatus myStatus;

    std::string myDummyString;
    // Set Pointers to each piece.

};

#endif
