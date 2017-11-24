/*
 *  Copyright (C) 2010  Regents of the University of Michigan
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

#ifndef __GLF_RECORD_H__
#define __GLF_RECORD_H__

#include <map>
#include <stdint.h>

#include "InputFile.h" 
#include "CharBuffer.h"

/// This class allows a user to easily get/set the fields in a GLF record. 
class GlfRecord
{
public:
    /// Constructor
    GlfRecord();

    /// Destructor
    ~GlfRecord();

//     // Copy Constructor   
//     GlfRecord(const GlfRecord& record);

//     // Overload operator = to copy the passed in record into this record.
//     GlfRecord & operator = (const GlfRecord& record);

//     // Overload operator = to copy the passed in record into this record.
//     bool copy(const GlfRecord& record);

    /// Clear this record back to the default setting.
    void reset();
   
    /// Read the record from the specified file (file MUST be in
    /// the correct position for reading a record).
    /// \param filePtr file to read from that is in the correct position.
    /// \return true if the record was successfully read from the file (even
    /// if it is an endMarker), false if it was not successfully read.
    bool read(IFILE filePtr);

    /// Write the record to the specified file.
    /// \param filePtr file to write to that is in the correct position.
    /// \return true if the record was successfully written to the 
    /// file, false if not.
    bool write(IFILE filePtr) const;

    /// Print the reference section in a readable format.
    void print() const;

    /// @name Generic Accessors for Record Types 1 & 2
    //@{
    /// Set the record type and reference base.
    /// \param rtypeRef record type & reference base. Formatted as:
    /// record_type<<4|numeric_ref_base.
    /// \return true if the record type and reference base were successfully
    /// set, false if not.
    bool setRtypeRef(uint8_t rtypeRef);

    /// Set the record type.
    /// \param recType record type: 1 - simple likelihood record, 
    /// 2 - indel likelihood record, 0 - end maker
    /// \return true if the record type was successfully set, false if not.
    bool setRecordType(uint8_t recType);

    /// Set the reference base from an integer value.
    /// \param refBase integer representation of the reference base.
    /// \anchor BaseCharacterIntMap
    /// <table>
    /// <tr><th>Int Value</th><td>0</td><td>1</td><td>2</td><td>3</td><td>4</td><td>5</td><td>6</td><td>7</td><td>8</td><td>9</td><td>10</td><td>11</td><td>12</td><td>13</td><td>14</td><td>15</td></tr>
    /// <tr><th>Character Base</th><td>X</td><td>A</td><td>C</td><td>M</td><td>G</td><td>R</td><td>S</td><td>V</td><td>T</td><td>W</td><td>Y</td><td>H</td><td>K</td><td>D</td><td>B</td><td>N</td></tr>
    /// </table>
    /// \return true if the reference base was successfully set, false if not.
    bool setRefBaseInt(uint8_t refBase);

    // TODO   bool setRefBaseChar(char refBase);

    /// Set the offset from the precedent record.
    /// 0-based coordinate of the record minus the coordinate of the
    /// precedent record. For the first record in a reference sequence,
    /// the previous coordinate is 0.
    /// For insertions between x & x+1, the coordinate is x.
    /// For deletions between x & y, the coordinate is x. 
    /// \param offset offset from the precedent record.
    /// \return true if successfully set, false if not.
    bool setOffset(uint32_t offset);

    /// Set the minimum likelihood and the read depth.
    /// \param minDepth minimum likelihood and read depth. Formatted as:
    /// min_lk<<24|read_dpeth. (min_lk capped at 255)
    /// \return true if successfully set, false if not.
    bool setMinDepth(uint32_t minDepth);

    /// Set the minimum likelihood.
    /// \param minLk minimum likelihood (capped at 255).
    /// \return true if successfully set, false if not.
    bool setMinLk(uint8_t minLk);

    /// Set the the read depth.
    /// \param readDepth read depth.
    /// \return true if successfully set, false if not.
    bool setReadDepth(uint32_t readDepth);

    /// Set the RMS of mapping qualities of reads covering the site.
    /// \param rmsMapQ RMS of mapping qualities
    /// \return true if successfully set, false if not.
    bool setRmsMapQ(uint8_t rmsMapQ);
 
    /// Return the record type.
    /// \return record type for this record: 0 - endMarker, 
    /// 1 - simple likelihood, 2 - indel likelihood
    inline int getRecordType() const
    {
        return(myRecTypeRefBase >> REC_TYPE_SHIFT);
    }

    /// Return the reference base as an integer.
    /// \return integer representation of the reference base.
    /// See: \ref BaseCharacterIntMap
    inline int getRefBase() const
    {
        return(myRecTypeRefBase & REF_BASE_MASK);
    }

    /// Return the reference base as a character.
    /// \return character representation of the reference base.
    char getRefBaseChar() const;

    /// Return the offset from the precedent record.
    /// \return offset from the precedent record.
    uint32_t getOffset() const;

    /// Return the minimum likelihood and read depth.  Formatted as:
    /// min_lk<<24|read_dpeth. (min_lk capped at 255)
    /// \return minimum likelihood and read depth
    uint32_t getMinDepth() const;

    /// Return the minimum likelihood
    /// \return minimum likelihood
    uint8_t getMinLk() const;

    /// Return the read depth.
    /// \return read depth
    uint32_t getReadDepth() const;

    /// Return the RMS of mapping qualities of reads covering the site.
    /// \return RMS of maping qualities.
    uint8_t getRmsMapQ() const;

    //@}
    
    /// @name Record Type 1 Accessors
    /// Record Type 1: Simple Likelihood Record
    //@{
    //bool setType1(all fields for type 1);

    /// Set the likelihood for the specified genotype.
    /// Throws an exception if index is out of range.
    /// \param index index for the genotype for which the likelihood is 
    /// being set.
    /// \anchor GenotypeIndexTable
    /// <table>
    /// <tr><th>Index</th><td>0</td><td>1</td><td>2</td><td>3</td><td>4</td><td>5</td><td>6</td><td>7</td><td>8</td><td>9</td></tr>
    /// <tr><th>Genotype</th><td>AA</td><td>AC</td><td>AG</td><td>AT</td><td>CC</td><td>CG</td><td>CT</td><td>GG</td><td>GT</td><td>TT</td></tr>
    /// </table>
    /// \param value likelihood for the genotype at the specified index.
    /// \return true if successfully set, false if not.
    bool setLk(int index, uint8_t value);

    //bool getType1(all fields for type 1);

    /// Get the likelihood for the specified genotype index.
    /// Throws an exception if index is out of range.
    /// \param index index of the genotype for which the likelihood should
    /// be returned.  See: \ref GenotypeIndexTable
    /// \return likelihood of the specified index.
    uint8_t getLk(int index);    
    //@}

    /// @name Record Type 2 Accessors
    /// Record Type2: Indel Likelihood Record
    //@{
//     bool setType2(all fields for type 2);

    /// Set the likelihood of the first homozygous indel allele.
    /// \param lk likelihood of the 1st homozygous indel allele (capped at 255)
    /// \return true if successfully set, false if not.
    bool setLkHom1(uint8_t lk);

    /// Set the likelihood of the 2nd homozygous indel allele.
    /// \param lk likelihood of the 2nd homozygous indel allele (capped at 255)
    /// \return true if successfully set, false if not.
    bool setLkHom2(uint8_t lk);

    /// Set the likelihood of a heterozygote.
    /// \param lk likelihood of a heterozygote (capped at 255)
    /// \return true if successfully set, false if not.
    bool setLkHet(uint8_t lk);

    /// Set the sequence of the first indel allele if the
    /// first indel is an insertion.
    /// \param indelSeq sequence of the first indel allele (insertion).
    /// \return true if successfully set, false if not.
    bool setInsertionIndel1(const std::string& indelSeq);

    /// Set the sequence of the first indel allele if the
    /// first indel is an deletion.
    /// \param indelSeq sequence of the first indel allele (deletion).
    /// \return true if successfully set, false if not.
    bool setDeletionIndel1(const std::string& indelSeq);

    /// Set the sequence of the 2nd indel allele if the
    /// 2nd indel is an insertion.
    /// \param indelSeq sequence of the 2nd indel allele (insertion).
    /// \return true if successfully set, false if not.
    bool setInsertionIndel2(const std::string& indelSeq);

    /// Set the sequence of the 2nd indel allele if the
    /// 2nd indel is an deletion.
    /// \param indelSeq sequence of the 2nd indel allele (deletion).
    /// \return true if successfully set, false if not.
    bool setDeletionIndel2(const std::string& indelSeq);

    //     bool setType2(all fields for type 2);

    /// Return the likelihood of the 1st homozygous indel allele.
    /// \return likelihood of the 1st homozygous indel allele.
    uint8_t getLkHom1();

    /// Return the likelihood of the 2nd homozygous indel allele.
    /// \return likelihood of the 2nd homozygous indel allele.
    uint8_t getLkHom2();

    /// Return the likelihood of a heterozygote.
    /// \return likelihood of a hetereozygote.
    uint8_t getLkHet();

    /// Get the sequence and length (+:ins, -:del) of the 1st indel allele.
    /// \param indelSeq string to set with the sequence of the 1st indel allele
    /// \return length of the 1st indel allele
    /// (positive=insertion; negative=deletion; 0=no-indel)
    int16_t getIndel1(std::string& indelSeq);

    /// Get the sequence and length (+:ins, -:del) of the 2nd indel allele.
    /// \param indelSeq string to set with the sequence of the 2nd indel allele
    /// \return length of the 2nd indel allele
    /// (positive=insertion; negative=deletion; 0=no-indel)
    int16_t getIndel2(std::string& indelSeq);
    //@}

private:
    // Read a record of record type 1.
    void readType1(IFILE filePtr);

    // Read a record of record type 2.
    void readType2(IFILE filePtr);


    // Write the rtyperef field.
    void writeRtypeRef(IFILE filePtr) const;


    // Write a record of record type 1.
    void writeType1(IFILE filePtr) const;

    // Write a record of record type 2.
    void writeType2(IFILE filePtr) const;

    // Contains record_type and ref_base.
    uint8_t myRecTypeRefBase;

    static const uint8_t REC_TYPE_SHIFT = 4;
    static const uint8_t REF_BASE_MASK = 0xF;
    static const uint8_t REC_TYPE_MASK = 0xF0;

    static const uint32_t MIN_LK_SHIFT = 24;
    static const uint32_t READ_DEPTH_MASK = 0xFFFFFF;
    static const uint32_t MIN_LK_MASK = 0xFF000000;

    static const char REF_BASE_MAX = 15;
    static std::string REF_BASE_CHAR;

    static const int NUM_REC1_LIKELIHOOD = 10;

    struct
    {
        uint32_t offset;
        uint32_t min_depth;
        uint8_t rmsMapQ;
        uint8_t lk[GlfRecord::NUM_REC1_LIKELIHOOD];
    } myRec1Base;

    static const int REC1_BASE_SIZE = 19;

    struct
    {
        uint32_t offset;
        uint32_t min_depth;
        uint8_t rmsMapQ;
        uint8_t lkHom1;
        uint8_t lkHom2;
        uint8_t lkHet;
        int16_t indelLen1;
        int16_t indelLen2;
    } myRec2Base;

    // TODO rest of rec 2.
    CharBuffer myIndelSeq1;
    CharBuffer myIndelSeq2;

    static const int REC2_BASE_SIZE = 16;

};

#endif
