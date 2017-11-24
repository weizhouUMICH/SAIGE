/*
 *  Copyright (C) 2010-2011  Regents of the University of Michigan
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

#if !defined(_CIGAR_ROLLER_H)
#define _CIGAR_ROLLER_H

#include "Cigar.h"

/// The purpose of this class is to provide accessors for setting, updating, modifying the CIGAR object. It is a child class of Cigar.

///
/// Docs from Sam1.pdf:
///
/// Clipped alignment. In Smith-Waterman alignment, a sequence may not be aligned from the first residue to the last one.
/// Subsequences at the ends may be clipped off. We introduce operation ʻSʼ to describe (softly) clipped alignment. Here is
/// an example. Suppose the clipped alignment is:
/// REF:  AGCTAGCATCGTGTCGCCCGTCTAGCATACGCATGATCGACTGTCAGCTAGTCAGACTAGTCGATCGATGTG
/// READ:        gggGTGTAACC-GACTAGgggg
/// where on the read sequence, bases in uppercase are matches and bases in lowercase are clipped off. The CIGAR for
/// this alignment is: 3S8M1D6M4S.
///
///
/// If the mapping position of the query is not available, RNAME and
/// CIGAR are set as “*”
///
/// A CIGAR string is comprised of a series of operation lengths plus the operations. The conventional CIGAR format allows
/// for three types of operations: M for match or mismatch, I for insertion and D for deletion. The extended CIGAR format
/// further allows four more operations, as is shown in the following table, to describe clipping, padding and splicing:
///
/// op   Description
/// --   -----------
/// M    Match or mismatch
/// I    Insertion to the reference
/// D    Deletion from the reference
/// N    Skipped region from the reference
/// S    Soft clip on the read (clipped sequence present in <seq>)
/// H    Hard clip on the read (clipped sequence NOT present in <seq>)
/// P    Padding (silent deletion from the padded reference sequence)
///



////////////////////////////////////////////////////////////////////////
///
/// CigarRoller is an aid to correctly generating the CIGAR strings
/// necessary to represent how a read maps to the reference.
///
/// It is called once a particular match candidate is being written
/// out, so it is far less performance sensitive than the Smith Waterman
/// code below.
///
class CigarRoller : public Cigar
{
public:

    ////////////////////////////////////////////////////////////////////////
    //
    // Cigar Roller Class
    //
    /// Writes all of the cigar operations contained in this roller to the
    /// passed in stream.
    friend std::ostream &operator << (std::ostream &stream, const CigarRoller& roller);

    /// Default constructor initializes as a CIGAR with no operations.
    CigarRoller()
    {
        clearQueryAndReferenceIndexes();
    }

    /// Constructor that initializes the object with the specified cigarString.
    CigarRoller(const char *cigarString)
    {
        Set(cigarString);
    }

    /// Add the contents of the specified CigarRoller to this object.
    CigarRoller & operator += (CigarRoller &rhs);

    /// Append the specified operator to this object.
    CigarRoller & operator += (const CigarOperator &rhs);

    /// Set this object to be equal to the specified CigarRoller.
    CigarRoller & operator = (CigarRoller &rhs);

    /// Append the specified operation with the specified count to this object.
    void Add(Operation operation, int count);

    /// Append the specified operation with the specified count to this object.
    void Add(char operation, int count);

    /// Append the specified cigarString to this object.
    void Add(const char *cigarString);

    /// Append the specified Cigar object to this object.
    void Add(CigarRoller &rhs)
    {
        (*this) += rhs;
    }

    /// Remove the operation at the specified index.
    /// \return true if successfully removed, false if not.
    bool Remove(int index);

    /// Increments the count for the operation at the specified index
    /// by the specified value, specify a negative value to decrement.
    /// \return true if it is successfully incremented, false if not.
    bool IncrementCount(int index, int increment);

    /// Updates the operation at the specified index to be the specified
    /// operation and have the specified count.
    /// \return true if it is successfully updated, false if not.
    bool Update(int index, Operation op, int count);

    /// Sets this object to the specified cigarString.
    void Set(const char *cigarString);

    /// Sets this object to the BAM formatted cigar found at the beginning
    /// of the specified buffer which is bufferLen long.
    void Set(const uint32_t* cigarBuffer, uint16_t bufferLen);

    //
    // when we examine CIGAR strings, we need to know how
    // many cumulative insert and delete positions there are
    // so that we can adjust the read location appropriately.
    //
    // Here, we iterate over the vector of CIGAR operations,
    // summaring the count for each insert or delete (insert
    // increases the offset, delete decreases it).
    //
    // The use case for this is when we have a genome match
    // position based on an index word other than the first one,
    // and there is also a insert or delete between the beginning
    // of the read and the index word.  We can't simply report
    // the match position without taking into account the indels,
    // otherwise we'll be off by N where N is the sum of this
    // indel count.
    //
    /// DEPRECATED - do not use, there are better ways to accomplish that by 
    /// using read lengths, reference lengths, span of the read, etc.
    int getMatchPositionOffset();

    /// Get the string reprentation of the Cigar operations in this object,
    /// caller must delete the returned value.
    const char *getString();

    /// Clear this object so that it has no Cigar Operations.
    void clear();

private:
};


inline std::ostream &operator << (std::ostream &stream, const CigarRoller& roller)
{
    stream << roller.cigarOperations;
    return stream;
}

#endif
