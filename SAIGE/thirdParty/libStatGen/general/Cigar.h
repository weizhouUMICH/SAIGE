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

#if !defined(_CIGAR_H)
#define _CIGAR_H

#include <string.h> // for inline use of strcat, etc
#include <limits.h> // for INT_MAX
#include <stdint.h> // for uint32_t and friends


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <utility>
#include <vector>

#include "Generic.h"
#include "StringBasics.h"

/// This class represents the CIGAR without any methods to set the cigar
/// (see CigarRoller for that).

//
// Docs from Sam1.pdf:
//
// Clipped alignment. In Smith-Waterman alignment, a sequence may not be aligned from the first residue to the last one.
// Subsequences at the ends may be clipped off. We introduce operation ʻSʼ to describe (softly) clipped alignment. Here is
// an example. Suppose the clipped alignment is:
// REF:  AGCTAGCATCGTGTCGCCCGTCTAGCATACGCATGATCGACTGTCAGCTAGTCAGACTAGTCGATCGATGTG
// READ:        gggGTGTAACC-GACTAGgggg
// where on the read sequence, bases in uppercase are matches and bases in lowercase are clipped off. The CIGAR for
// this alignment is: 3S8M1D6M4S.
//
//
// If the mapping position of the query is not available, RNAME and
// CIGAR are set as “*”
//
// A CIGAR string is comprised of a series of operation lengths plus the operations. The conventional CIGAR format allows
// for three types of operations: M for match or mismatch, I for insertion and D for deletion. The extended CIGAR format
// further allows four more operations, as is shown in the following table, to describe clipping, padding and splicing:
//
// op   Description
// --   -----------
// M    Match or mismatch
// I    Insertion to the reference
// D    Deletion from the reference
// N    Skipped region from the reference
// S    Soft clip on the read (clipped sequence present in <seq>)
// H    Hard clip on the read (clipped sequence NOT present in <seq>)
// P    Padding (silent deletion from the padded reference sequence)
//




////////////////////////////////////////////////////////////////////////
///
/// This class represents the CIGAR.  It contains methods for converting
/// to strings and extracting information from the cigar on how a read
/// maps to the reference.
///
/// It only contains read only methods.  There are no ways to set
/// values.  To set a value, a child class must be used.
///
class Cigar
{
public:
    /// Enum for the cigar operations.
    enum Operation {
        none=0, ///< no operation has been set.
        match, ///< match/mismatch operation.  Associated with CIGAR Operation "M"
        mismatch, ///< mismatch operation.  Associated with CIGAR Operation "M"
        insert,  ///< insertion to the reference (the query sequence contains bases that have no corresponding base in the reference).  Associated with CIGAR Operation "I"
        del,  ///< deletion from the reference (the reference contains bases that have no corresponding base in the query sequence).  Associated with CIGAR Operation "D"
        skip,  ///< skipped region from the reference (the reference contains bases that have no corresponding base in the query sequence).  Associated with CIGAR Operation "N"
        softClip,  ///< Soft clip on the read (clipped sequence present in the query sequence, but not in reference).  Associated with CIGAR Operation "S"
        hardClip,  ///< Hard clip on the read (clipped sequence not present in the query sequence or reference).  Associated with CIGAR Operation "H"
        pad ///< Padding (not in reference or query).  Associated with CIGAR Operation "P"
    };

    // The maximum value in the operation enum (used for setting up a bitset of
    // operations.
    static const int MAX_OP_VALUE = pad;

    ////////////////////////////////////////////////////////////////////////
    //
    // Nested Struct : CigarOperator
    //
    struct CigarOperator
    {

        CigarOperator()
        {
            operation = none;
            count = 0;
        }

        /// Set the cigar operator with the specified operation and
        /// count length.
        CigarOperator(Operation operation, uint32_t count)
            : operation(operation), count(count) {};

        Operation operation;

        uint32_t count;

        /// Get the character code (M, I, D, N, S, H, or P) associated with
        /// this operation.
        char getChar() const
        {
            switch (operation)
            {
                case none:
                    return '?';  // error
                case match:
                case mismatch:
                    return'M';
                case insert:
                    return 'I';
                case del:
                    return'D';
                case skip:
                    return 'N';
                case softClip:
                    return 'S';
                case hardClip:
                    return 'H';
                case pad:
                    return 'P';
            }
            return '?'; // actually it is an error to get here
        }

        /// Compare only on the operator, true if they are the same, false if not.  Match and mismatch are considered the same for CIGAR strings.
        bool operator == (const CigarOperator &rhs) const
        {
            if (operation==rhs.operation)
                return true;
            if ((operation == mismatch || operation == match) && (rhs.operation == mismatch || rhs.operation == match))
                return true;
            return false;
        }

        /// Compare only on the operator, false if they are the same, true if not.  Match and mismatch are considered the same for CIGAR strings.
        bool operator != (const CigarOperator &rhs) const
        {
            return !((*this) == rhs) ;
        }

    };

    ////////////////////////////////////////////////////////////////////////
    //
    // Cigar  Class statics
    //

    /// Return true if the specified operation is found in the
    /// reference sequence, false if not.
    static bool foundInReference(Operation op)
    {
        switch(op)
        {
            case match:
            case mismatch:
            case del:
            case skip:
                return true;
            default:
                return false;
        }
        return false;
    }

    /// Return true if the specified operation is found in the
    /// reference sequence, false if not.
    static bool foundInReference(char op)
    {
        switch(op)
        {
            case 'M':
            case '=':
            case 'X':
            case 'D':
            case 'N':
                return true;
            default:
                return false;
        }
        return false;
    }

    /// Return true if the specified operation is found in the
    /// reference sequence, false if not.
    static bool foundInReference(const CigarOperator &op)
    {
        return(foundInReference(op.operation));
    }

    /// Return true if the specified operation is found in the
    /// query sequence, false if not.
    static bool foundInQuery(Operation op)
    {
        switch(op)
        {
            case match:
            case mismatch:
            case insert:
            case softClip:
                return true;
            default:
                return false;
        }
        return false;
    }
    
    /// Return true if the specified operation is found in the
    /// query sequence, false if not.
    static bool foundInQuery(char op)
    {
        switch(op)
        {
            case 'M':
            case '=':
            case 'X':
            case 'I':
            case 'S':
                return true;
            default:
                return false;
        }
        return false;
    }

    /// Return true if the specified operation is found in the
    /// query sequence, false if not.
    static bool foundInQuery(const CigarOperator &op)
    {
        return(foundInQuery(op.operation));
    }
    
    /// Return true if the specified operation is a clipping operation,
    /// false if not.
    static bool isClip(Operation op)
    {
        switch(op)
        {
            case softClip:
            case hardClip:
                return true;
            default:
                return false;
        }
        return false;
    }

    /// Return true if the specified operation is a clipping operation,
    /// false if not.
    static bool isClip(char op)
    {
        switch(op)
        {
            case 'S':
            case 'H':
                return true;
            default:
                return false;
        }
        return false;
    }

    /// Return true if the specified operation is a clipping operation,
    /// false if not.
    static bool isClip(const CigarOperator &op)
    {
        return(isClip(op.operation));
    }

    /// Return true if the specified operation is a match/mismatch operation,
    /// false if not.
    static bool isMatchOrMismatch(Operation op)
    {
        switch(op)
        {
            case match:
            case mismatch:
                return true;
            default:
                return false;
        }
        return false;
    }
    
    /// Return true if the specified operation is a match/mismatch operation,
    /// false if not.
    static bool isMatchOrMismatch(const CigarOperator &op)
    {
        return(isMatchOrMismatch(op.operation));
    }


    ////////////////////////////////////////////////////////////////////////
    //
    // Cigar  Class non static
    //
    friend std::ostream &operator << (std::ostream &stream, const Cigar& cigar);

    /// Default constructor initializes as a CIGAR with no operations.
    Cigar()
    {
        clearQueryAndReferenceIndexes();
    }


    /// Set the passed in String to the string reprentation of the Cigar
    /// operations in this object.
    void getCigarString(String& cigarString) const;

    /// Set the passed in std::string to the string reprentation of the Cigar
    /// operations in this object.
    void getCigarString(std::string& cigarString) const;

    /// Sets the specified string to a valid CIGAR string of characters that
    /// represent the cigar with no digits (a CIGAR of "3M" would return "MMM").
    /// The returned string is actually also a valid CIGAR string.
    /// In theory this makes it easier to parse some reads.
    /// \return s the string to populate
    void getExpandedString(std::string &s) const;

    /// Return the Cigar Operation at the specified index (starting at 0).
    const CigarOperator & operator [](int i) const
    {
        return cigarOperations[i];
    }

    /// Return the Cigar Operation at the specified index (starting at 0).
    const CigarOperator & getOperator(int i) const
    {
        return cigarOperations[i];
    }

    /// Return true if the 2 Cigars are the same
    /// (the same operations of the same sizes).
    bool operator == (Cigar &rhs) const;

    /// Return the number of cigar operations
    int size()  const
    {
        return cigarOperations.size();
    }

    /// Write this object as a string to cout.
    void Dump() const
    {
        String cigarString;
        getCigarString(cigarString);
        std::cout << cigarString ;
    }

    /// Return the length of the read that corresponds to
    /// the current CIGAR string.
    ///
    /// For validation, we should expect that a sequence
    /// read in a SAM file will be the same length as the
    /// value returned by this method.
    ///
    /// Example: 3M2D3M describes a read with three bases
    /// matching the reference, then skips 2 bases, then has
    /// three more bases that match the reference (match/mismatch).
    /// In this case, the read length is expected to be 6.
    ///
    /// Example: 3M2I3M describes a read with 3 match/mismatch
    /// bases, two extra bases, and then 3 more match/mistmatch
    /// bases.  The total in this example is 8 bases.
    ///
    /// \return returns the expected read length
    int getExpectedQueryBaseCount() const;

    /// Return the number of bases in the reference that
    /// this CIGAR "spans".
    ///
    /// When doing range checking, we occassionally need to know
    /// how many total bases the CIGAR string represents as compared
    /// to the reference.
    ///
    /// Examples: 3M2D3M describes a read that overlays 8 bases in
    /// the reference.  3M2I3M describes a read with 3 bases that
    /// match the reference, two additional bases that aren't in the
    /// reference, and 3 more bases that match the reference, so it
    /// spans 6 bases in the reference.
    ///
    /// \return how many bases in the reference are spanned
    /// by the given CIGAR string
    ///
    int getExpectedReferenceBaseCount() const;

    /// Return the number of clips that are at the beginning of the cigar.
    int getNumBeginClips() const;

    /// Return the number of clips that are at the end of the cigar.
    int getNumEndClips() const;

    /// Return the reference offset associated with the specified
    /// query index or INDEX_NA based on this cigar.
    int32_t getRefOffset(int32_t queryIndex);

    /// Return the query index associated with the specified
    /// reference offset or INDEX_NA based on this cigar.
    int32_t getQueryIndex(int32_t refOffset);

    /// Return the reference position associated with the specified query index
    /// or INDEX_NA based on this cigar and the specified queryStartPos which
    /// is the leftmost mapping position of the first matching base in the
    /// query.
    int32_t getRefPosition(int32_t queryIndex, int32_t queryStartPos);

    /// Return the query index or INDEX_NA associated with the specified 
    /// reference offset when the query starts at the specified reference
    /// position.
    int32_t getQueryIndex(int32_t refPosition, int32_t queryStartPos);

    /// Returns the index into the expanded cigar for the cigar 
    /// associated with the specified queryIndex.
    /// INDEX_NA returned if the index is out of range.
    int32_t getExpandedCigarIndexFromQueryIndex(int32_t queryIndex);

    /// Returns the index into the expanded cigar for the cigar 
    /// associated with the specified reference offset.
    /// INDEX_NA returned if the offset is out of range.
    int32_t getExpandedCigarIndexFromRefOffset(int32_t refOffset);

    /// Returns the index into the expanded cigar for the cigar 
    /// associated with the specified reference position and queryStartPos.
    /// INDEX_NA returned if the position is out of range.
    int32_t getExpandedCigarIndexFromRefPos(int32_t refPosition, 
                                            int32_t queryStartPos);

    /// Return the character code of the cigar operator associated with the
    /// specified expanded CIGAR index.  '?' is returned for an out of range
    /// index.
    char getCigarCharOp(int32_t expandedCigarIndex);

    /// Return the character code of the cigar operator associated with
    /// the specified queryIndex.  '?' is returned for an out of range index.
    char getCigarCharOpFromQueryIndex(int32_t queryIndex);

    /// Return the character code of the cigar operator associated with
    /// the specified reference offset.  '?' is returned for an out of range offset.
    char getCigarCharOpFromRefOffset(int32_t refOffset);

    /// Return the character code of the cigar operator associated with
    /// the specified reference position.  '?' is returned for an out of
    /// range reference position.
    char getCigarCharOpFromRefPos(int32_t refPosition, int32_t queryStartPos);

    /// Return the number of bases that overlap the reference and the
    /// read associated with this cigar that falls within the specified region.
    /// \param start : inclusive 0-based start position (reference position) of
    ///         the region to check for overlaps in
    ///         (-1 indicates to start at the beginning of the reference.)
    /// \param end  : exclusive 0-based end position (reference position) of the
    ///          region to check for overlaps in
    ///         (-1 indicates to go to the end of the reference.)
    /// \param queryStartPos : 0-based leftmost mapping position of the first 
    ///                 matcihng base in the query.
    uint32_t getNumOverlaps(int32_t start, int32_t end, int32_t queryStartPos);

    /// Return whether or not the cigar has indels (insertions or delections)
    /// \return true if it has an insertion or deletion, false if not.
    bool hasIndel();

    /// Value associated with an index that is not applicable/does not exist,
    /// used for converting between query and reference indexes/offsets when
    /// an associated index/offset does not exist.
    static const int32_t INDEX_NA;

protected:
    // Clear the query index/reference offset index vectors.
    void clearQueryAndReferenceIndexes();

    // Set the query index/reference offset index vectors.
    void setQueryAndReferenceIndexes();

    // Container for the cigar operations in this cigar.
    std::vector<CigarOperator> cigarOperations;

private:
    // The vector is indexed by query index and contains the reference
    // offset associated with that query index.
    // The vector is reset each time a new cigar operation is added, and
    // is calculated when accessed if it is not already set.
    std::vector<int32_t> queryToRef;

    // The vector is indexed by reference offset and contains the query
    // index associated with that reference offset.
    // The vector is reset each time a new cigar operation is added, and
    // is calculated when accessed if it is not already set.
    std::vector<int32_t> refToQuery;

    // The vector is indexed by reference offset and contains the offset into
    // the expanded cigar associated with that reference offset.
    // The vector is reset each time a new cigar operation is added, and
    // is calculated when accessed if it is not already set.
    std::vector<int32_t> refToCigar;

    // The vector is indexed by query index and contains the offset into
    // the expanded cigar associated with that query index.
    // The vector is reset each time a new cigar operation is added, and
    // is calculated when accessed if it is not already set.
    std::vector<int32_t> queryToCigar;

    std::string myExpandedCigar;
};

/// Writes the specified cigar operation to the specified stream as <count><char> (3M).
inline std::ostream &operator << (std::ostream &stream, const Cigar::CigarOperator& o)
{
    stream << o.count << o.getChar();
    return stream;
}

/// Writes all of the cigar operations contained in the cigar to the passed in stream.
inline std::ostream &operator << (std::ostream &stream, const Cigar& cigar)
{
    stream << cigar.cigarOperations;
    return stream;
}

#endif
