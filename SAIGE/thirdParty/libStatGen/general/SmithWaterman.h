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

#if !defined(_SMITH_WATERMAN_)
#define _SMITH_WATERMAN_

#include <string.h> // for inline use of strcat, etc
#include <limits.h> // for INT_MAX
#include <stdint.h> // for uint32_t and friends

//
// This file implements a bi-directional, banded Smith Waterman matching
// algorithm.
//
// The design is dictated by several observations:
//
//  - building the full matrix H for the read and reference is slow,
//    so we perform it only for a band down the diagonal, thus speeding
//    up the algorithm at the cost of reduced detection of insert/deletes.
//
//  - we must minimize data copying
//
//  - we must have the ability to test the algorithm quickly and simply,
//    hence we implement the class as a template so that this file doesn't
//    have to depend on Karma's GenomeSequence object, which is a relatively
//    heavy weight object to do testing against.
//
//  - because Karma uses index words to determine match candidate locations
//    across the genome, and because we use the banded Smith Waterman approach,
//    we must provide bi-directional Smith Waterman matching.  See example below.
//
// To fully understand the examples below, make sure you understand
// Phred quality scores, CIGAR strings, and pattern matching.
//
// Simple Functional Examples:
//
//   Given a read, a read quality and a reference, we want to obtain some
//   measure of how well that read matches the reference at that given location.
//   So, we have:
//          Read: "ATCG"
//       Quality: "$$$$" (Phred scores - all are decimal 3)
//     Reference: "ATCG"
//     We expect: sumQ = 0, and Cigar "4M"
//
// Complex Functional Examples:
//          Read: "AATCG"
//       Quality: "$$$$$" (Phred scores - all are decimal 3)
//     Reference: "ATCG"
//     We expect: sumQ = 3, and Cigar "1M1I3M"
//
// Backwards matching:
// It is harder for me to construct a clear example, so imagine a read with
// cumulative inserts or deletes that sum up to a number larger than the
// width of the band along the diagonal.  If we perform SW in the 'wrong'
// direction, we will lose the read entirely, whereas if we start from the
// end that is known to be matching, we may obtain a good match for the bulk
// of the read.
//
// For example, you could imagine a read where the first 10 bases had a mess
// of inserts, but then was clean for the next 100 bases.  You'd want it
// to match as many of the good bases as practical, even if you knew you were
// losing information at the end.
//

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <utility>
#include <vector>

#include "CigarRoller.h"
#include "Generic.h"

using std::cout;
using std::cin;
using std::cerr;
using std::setw;
using std::endl;
using std::pair;
using std::vector;

#if !defined(MAX)
#define MAX(x,y) ((x)>(y) ? (x) : (y))
#endif
#if !defined(MIN)
#define MIN(x,y) ((x)<(y) ? (x) : (y))
#endif

//
// Implement the core of Smith Waterman as described in:
// http://en.wikipedia.org/wiki/Smith_waterman
//
// The only variation from the basic SW algorithm is the
// use of a banded approach - to limit the algorithm to
// a band along the diagonal.  This limits the maximum
// additive number of indels, but allows an O(c*M) approach,
// where c is the constant max number of indels allowed.
//
// This is implemented as function templates because for testing, it is easier
// to use character arrays.  In our mapper, we will be using a
// combination of a String object for the read and the genome object
// as the reference.  Both can be indexed, and give a character (base)
// value, but the code would be duplicated if we implement SW for
// each type of argument.
//
// Htype -> the type of the array H cell (8 or 16 bit unsigned int)
// Atype -> the read string type (must have Atype::operator [] defined)
//
template <int maxReadLengthH, int maxReferenceLengthH, typename HCellType, typename Atype, typename Btype, typename QualityType, typename readIndexType, typename referenceIndexType>
class SmithWaterman
{
public:

    //
    // XXX in theory, this weight should be sensitive
    // to the quality of the base, and should have
    // an appropriate cost for an indel, as well.
    //
    // I think we need to get rid of this, since it is
    // basically wrong for our needs.
    //
    struct weight
    {
        weight()
        {
            match=2;
            misMatch=-1;
            insert=-1;
            del=-1;
        };

        int match;
        int misMatch;
        int insert;
        int del;
    };

    HCellType   H[maxReadLengthH][maxReferenceLengthH];
    Atype   *A;
    Btype   *B;
    QualityType *qualities;

    int     m,n;
    readIndexType MOffset; // constant offset to m (read)
    referenceIndexType NOffset; // constant offset to n (reference)
    weight  w;
    int     allowedInsertDelete;
    int     direction;
    int     gapOpenCount;
    int     gapCloseCount;
    int     gapExtendCount;
    vector<pair<int,int> > alignment;
    void clearAlignment()
    {
        alignment.clear();
    }

    HCellType maxCostValue;    // max Cost value in H
    pair<int,int>   maxCostPosition;    // row/col of max cost value in H

    //
    // Clear the member variables only.
    // To clear H, call clearH().
    //
    // In theory, clear() plus set() should be sufficient to
    // get a clean run, but I haven't tested this extensively.
    //
    void clear()
    {
        maxCostPosition.first = 0;
        maxCostPosition.second = 0;
        A = NULL;
        B = NULL;
        qualities = NULL;
        m = 0;
        n = 0;
        MOffset = 0;
        NOffset = 0;
        allowedInsertDelete = 0;
        direction = 0;
        gapOpenCount = 0;
        gapCloseCount = 0;
        gapExtendCount = 0;
    }

    // caller will be using set* methods to set everything up.
    SmithWaterman()
    {
        clear();
        clearH();
    }

    // construct with everything and the kitchen sink:
    SmithWaterman(
        Atype *A,
        QualityType *qualities,
        Btype *B,
        int m,
        int n,
        int allowedInsertDelete = INT_MAX,
        int direction = 1,
        readIndexType MOffset = 0,
        referenceIndexType NOffset = 0):
            A(A),
            qualities(qualities),
            B(B),
            m(m),
            n(n),
            allowedInsertDelete(allowedInsertDelete),
            direction(direction),
            MOffset(MOffset),
            NOffset(NOffset),
            maxCostValue((HCellType) 0)
    {
    }

    void setRead(Atype *A)
    {
        this->A = A;
    }
    void setReadQuality(QualityType *qualities)
    {
        this->qualities = qualities;
    }
    void setReference(Btype *B)
    {
        this->B = B;
    }

    // Caller may wish to index into the read to do the matching against
    // only part of the read.
    // NB: the quality length and offset are the same as the read.
    void setReadLength(int m)
    {
        this->m = m;
    }
    void setReadOffset(readIndexType MOffset)
    {
        this->MOffset = MOffset;
    }

    // The reference is typically long, and not necessarily a char *,
    // so we provide an offset here.  If it were always a char *,
    // we'd just modify the caller to point directly at the reference
    // location.
    void setReferenceLength(int n)
    {
        this->n = n;
    }
    void setReferenceOffset(referenceIndexType NOffset)
    {
        this->NOffset = NOffset;
    }

    //
    // Configuration: how wide is the band on the diagonal?
    // We should keep this small -- 1, 2, 3 or similar.  If
    // the value is default (INT_MAX), then the full matrix
    // will be built, which is fine, but quite slow.
    //
    // If this paramater is made smaller than when a previous
    // call to populateH was made, clearH will also need to be called.
    //
    void setAllowedInsertDelete(int allowedInsertDelete = INT_MAX)
    {
        this->allowedInsertDelete = allowedInsertDelete;
    }

    //
    // Configuration: which end do we begin performing SW matching
    // from?  We need this because of index 'anchors' in the karma
    // matcher.
    void setDirection(int direction)
    {
        this->direction = direction;
    }

    void clearH()
    {
        memset(H, 0, sizeof(H));
    }

    void populateH()
    {

        maxCostValue = 0;

        for (int i=1; i<=m ; i++)
        {

            // implement a banded Smith-Waterman approach:
            int low = MAX(1, i - allowedInsertDelete);
            int high = MIN(n, i + allowedInsertDelete);

            for (int j=low; j<=high ; j++)
            {
                HCellType c;
                c = 0;
                if (direction>0) c = MAX(c, H[i-1][j-1] + (((*A)[MOffset + i-1]==(*B)[NOffset + j-1]) ? w.match : w.misMatch));
                else c = MAX(c, H[i-1][j-1] + (((*A)[MOffset + m-i+0]==(*B)[NOffset + n-j+0]) ? w.match : w.misMatch));
                c = MAX(c, H[i-1][j] + w.del);
                c = MAX(c, H[i][j-1] + w.insert);
                H[i][j] = c;
                if (c>maxCostValue)
                {
                    maxCostValue = c;
                    maxCostPosition.first = i;
                    maxCostPosition.second = j;
                }
            }
        }
    }

//
// Given the matrix H as filled in by above routine, print it out.
//
    void printH(bool prettyPrint = true)
    {
        // print the scoring matrix:
        for (int i=-1; i<=m ; i++)
        {
            for (int j=-1; j<=n ; j++)
            {
                if (prettyPrint) cout << setw(3);
                if (i==-1 && j==-1)
                {
                    if (prettyPrint) cout << " ";
                    else cout << "\t";
                }
                else if (j==-1)
                {
                    if (!prettyPrint) cout << "\t";
                    if (i==0) cout << "-";
                    else cout << (*A)[MOffset + direction>0 ? i-1 : m - i];
                }
                else if (i==-1)
                {
                    if (!prettyPrint) cout << "\t";
                    if (j==0) cout << "-";
                    else cout << (*B)[NOffset + direction>0 ? j-1 : n - j];
                }
                else
                {
                    if (!prettyPrint) cout << "\t";
                    cout << H[i][j];
                }
            }
            cout << endl;
        }
    }

    void debugPrint(bool doPrintH = true)
    {
        if (doPrintH) printH();
        cout << "maxCostPosition = " << maxCostPosition << std::endl;
        if (alignment.empty()) cout << "alignment vector is empty.\n";
        else
        {
            cout << "alignment vector:\n";
            for (vector<pair<int,int> >::iterator i=alignment.begin(); i < alignment.end(); i++)
            {
                cout << (i - alignment.begin()) << ": " << *i << "\n";
            }
        }
        cout << std::endl;
    }

//
// Given the Matrix H as filled in by populateH, fill in the
// alignment vector with the indeces of the optimal match.
//
    void populateAlignment()
    {
        alignment.clear();
        int i = m, j = n;

        i = maxCostPosition.first;
        j = maxCostPosition.second;

        //
        // Stop when we either reach zero cost cell or
        // when we reach the upper left corner of H.
        // A zero cost cell to the lower right means we
        // are soft clipping that end.
        //
        while (H[i][j] > 0 || (i>0 && j>0))
        {
// #define DEBUG_ALIGNMENT_VECTOR
#if defined(DEBUG_ALIGNMENT_VECTOR)
            cout << "alignment.push_back(" << i << ", " << j << ")" << endl;
#endif
            alignment.push_back(pair<int,int>(i,j));
            if (H[i-1][j-1]>=H[i-1][j] && H[i-1][j-1]>=H[i][j-1])
            {
                // diagonal upper left cell is biggest
                i--;
                j--;
            }
            else if (H[i-1][j] < H[i][j-1])
            {
                // upper cell is biggest
                j--;
            }
            else
            {
                // left cell is biggest
                i--;
            }
        }
        alignment.push_back(pair<int,int>(i,j));
#if defined(DEBUG_ALIGNMENT_VECTOR)
        cout << "alignment.push_back(" << i << ", " << j << ")" << endl;
        cout << "alignment.size(): " << alignment.size() << endl;
#endif
    }

    //
    // Compute the sumQ for a read that has been mapped using populateH().
    //
    // In the simplest case, the read lies on the diagonal of the
    // matrix H, which means it has only matches and mismatches:
    // no inserts or deletes.
    //
    // However, in general, it is possible to have 0 or more insert,
    // delete, mismatch and soft clipped bases in the read, so we
    // need to accomodate all of those variations.
    //
    // XXX finish this.
    //

    int getSumQ()
    {
        if (direction>0) return getSumQForward();
        else return getSumQBackward();
    }

    int getSumQForward()
    {
        int sumQ = 0;
        vector<pair<int,int> >::reverse_iterator i;

        for (i=alignment.rbegin(); i < alignment.rend() - 1; i++)
        {
// #define DEBUG_GETSUMQ
#if defined(DEBUG_GETSUMQ)
            cout << *i << ": ";
#endif
            if ((*(i+1)).first == ((*i).first+1) && (*(i+1)).second == ((*i).second + 1))
            {
                // match/mismatch
#if defined(DEBUG_GETSUMQ)
                cout << "Match/Mismatch";
#endif
                if ((*A)[MOffset + (*i).first] != (*B)[NOffset + (*i).second])
                    sumQ += (*qualities)[MOffset + (*i).first] - '!';
            }
            else if ((*(i+1)).first == ((*i).first+1) && (*(i+1)).second == ((*i).second))
            {
                // insert?
#if defined(DEBUG_GETSUMQ)
                cout << "Insert";
#endif
                sumQ += 50;
            }
            else if ((*(i+1)).first == ((*i).first) && (*(i+1)).second == ((*i).second + 1))
            {
                // delete?
#if defined(DEBUG_GETSUMQ)
                cout << "Delete";
#endif
                sumQ += 50;
            }
        }
#if defined(DEBUG_GETSUMQ)
        cout << endl;
#endif
        return sumQ;
    }

    int getSumQBackward()
    {
        int sumQ = 0;
        vector<pair<int,int> >::iterator i;

        for (i=alignment.begin(); i < alignment.end() - 1; i++)
        {
#if defined(DEBUG_GETSUMQ)
            cout << *i << ": ";
#endif
            if ((*(i+1)).first == ((*i).first-1) && (*(i+1)).second == ((*i).second - 1))
            {
                // match/mismatch
#if defined(DEBUG_GETSUMQ)
                cout << "Match/Mismatch";
#endif
                if ((*A)[MOffset + m - (*i).first] != (*B)[NOffset + n - (*i).second])
                    sumQ += (*qualities)[MOffset + m - (*i).first] - '!';
            }
            else if ((*(i+1)).first == ((*i).first-1) && (*(i+1)).second == ((*i).second))
            {
                // insert?
#if defined(DEBUG_GETSUMQ)
                cout << "Insert?";
#endif
                sumQ += 50;
            }
            else if ((*(i+1)).first == ((*i).first) && (*(i+1)).second == ((*i).second - 1))
            {
                // delete?
#if defined(DEBUG_GETSUMQ)
                cout << "Delete?";
#endif
                sumQ += 50;
            }
        }
#if defined(DEBUG_GETSUMQ)
        cout << endl;
#endif
        return sumQ;
    }

#if 0
    int getSumQ()
    {
        vector<pair<int,int> >::reverse_iterator i;
        int sumQ = 0;
        for (i=alignment.rbegin(); i < alignment.rend() - 1; i++)
        {
#if defined(DEBUG_ALIGNMENT_VECTOR)
            cout << "i: " << i - alignment.rbegin() << *i << endl;
#endif
            // XXX NOT THIS SIMPLE - need to account for indels
            if (direction>0)
            {
                if ((*A)[MOffset + (*i).first] != (*B)[NOffset + (*i).second])
                    sumQ += (*qualities)[MOffset + (*i).first] - '!';
            }
            else
            {
                // m and n are sizes, first and second are 1 based offsets
                if ((*A)[MOffset + m - (*i).first] != (*B)[NOffset + n - (*i).second])
                    sumQ += (*qualities)[MOffset + m - (*i).first] - '!';
            }
        }
        return sumQ;
    }
#endif

    //
    // Append cigar operations to an existing cigar list.
    //
    // XXX we no longer need the CigarRoller += methods.
    //
    // In this case, the Smith Waterman array H was created from
    // the read and reference in the forward direction.
    //
    void rollCigarForward(CigarRoller &cigar)
    {
        vector<pair<int,int> >::reverse_iterator i;

        for (i=alignment.rbegin(); i < alignment.rend() - 1; i++)
        {
// #define DEBUG_CIGAR
#if defined(DEBUG_CIGAR)
            cout << *i << ": ";
#endif
            if ((*(i+1)).first == ((*i).first+1) && (*(i+1)).second == ((*i).second + 1))
            {
                // match/mismatch
#if defined(DEBUG_CIGAR)
                cout << "Match/Mismatch";
#endif
                cigar.Add(CigarRoller::match, 1);
            }
            else if ((*(i+1)).first == ((*i).first+1) && (*(i+1)).second == ((*i).second))
            {
                // insert?
#if defined(DEBUG_CIGAR)
                cout << "Insert";
#endif
                cigar.Add(CigarRoller::insert, 1);
            }
            else if ((*(i+1)).first == ((*i).first) && (*(i+1)).second == ((*i).second + 1))
            {
                // delete?
#if defined(DEBUG_CIGAR)
                cout << "Delete";
#endif
                cigar.Add(CigarRoller::del, 1);
            }
        }
        // if there is soft clipping, allow for it (::Add will
        // ignore if the count is 0):
        cigar.Add(CigarRoller::softClip, getSoftClipCount());
#if defined(DEBUG_CIGAR)
        cout << endl;
#endif
    }

    //
    // Append cigar operations to an existing cigar list.
    //
    // XXX we no longer need the CigarRoller += methods.
    //
    // In this case, the Smith Waterman array H was created from
    // the read and reference in the reverse direction.
    //
    void rollCigarBackward(CigarRoller &cigar)
    {
        vector<pair<int,int> >::iterator i;

        // if there is soft clipping, allow for it (::Add will
        // ignore if the count is 0):
        cigar.Add(CigarRoller::softClip, getSoftClipCount());

        i = alignment.begin();

        for (i=alignment.begin();
                i < alignment.end() - 1;
                i++)
        {
#if defined(DEBUG_CIGAR)
            cout << *i << ": ";
#endif
            if ((*(i+1)).first == ((*i).first-1) && (*(i+1)).second == ((*i).second - 1))
            {
                // match/mismatch
#if defined(DEBUG_CIGAR)
                cout << "Match/Mismatch";
#endif
                cigar.Add(CigarRoller::match, 1);
            }
            else if ((*(i+1)).first == ((*i).first-1) && (*(i+1)).second == ((*i).second))
            {
                // insert?
#if defined(DEBUG_CIGAR)
                cout << "Insert?";
#endif
                cigar.Add(CigarRoller::insert, 1);
            }
            else if ((*(i+1)).first == ((*i).first) && (*(i+1)).second == ((*i).second - 1))
            {
                // delete?
#if defined(DEBUG_CIGAR)
                cout << "Delete?";
#endif
                cigar.Add(CigarRoller::del, 1);
            }
        }
#if defined(DEBUG_CIGAR)
        cout << endl;
#endif
    }

    //
    // Given the direction, and the alignment vector, obtain
    // the soft clip (the mismatches at the end of the string which
    // can in Smith Waterman matching be considered as a separate case).
    //
    // NB: be careful that the backward case is correct - it passes
    // all of two built in tests, but it may not be generally correct.
    //
    int getSoftClipCount()
    {
        if (direction>0)
        {
            // invariant: assert(maxCostPosition == alignment.front());
            return m - maxCostPosition.first;
        }
        else
        {
//          return alignment.back().first;  // nope, this always returns 0
            // XXX BE CAREFUL... not sure this is right, either.
//          return n - maxCostPosition.second;
            return m - maxCostPosition.first;
        }
    }

    void rollCigar(CigarRoller &cigar)
    {
        if (direction>0) rollCigarForward(cigar);
        else rollCigarBackward(cigar);
    }

    //
    // all in one local alignment:
    //
    // Steps:
    //   1 - do internal setup
    //   2 - populate H
    //   3 - create alignment vector (this chooses the best path)
    //   4 - compute sumQ
    //   5 - compute the cigar string
    //   6 - compute and update the softclip for the read
    //
    bool localAlignment(
        uint32_t bandSize,
        Atype &read,
        readIndexType readLength,
        QualityType &quality,
        Btype &reference,
        referenceIndexType referenceLength,
        referenceIndexType referenceOffset,
        CigarRoller &cigarRoller,
        uint32_t &softClipCount,
        referenceIndexType &cigarStartingPoint,
        int             &sumQ
    )
    {

        clear();

        cigarRoller.clear();

        setDirection(+1);
        setAllowedInsertDelete(bandSize);

        setRead(&read);
        setReadOffset(0);
        setReadLength(readLength);

        setReadQuality(&quality);

        setReference(&reference);
        setReferenceOffset(referenceOffset);
        setReferenceLength(referenceLength);

        populateH();

        softClipCount = getSoftClipCount();

        populateAlignment();

        rollCigar(cigarRoller);

        sumQ = getSumQ();

        return false;

    };

};

#endif
