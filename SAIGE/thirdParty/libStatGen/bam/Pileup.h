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

#ifndef __PILEUP_H__
#define __PILEUP_H__

#include <stdexcept>
#include "SamFile.h"
#include "PosList.h"


class PileupHelper
{
public:
    static const int DEFAULT_WINDOW_SIZE = 1024;
};

template <class PILEUP_TYPE>
class defaultPileup
{
public:
    bool operator() (PILEUP_TYPE& element)
    {
        element.analyze();
        return(true);
    }
    void analyze(PILEUP_TYPE element)
    {
        element.analyze();
    }
};

template <class PILEUP_TYPE>
void defaultPileupAnalyze(PILEUP_TYPE& element)
{
    element.analyze();
}


/// Class to perform a pileup of all reads by position, assuming
/// the reads are coordinate sorted.
template <class PILEUP_TYPE, 
          class FUNC_CLASS = defaultPileup<PILEUP_TYPE> >
class Pileup
{
public:
    /// Constructor using the default maximum number of bases a read spans.
    Pileup(const FUNC_CLASS& fp = FUNC_CLASS());

    /// Constructor that sets the maximum number of bases a read spans.
    /// This is the "window" the length of the buffer that holds
    /// the pileups for each position until the read start has moved
    /// past the position.
    Pileup(int window, const FUNC_CLASS& fp = FUNC_CLASS());

    /// Perform pileup with a reference.
    Pileup(const std::string& refSeqFileName,
           const FUNC_CLASS& fp = FUNC_CLASS());

    /// Perform pileup with a reference and a specified window size.
    Pileup(int window, const std::string& refSeqFileName, 
           const FUNC_CLASS& fp = FUNC_CLASS());

    /// Destructor
    virtual ~Pileup();

    /// Performs a pileup on the specified file.
    /// \param excludeFlag if specified, if any bit set in the exclude flag
    ///                    is set in the record's flag, it will be dropped.
    /// Defaulted to exclude:
    ///    * unmapped, 
    ///    * not primary alignment
    ///    * failed platform/vendor quality check
    ///    * PCR or optical duplicate
    /// \param includeFlag if specified, every bit must be set in the 
    ///                    record's flag for it to be included - 
    ///                    defaulted to 0, no bits are required to be set.
    /// \return 0 for success and non-zero for failure.
    virtual int processFile(const std::string& fileName, 
                            uint16_t excludeFlag = 0x0704, 
                            uint16_t includeFlag = 0);

    /// Add an alignment to the pileup.
    virtual void processAlignment(SamRecord& record);
   
    /// Add only positions that fall within the specified region of the
    /// alignment to the pileup and outside of the specified excluded positions.
    /// \param record alignment to be added to the pileup.
    /// \param startPos 0-based start position of the bases that should be
    ///                 added to the pileup.
    /// \param endPos   0-based end position of the bases that should be added
    ///                 to the pileup (this position is not added).
    ///                 Set to -1 if there is no end position to the region.
    /// \param excludeList list of refID/positions to exclude from processing.
    virtual void processAlignmentRegion(SamRecord& record,
                                        int startPos,
                                        int endPos,
                                        PosList* excludeList = NULL);
   
    /// Done processing, flush every position that is currently being stored
    /// in the pileup.
    void flushPileup();

protected:
    FUNC_CLASS myAnalyzeFuncPtr;

    // Always need the reference position.
    void addAlignmentPosition(int refPosition, SamRecord& record);


    virtual void flushPileup(int refID, int refPosition);
    void flushPileup(int refPosition);
    
    // Get the position in the myElements container that is associated
    // with the specified position.  If the specified position cannot
    // fit within the myElements container, -1 is returned.
    int pileupPosition(int refPosition);

    virtual void resetElement(PILEUP_TYPE& element, int position);
    virtual void addElement(PILEUP_TYPE& element, SamRecord& record);
    virtual void analyzeElement(PILEUP_TYPE& element);
    virtual void analyzeHead();

    std::vector<PILEUP_TYPE> myElements;

    int    pileupStart;
    int    pileupHead;
    int    pileupTail; // last used position
    int    pileupWindow;

    int myCurrentRefID;

    GenomeSequence* myRefPtr;
};


template <class PILEUP_TYPE, class FUNC_CLASS>
Pileup<PILEUP_TYPE, FUNC_CLASS>::Pileup(const FUNC_CLASS& fp)
    : myAnalyzeFuncPtr(fp),
      myElements(),
      pileupStart(0),
      pileupHead(0),
      pileupTail(-1),
      pileupWindow(PileupHelper::DEFAULT_WINDOW_SIZE),
      myCurrentRefID(-2),
      myRefPtr(NULL)
{
    // Not using pointers since this is templated.
    myElements.resize(pileupWindow);
}


template <class PILEUP_TYPE, class FUNC_CLASS>
Pileup<PILEUP_TYPE, FUNC_CLASS>::Pileup(int window, const FUNC_CLASS& fp)
    : myAnalyzeFuncPtr(fp),
      myElements(),
      pileupStart(0),
      pileupHead(0),
      pileupTail(-1),
      pileupWindow(window),
      myCurrentRefID(-2),
      myRefPtr(NULL)
{
    // Not using pointers since this is templated.
    myElements.resize(window);
}


template <class PILEUP_TYPE, class FUNC_CLASS>
Pileup<PILEUP_TYPE, FUNC_CLASS>::Pileup(const std::string& refSeqFileName, const FUNC_CLASS& fp)
    : myAnalyzeFuncPtr(fp),
      myElements(),
      pileupStart(0),
      pileupHead(0),
      pileupTail(-1),
      pileupWindow(PileupHelper::DEFAULT_WINDOW_SIZE),
      myCurrentRefID(-2),
      myRefPtr(NULL)
{
    myRefPtr = new GenomeSequence(refSeqFileName.c_str());

    // Not using pointers since this is templated.
    myElements.resize(pileupWindow);

    PILEUP_TYPE::setReference(myRefPtr);
}


template <class PILEUP_TYPE, class FUNC_CLASS>
Pileup<PILEUP_TYPE, FUNC_CLASS>::Pileup(int window, const std::string& refSeqFileName, const FUNC_CLASS& fp)
    : myAnalyzeFuncPtr(fp),
      myElements(),
      pileupStart(0),
      pileupHead(0),
      pileupTail(-1),
      pileupWindow(window),
      myCurrentRefID(-2),
      myRefPtr(NULL)
{
    myRefPtr = new GenomeSequence(refSeqFileName.c_str());

    // Not using pointers since this is templated.
    myElements.resize(window);

    PILEUP_TYPE::setReference(myRefPtr);
}


template <class PILEUP_TYPE, class FUNC_CLASS>
Pileup<PILEUP_TYPE, FUNC_CLASS>::~Pileup()
{
    flushPileup();
    if(myRefPtr != NULL)
    {
        delete myRefPtr;
        myRefPtr = NULL;
    }
}

template <class PILEUP_TYPE, class FUNC_CLASS>
int Pileup<PILEUP_TYPE, FUNC_CLASS>::processFile(const std::string& fileName, 
                                                 uint16_t excludeFlag,
                                                 uint16_t includeFlag)
{
    SamFile samIn;
    SamFileHeader header;
    SamRecord record;
    
    if(myRefPtr != NULL)
    {
        samIn.SetReference(myRefPtr);
    }

    if(!samIn.OpenForRead(fileName.c_str()))
    {
        fprintf(stderr, "%s\n", samIn.GetStatusMessage());
        return(samIn.GetStatus());
    }
    
    if(!samIn.ReadHeader(header))
    {
        fprintf(stderr, "%s\n", samIn.GetStatusMessage());
        return(samIn.GetStatus());
    }

    // The file needs to be sorted by coordinate.
    samIn.setSortedValidation(SamFile::COORDINATE);

    // Iterate over all records
    while (samIn.ReadRecord(header, record))
    {
        uint16_t flag = record.getFlag();
        if(flag & excludeFlag)
        {
            // This record has an excluded flag set, 
            // so continue to the next one.
            continue;
        }
        if((flag & includeFlag) != includeFlag)
        {
            // This record does not have all required flags set, 
            // so continue to the next one.
            continue;
        }
        processAlignment(record);
    }

    flushPileup();

    int returnValue = 0;
    if(samIn.GetStatus() != SamStatus::NO_MORE_RECS)
    {
        // Failed to read a record.
        fprintf(stderr, "%s\n", samIn.GetStatusMessage());
        returnValue = samIn.GetStatus();
    }
    return(returnValue);  
}


template <class PILEUP_TYPE, class FUNC_CLASS>
void Pileup<PILEUP_TYPE, FUNC_CLASS>::processAlignment(SamRecord& record)
{
    int refPosition = record.get0BasedPosition();
    int refID = record.getReferenceID();

    // Flush any elements from the pileup that are prior to this record
    // since the file is sorted, we are done with those positions.
    flushPileup(refID, refPosition);
    
    // Loop through for each reference position covered by the record.
    // It is up to the PILEUP_TYPE to handle insertions/deletions, etc
    // that are related with the given reference position.
    for(; refPosition <= record.get0BasedAlignmentEnd(); ++refPosition)
    {
        addAlignmentPosition(refPosition, record);
    }
}


template <class PILEUP_TYPE, class FUNC_CLASS>
void Pileup<PILEUP_TYPE, FUNC_CLASS>::processAlignmentRegion(SamRecord& record,
                                                             int startPos,
                                                             int endPos,
                                                             PosList* excludeList)
{
    int refPosition = record.get0BasedPosition();
    int refID = record.getReferenceID();

    // Flush any elements from the pileup that are prior to this record
    // since the file is sorted, we are done with those positions.
    flushPileup(refID, refPosition);
    
    // Check if the region starts after this reference starts.  If so,
    // we only want to start adding at the region start position.
    if(startPos > refPosition)
    {
        refPosition = startPos;
    }

    // Loop through for each reference position covered by the record.
    // It is up to the PILEUP_TYPE to handle insertions/deletions, etc
    // that are related with the given reference position.
    for(; refPosition <= record.get0BasedAlignmentEnd(); ++refPosition)
    {
        // Check to see if we have gone past the end of the region, in which
        // case we can stop processing this record.  Check >= since the
        // end position is not in the region.
        if((endPos != -1) && (refPosition >= endPos))
        {
            break;
        }

        // Check to see if this position is in the exclude list.
        bool addPos = true;
        if(excludeList != NULL)
        {
            // There is an exclude list, so lookup the position.
            if(excludeList->hasPosition(refID, refPosition))
            {
                // This position is in the exclude list, so don't add it.
                addPos = false;
            }
        }
        if(addPos)
        {
            addAlignmentPosition(refPosition, record);
        }
    }
}


template <class PILEUP_TYPE, class FUNC_CLASS>
void Pileup<PILEUP_TYPE, FUNC_CLASS>::flushPileup()
{
    // while there are still entries between the head and tail, flush,
    // but no need to flush if pileupTail == -1 because in that case 
    // no entries have been added
    while ((pileupHead <= pileupTail) && (pileupTail != -1))
    {
        flushPileup(pileupHead+1);
    }
    pileupStart = pileupHead = 0;
    pileupTail = -1;
}


// Always need the reference position.
template <class PILEUP_TYPE, class FUNC_CLASS>
void Pileup<PILEUP_TYPE, FUNC_CLASS>::addAlignmentPosition(int refPosition,
                                                           SamRecord& record)
{
    int offset = 0;
    try{
        offset = pileupPosition(refPosition);
    }
    catch(std::runtime_error& err)
    {
        const char* overflowErr = "Overflow on the pileup buffer:";
        String errorMessage = err.what();
        if(strncmp(err.what(), overflowErr, strlen(overflowErr)) == 0)
        {
            
            errorMessage += "\n\tPileup Buffer Overflow: recordName = ";
            errorMessage += record.getReadName();
            errorMessage += "; Cigar = ";
            errorMessage += record.getCigar();
        }
        throw std::runtime_error(errorMessage.c_str());
    }
    
    if((offset < 0) || (offset >= pileupWindow))
    {
        std::cerr << "Pileup Buffer Overflow: position = " << refPosition
                  << "; refID = " << record.getReferenceID() 
                  << "; recStartPos = " << record.get1BasedPosition()
                  << "; pileupStart = " << pileupStart
                  << "; pileupHead = " << pileupHead
                  << "; pileupTail = " << pileupTail;
    }

    addElement(myElements[offset], record);
}


template <class PILEUP_TYPE, class FUNC_CLASS>
void Pileup<PILEUP_TYPE, FUNC_CLASS>::flushPileup(int refID, int position)
{
    // if new chromosome, flush the entire pileup.
    if(refID != myCurrentRefID)
    {
        // New chromosome, flush everything.
        flushPileup();
        myCurrentRefID = refID;
    }
    else
    {
        // on the same chromosome, so flush just up to this new position.
        flushPileup(position);
    }
}


template <class PILEUP_TYPE, class FUNC_CLASS>
void Pileup<PILEUP_TYPE, FUNC_CLASS>::flushPileup(int position)
{
    // Flush up to this new position, but no reason to flush if
    // pileupHead has not been set.
    while((pileupHead < position) && (pileupHead <= pileupTail))
    {
        analyzeHead();

        pileupHead++;
        
        if(pileupHead - pileupStart >= pileupWindow)
            pileupStart += pileupWindow;
    }

    if(pileupHead > pileupTail)
    {
        // All positions have been flushed, so reset pileup info
        pileupHead = pileupStart = 0;
        pileupTail = -1;
    }
}


// Get the position in the myElements container that is associated
// with the specified position.  If the specified position cannot
// fit within the myElements container, -1 is returned.
template <class PILEUP_TYPE, class FUNC_CLASS>
int Pileup<PILEUP_TYPE, FUNC_CLASS>::pileupPosition(int position)
{
    // Check to see if this is the first position (pileupTail == -1)
    if(pileupTail == -1)
    {
        pileupStart = pileupHead = position;
        // This is the first time this position is being used, so
        // reset the element.
        resetElement(myElements[0], position);
        pileupTail = position;
        return(0);
    }


    if((position < pileupHead) || (position > (pileupHead + pileupWindow)))
    {
        String errorMessage =
            "Overflow on the pileup buffer: specifiedPosition = ";
        errorMessage += position;
        errorMessage += ", pileup buffer start position: ";
        errorMessage += pileupHead;
        errorMessage += ", pileup buffer end position: ";
        errorMessage += pileupHead + pileupWindow;

        throw std::runtime_error(errorMessage.c_str());
    }    

    //   int offset = position - pileupStart;
    int offset = position - pileupStart;
    
    if(offset >= pileupWindow)
    {
        offset -= pileupWindow;
    }

    // Check to see if position is past the end of the currently
    // setup pileup positions.
    while(position > pileupTail)
    {
        // Increment pileupTail to the next position since the current
        // pileupTail is already in use.
        ++pileupTail;

        // Figure out the offset for this next position.
        offset = pileupTail - pileupStart;
        if(offset >= pileupWindow)
        {
            offset -= pileupWindow;
        }

        // This is the first time this position is being used, so
        // reset the element.
        resetElement(myElements[offset], pileupTail);
    }

    return(offset);
}


template <class PILEUP_TYPE, class FUNC_CLASS>
void Pileup<PILEUP_TYPE, FUNC_CLASS>::resetElement(PILEUP_TYPE& element,
                                                   int position)
{
    element.reset(position);
}


template <class PILEUP_TYPE, class FUNC_CLASS>
void Pileup<PILEUP_TYPE, FUNC_CLASS>::addElement(PILEUP_TYPE& element,
                                                 SamRecord& record)
{
    element.addEntry(record);
}


template <class PILEUP_TYPE, class FUNC_CLASS>
void Pileup<PILEUP_TYPE, FUNC_CLASS>::analyzeElement(PILEUP_TYPE& element)
{
    myAnalyzeFuncPtr(element);
}


template <class PILEUP_TYPE, class FUNC_CLASS>
void Pileup<PILEUP_TYPE, FUNC_CLASS>::analyzeHead()
{
    myAnalyzeFuncPtr(myElements[pileupHead - pileupStart]);
}


#endif
