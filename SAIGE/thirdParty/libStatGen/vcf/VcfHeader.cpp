/*
 *  Copyright (C) 2010-2011  Regents of the University of Michigan,
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

#include "VcfHeader.h"

VcfHeader::VcfHeader()
    : myHeaderLines()
{
    reset();
}


VcfHeader::~VcfHeader()
{
}


bool VcfHeader::read(IFILE filePtr)
{
    // Reading, so clean out this header.
    reset();

    if(filePtr == NULL)
    {
        // No file was passed in.
        myStatus.setStatus(StatGenStatus::FAIL_ORDER, 
                           "Need to pass in an open file ptr to VcfHeader::read.");
        return(false);
    }

    // Read until the header line has been read (after the meta lines).
    while(!myHasHeaderLine)
    {
        // Increase the size of headerlines by 1 to fit the new line.
        myHeaderLines.resize(myHeaderLines.size() + 1);
        
        // Read the next line from the file into the header structure.
        String& newStr = myHeaderLines.back();
        if(newStr.ReadLine(filePtr) < 0)
        {
            // Error, unable to read an entire header from the file.
            myStatus.setStatus(StatGenStatus::INVALID, 
                               "Error reading VCF Meta/Header, EOF found before the header line.");
            return(false);
        }
        if(newStr.Length() <= 2)
        {
            // A header/meta line must have at least 2 characters
            // ## or # and 8 fields, so if less than 2 characters,
            // error.
            myStatus.setStatus(StatGenStatus::INVALID, 
                               "Error reading VCF Meta/Header, line without at least 2 characters found before the header line.");
            return(false);
        }

        // Check if it is a header (first char is # and 2nd one is not).
        if((newStr[0] == '#') && (newStr[1] != '#'))
        {
            myHasHeaderLine = true;

            // Parse the header line to get the sample information.
            myParsedHeaderLine.ReplaceColumns(newStr, '\t');
        }
        else if((newStr[0] != '#') || (newStr[1] != '#'))
        {
            // A meta line must start with "##", we expect meta lines until
            // the header line is found.
            myStatus.setStatus(StatGenStatus::INVALID, 
                               "Error reading VCF Meta/Header, line not starting with '##' found before the header line.");
            return(false);
        }
    }
    return(true);
}


bool VcfHeader::write(IFILE filePtr)
{
    if(filePtr == NULL)
    {
        // No file was passed in.
        myStatus.setStatus(StatGenStatus::FAIL_ORDER, 
                           "Need to pass in an open file ptr to VcfHeader::write.");
        return(false);
    }
    
    // Make sure the last header line is synced with the parsed header line.
    syncHeaderLine();

    int numWritten = 0;
    int numExpected = 0;
    for(std::vector<String>::iterator iter = myHeaderLines.begin(); 
        iter != myHeaderLines.end(); iter++)
    {
        numWritten += ifprintf(filePtr, "%s\n", iter->c_str());
        // expected to write string + new line.
        numExpected += iter->Length();
        numExpected += 1;
    }
    if(numWritten != numExpected)
    {
        myStatus.setStatus(StatGenStatus::FAIL_IO, 
                           "Failed writing VCF Meta/Header.");
    }
    return(numWritten == numExpected);
}


void VcfHeader::reset()
{
    myHasHeaderLine = false;
    myHeaderLines.clear();
}


// Return the error after a failed call.
const StatGenStatus& VcfHeader::getStatus()
{
    return(myStatus);
}


int VcfHeader::getNumMetaLines()
{
    int numHeaderLines = myHeaderLines.size();
    if((numHeaderLines >= 1) && (myHasHeaderLine))
    {
        // Remove the header line from the count.
        return(numHeaderLines-1);
    }
    return(numHeaderLines);
}


const char* VcfHeader::getMetaLine(unsigned int index)
{
    if(index >= myHeaderLines.size())
    {
        return(NULL);
    }
    else
    {
        return(myHeaderLines[index].c_str());
    }
    return(NULL);
}


const char* VcfHeader::getHeaderLine()
{
    // Make sure the last header line is synced with the parsed header line.
    syncHeaderLine();
    if(myHasHeaderLine)
    {
        return(myHeaderLines.back().c_str());
    }
    return(NULL);
}


int VcfHeader::getNumSamples() const
{
    if(!myHasHeaderLine)
    {
        return(0);
    }
    
    int numFields = myParsedHeaderLine.Length();

    if(numFields > NUM_NON_SAMPLE_HEADER_COLS)
    {
        // There are samples.
        return(numFields - NUM_NON_SAMPLE_HEADER_COLS);
    }

    // No sample fields
    return(0);
}


const char* VcfHeader::getSampleName(unsigned int index) const
{
    if(!myHasHeaderLine)
    {
        // No header.
        return(NULL);
    }
    int position = index + NUM_NON_SAMPLE_HEADER_COLS;

    if(position >= myParsedHeaderLine.Length())
    {
        // Out of range.
        return(NULL);
    }

    return(myParsedHeaderLine[position].c_str());
}


int VcfHeader::getSampleIndex(const char* sampleName) const
{
    if(!myHasHeaderLine)
    {
        // No header.
        return(-1);
    }
    for(int index = NUM_NON_SAMPLE_HEADER_COLS; 
        index < myParsedHeaderLine.Length(); index++)
    {
        if(myParsedHeaderLine[index] == sampleName)
        {
            // Found.
            return(index - NUM_NON_SAMPLE_HEADER_COLS);
        }
    }
    // Not found.
    return(-1);
}


void VcfHeader::removeSample(unsigned int index)
{
    int position = index + NUM_NON_SAMPLE_HEADER_COLS;

    if(position >= myParsedHeaderLine.Length())
    {
        // Out of range, so just return, nothing to remove.
        return;
    }

    // Remove it from the parsed header line.
    myParsedHeaderLine.Delete(position);

    // Removed a sample, so clear the header line so the next time it is
    // accessed it will be reset based on the existing samples.
    String& hdrLine = myHeaderLines.back();
    hdrLine.Clear();
}


bool VcfHeader::appendMetaLine(const char* metaLine)
{
    // Check that the line starts with "##".
    if(strncmp(metaLine, "##", 2) != 0)
    {
        // Does not start with "##"
        return(false);
    }
    if(!myHasHeaderLine)
    {
        // No header line, so just add to the end of the vector.
        myHeaderLines.push_back(metaLine);
        return(true);
    }
    // There is a header line, so insert this just before that line.
    // The headerLine is one position before "end".
    std::vector<String>::iterator headerLine = myHeaderLines.end();
    --headerLine;
    // Insert just before the header line.
    myHeaderLines.insert(headerLine, metaLine);
    return(true);
}


bool VcfHeader::addHeaderLine(const char* headerLine)
{
    // Check that the line starts with "#".
    if(strncmp(headerLine, "#", 1) != 0)
    {
        // Does not start with "#"
        return(false);
    }

    if(myHasHeaderLine)
    {
        // There is a header line, so replace the current line.
        myHeaderLines.back() = headerLine;
    }
    else
    {
        // There is not a header line, so add it 
        myHeaderLines.push_back(headerLine);
    }

    myHasHeaderLine = true;
    // Parse the header line to get the sample information.
    myParsedHeaderLine.ReplaceColumns(headerLine, '\t');
    return(true);
}


void VcfHeader::syncHeaderLine()
{
    if(!myHasHeaderLine)
    {
        // No header line, so nothing to sync.
        return;
    }
    // Get the last header line and see if it is set.
    String& hdrLine = myHeaderLines.back();

    if(hdrLine.IsEmpty())
    {
        // The header line is not set, so set it.
        for(int i  = 0; i < myParsedHeaderLine.Length(); i++)
        {
            if(i != 0)
            {
                hdrLine += '\t';
            }
            hdrLine += myParsedHeaderLine[i];
        }
    }
}
