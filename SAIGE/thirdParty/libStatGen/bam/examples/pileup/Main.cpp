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
#include "SamFile.h"
#include "Pileup.h"
#include  "PileupElementBaseQual.h"


void newAnalyze(PileupElementBaseQual& element)
{
    std::cout << "newAnalyze: ";
    element.analyze();
}


class AnalyzeClass
{
public:
    AnalyzeClass()
    {
        myCounter = 33;
    }
    bool operator() (PileupElementBaseQual& element)
    {
        std::cout << "Class Analyze: Counter = " << myCounter << ": ";
        element.analyze();
        ++myCounter;
        return(true);
    }
    int myCounter;
private:
};


int main(int argc, char ** argv)
{
    const char* fileName = "../../test/testFiles/sortedBam.bam";
    const char* indexName = "../../test/testFiles/sortedBam.bam.bai";

    printf("\nPileup<PileupElementBaseQual> on entire file: %s\n", fileName);
    Pileup<PileupElementBaseQual> pileup(1024);
    pileup.processFile(fileName);

    printf("\nPileup<PileupElement> on entire file: %s\n", fileName);
    Pileup<PileupElement> pileup1(1024);
    pileup1.processFile(fileName);

    printf("\nPileup<PileupElementBaseQual> on a section of file: %s\n", fileName);
    // Read a sorted & indexed BAM file.
    Pileup<PileupElementBaseQual> pileup2(100);
    
    SamFile samIn;
    SamFileHeader header;
    SamRecord record;
    
    if(!samIn.OpenForRead(fileName))
    {
        fprintf(stderr, "%s\n", samIn.GetStatusMessage());
        return(samIn.GetStatus());
    }
    
    // Open the bam index file for reading.
    if(!samIn.ReadBamIndex(indexName))
    {
        fprintf(stderr, "%s\n", samIn.GetStatusMessage());
        return(samIn.GetStatus());
    }

    if(!samIn.ReadHeader(header))
    {
        fprintf(stderr, "%s\n", samIn.GetStatusMessage());
        return(samIn.GetStatus());
    }

    const char* refName = "1";
    int start = 1000;
    int end = 1500;
    if(!samIn.SetReadSection(refName, start, end))
    {
        fprintf(stderr, "%s\n", samIn.GetStatusMessage());
        return(samIn.GetStatus());
    }

    // Iterate over all records
    while (samIn.ReadRecord(header, record))
    {
        pileup2.processAlignment(record);
    }

    pileup2.flushPileup();

    int returnValue = 0;
    if(samIn.GetStatus() != SamStatus::NO_MORE_RECS)
    {
        // Failed to read a record.
        fprintf(stderr, "%s\n", samIn.GetStatusMessage());
        returnValue = samIn.GetStatus();
    }

    printf("\nPileup<PileupElementBaseQual> on entire file, newAnalyze: %s\n", fileName);

    void (*fnPtr)(PileupElementBaseQual&) = newAnalyze;

    Pileup<PileupElementBaseQual, void (*)(PileupElementBaseQual&)> pileup3(1024, fnPtr);
    pileup3.processFile(fileName);

    printf("\nPileup<PileupElementBaseQual> on entire file, newAnalyze: %s\n", fileName);
    
    AnalyzeClass myAnalyzeClass;
    myAnalyzeClass.myCounter = 2;
    Pileup<PileupElementBaseQual, AnalyzeClass> pileup4(1024, myAnalyzeClass);
    pileup4.processFile(fileName);

    return(0);
}
