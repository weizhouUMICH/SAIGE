#include <iostream>
#include <time.h>
#include "GenomeSequence.h"
#include "InputFile.h"

void readDbsnp(mmapArrayBool_t& dbSNP, const char* fileName, GenomeSequence& ref);

int main(int argc, char ** argv)
{
    //    time_t startTime;
    //    time_t endTime;
    //    startTime = time(NULL);
    GenomeSequence* refPtr = new GenomeSequence("testFiles/chr1_partial.fa");
    //    endTime = time(NULL);

    //    std::cerr << "Time to read reference: " << endTime - startTime << std::endl;
    if(refPtr == NULL)
    {
        std::cerr << "Failed to read the reference\n";
        return(-1);
    }
    std::cerr << "\nStandard VCF DBSNP test\n";
    mmapArrayBool_t dbsnpArray1;
    const char* dbsnpFileName = "testFiles/dbsnp.vcf";

    //    startTime = time(NULL);
    refPtr->loadDBSNP(dbsnpArray1, dbsnpFileName);
    //    endTime = time(NULL);
    //    std::cerr << "Time to read dbsnp through reference: " << endTime - startTime << std::endl;


    genomeIndex_t mapPos = 
        refPtr->getGenomePosition("1", 10233);
    std::cerr << "dbsnp " << mapPos   << ": " 
              << dbsnpArray1[mapPos] << std::endl;
    std::cerr << "dbsnp " << mapPos+1 << ": " 
              << dbsnpArray1[mapPos+1] << std::endl;
    std::cerr << "dbsnp " << mapPos+2 << ": " 
              << dbsnpArray1[mapPos+2] << std::endl;


    std::cerr << "\nGZIP VCF DBSNP test\n";

    mmapArrayBool_t dbsnpArray2;
    dbsnpFileName = "testFiles/dbsnp.vcf.gz";

    //    startTime = time(NULL);
    refPtr->loadDBSNP(dbsnpArray2, dbsnpFileName);
    //    endTime = time(NULL);
    //    std::cerr << "Time to read dbsnp through reference: " << endTime - startTime << std::endl;


    mapPos = refPtr->getGenomePosition("1", 10233);
    std::cerr << "dbsnp " << mapPos   << ": " 
              << dbsnpArray2[mapPos] << std::endl;
    std::cerr << "dbsnp " << mapPos+1 << ": " 
              << dbsnpArray2[mapPos+1] << std::endl;
    std::cerr << "dbsnp " << mapPos+2 << ": " 
              << dbsnpArray2[mapPos+2] << std::endl;
    return(0);
}
