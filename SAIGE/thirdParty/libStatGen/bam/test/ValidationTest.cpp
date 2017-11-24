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

#include "SamRecord.h"
#include "SamValidation.h"
#include "ValidationTest.h"
#include <assert.h>

void testSamQNAME()
{
    // This method tests:
    //   QNAME.Length() > 0 and <= 254
    //   QNAME does not contain [ \t\n\r]

    char qname[256];
    SamFileHeader samHeader;
    SamRecord testRecord(ErrorHandler::RETURN);
    // Error list
    SamValidationErrors errorList;
   
    // Test Length == 0 by setting qname[0] to 0 (end of char*)
    qname[0] = 0;
    // It fails, because it is a required field.
    assert(testRecord.setReadName(qname) == false);
    assert(strcmp(testRecord.getReadName(), "UNKNOWN") == 0);
    // It was reset to the default which is valid.
    assert(SamValidator::isValid(samHeader, testRecord, errorList) == true);
    assert(errorList.numErrors() == 0);
    assert(errorList.getNextError() == NULL);

    // Test too long of a read name.
    memset(qname, '.', 255);
    qname[255] = 0;
    assert(testRecord.setReadName(qname) == true);
    assert(strcmp(testRecord.getReadName(), qname) == 0);
    assert(SamValidator::isValid(samHeader, testRecord, errorList) == false);
    // 2 errors - 1 since the qname is longer than 254 (it is 255).
    // and the qname length including the null is 256, but the 
    // read name length is only 8 bits, so that is a 1.
    assert(errorList.numErrors() == 2);
    assert(errorList.getNextError()->getType() == 
           SamValidationError::INVALID_QNAME);
    assert(errorList.getNextError()->getType() == 
           SamValidationError::INVALID_QNAME);
    assert(errorList.getNextError() == NULL);

    // Clear the error list
    errorList.clear();

    // Setup a buffer to set the record to.
    int bufferBlockSize = 32;

    bamRecordStruct* bufferRecordPtr =
        (bamRecordStruct *) malloc(bufferBlockSize + sizeof(int));

    bufferRecordPtr->myBlockSize = bufferBlockSize;
    bufferRecordPtr->myReferenceID = -1;
    bufferRecordPtr->myPosition = 1010;
    // Set the read name length to 0.
    bufferRecordPtr->myReadNameLength = 0;
    bufferRecordPtr->myMapQuality = 0;
    bufferRecordPtr->myBin = 4681;
    bufferRecordPtr->myCigarLength = 0;
    bufferRecordPtr->myFlag = 73;
    bufferRecordPtr->myReadLength = 0;
    bufferRecordPtr->myMateReferenceID = -1;
    bufferRecordPtr->myMatePosition = 1010;
    bufferRecordPtr->myInsertSize = 0;

    assert(testRecord.setBuffer((const char*)bufferRecordPtr, 
                                bufferBlockSize + sizeof(int), 
                                samHeader) == SamStatus::SUCCESS);
    // 1 error - the read name length is 0.
    assert(SamValidator::isValid(samHeader, testRecord, errorList) == false);
    assert(errorList.numErrors() == 1);
    assert(errorList.getNextError()->getType() == 
           SamValidationError::INVALID_QNAME);
    assert(errorList.getNextError() == NULL);

    // Clear the error list
    errorList.clear();

    // Test a buffer that has a read name, but the length specified is
    // longer than the first null.
    bufferBlockSize = 40;
    bufferRecordPtr->myBlockSize = bufferBlockSize;
    // Set the read name length to 8 - longer than 3 - "HI\0".
    bufferRecordPtr->myReadNameLength = 8;
    bufferRecordPtr->myData[0] = 'H';
    bufferRecordPtr->myData[1] = 'I';
    bufferRecordPtr->myData[2] = 0;

    assert(testRecord.setBuffer((const char*)bufferRecordPtr,
                                bufferBlockSize + sizeof(int), 
                                samHeader) == SamStatus::SUCCESS);
    // 1 error - the read name length in the buffer does not match the
    // length of the read name to the first null.
    assert(SamValidator::isValid(samHeader, testRecord, errorList) == false);
    assert(errorList.numErrors() == 1);
    assert(errorList.getNextError()->getType() ==
           SamValidationError::INVALID_QNAME);
    assert(errorList.getNextError() == NULL);

    // Clear the error list
    errorList.clear();

    // Test a buffer that has a read name, but the length specified is
    // shorter than the first null.
    bufferBlockSize = 34;
    bufferRecordPtr->myBlockSize = bufferBlockSize;
    // Set the read name length to 2 - longer than 3 - "HI\0"..
    bufferRecordPtr->myReadNameLength = 2;
    bufferRecordPtr->myData[0] = 'H';
    bufferRecordPtr->myData[1] = 'I';
    bufferRecordPtr->myData[2] = 0;

    assert(testRecord.setBuffer((const char*)bufferRecordPtr, 
                                bufferBlockSize + sizeof(int), 
                                samHeader) == SamStatus::SUCCESS);
    // 1 error - the read name length in the buffer does not match
    // the length of the read name to the first null.
    assert(SamValidator::isValid(samHeader, testRecord, errorList) == false);
    assert(errorList.numErrors() == 1);
    assert(errorList.getNextError()->getType() ==
           SamValidationError::INVALID_QNAME);
    assert(errorList.getNextError() == NULL);

    // Clear the error list
    errorList.clear();
}


void testBamRID()
{
    // BAM
    SamRecord testRecord(ErrorHandler::RETURN);
    // Error list
    SamValidationErrors errorList;
    SamFileHeader samHeader;

    // Clear the error list
    errorList.clear();

    // Setup a buffer to set the record to.
    int bufferBlockSize = 35;

    bamRecordStruct* bufferRecordPtr =
        (bamRecordStruct *) malloc(bufferBlockSize + sizeof(int));

    bufferRecordPtr->myBlockSize = bufferBlockSize;
    bufferRecordPtr->myPosition = 1010;
    bufferRecordPtr->myReferenceID = -1;
    // Set the read name length to 0.
    bufferRecordPtr->myReadNameLength = 3;
    bufferRecordPtr->myMapQuality = 0;
    bufferRecordPtr->myBin = 4681;
    bufferRecordPtr->myCigarLength = 0;
    bufferRecordPtr->myFlag = 73;
    bufferRecordPtr->myReadLength = 0;
    bufferRecordPtr->myMateReferenceID = -1;
    bufferRecordPtr->myMatePosition = 1010;
    bufferRecordPtr->myInsertSize = 0;
    bufferRecordPtr->myData[0] = 'H';
    bufferRecordPtr->myData[1] = 'I';
    bufferRecordPtr->myData[2] = 0;

    ////////////////////////////////////////////
    // Test out of range reference sequence id.
    bufferRecordPtr->myReferenceID = 100;

    assert(testRecord.setBuffer((const char*)bufferRecordPtr, 
                                bufferBlockSize + sizeof(int),
                                samHeader) == SamStatus::SUCCESS);
    // 1 error - the read name length is 0.
    assert(SamValidator::isValid(samHeader, testRecord, errorList) == false);
    assert(errorList.numErrors() == 1);
    assert(errorList.getNextError()->getType() == 
           SamValidationError::INVALID_REF_ID);
    assert(errorList.getNextError() == NULL);

    // Clear the error list
    errorList.clear();

    ////////////////////////////////////////////
    // Test out of range reference sequence id.
    bufferRecordPtr->myReferenceID = -100;

    assert(testRecord.setBuffer((const char*)bufferRecordPtr,
                                bufferBlockSize + sizeof(int),
                                samHeader) == SamStatus::SUCCESS);
    // 1 error - the read name length is 0.
    assert(SamValidator::isValid(samHeader, testRecord, errorList) == false);
    assert(errorList.numErrors() == 1);
    assert(errorList.getNextError()->getType() == 
           SamValidationError::INVALID_REF_ID);
    assert(errorList.getNextError() == NULL);

    // Clear the error list
    errorList.clear();
}


void testEmptyQual()
{
   
}

