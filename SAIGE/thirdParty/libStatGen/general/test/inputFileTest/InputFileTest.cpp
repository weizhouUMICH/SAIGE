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
#include "InputFileTest.h"
#include <assert.h>
#include <iostream>
#include "StringBasics.h"

void testAdditional(const char *extension);
void testWrite();


int main(int argc, char ** argv)
{
   IFILE_Test myFile;

   myFile.test();

   testWrite();

   std::cout << "\nAdditional Tests: " << std::endl;

   testAdditional("txt");
#ifdef __ZLIB_AVAILABLE__
   testAdditional("gz");
#endif
}


const int IFILE_Test::TEST_FILE_SIZE = 37;
const int IFILE_Test::BGZF_TEST_FILE_SIZE = 93;
const std::string IFILE_Test::TEST_FILE_CONTENTS = "ABCDabcd1234\nEFGefg567\nhijklHIJKL8910";

void IFILE_Test::test()
{
   std::cout << "\nUncompressedFileType Tests:" << std::endl;
   testAll("txt");

#ifdef __ZLIB_AVAILABLE__
   std::cout << "\nGzipFileType Tests:" << std::endl;
   testAll("gz");

   std::cout << "\nBgzfFileType Tests:" << std::endl;
   testAll("bam");

   std::cout << "\n.glf file Tests:" << std::endl;
   testAll("glf");
#endif
}


void IFILE_Test::testAll(const char* extension)
{
    test_readFromFile(extension);
    test_readTilChar(extension);
    test_ifeof_ifrewind(extension);
    test_ifread_ifgetc(extension);
    test_ifclose(extension);
    test_ifseek(extension);
    test_noExistRead(extension);
}


void IFILE_Test::test_readFromFile(const char* extension)
{
   // First open the test file.
   openFile(extension);

   // Verify the file successfully opened.
   assert(myFileTypePtr != NULL);
   assert(isOpen());
   assert(myFileTypePtr->isOpen());

   // Track how many bytes are read by each call.
   int numBytesRead = 0;

   // Track the total number of the bytes that have been read from the file
   // at any given point.
   int totalBytesPreviouslyRead = 0;

   // Test readFromFile.
   numBytesRead = readFromFile(myTestBuffer, 4);
   assert(myTestBuffer[0] == TEST_FILE_CONTENTS[0]);
   assert(myTestBuffer[1] == TEST_FILE_CONTENTS[1]);
   assert(myTestBuffer[2] == TEST_FILE_CONTENTS[2]);
   assert(myTestBuffer[3] == TEST_FILE_CONTENTS[3]);
   assert(numBytesRead == 4);
   totalBytesPreviouslyRead += numBytesRead;
   // This read should not have affected the internal buffer.
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   // Should not be at eof
   assert(myFileTypePtr->eof() == false);
   assert(ifeof() == false);

   // Read again to verify that the next characters could be read.
   numBytesRead = readFromFile(myTestBuffer, 2);
   // Read 2 more characters from the test file.
   assert(myTestBuffer[0] == TEST_FILE_CONTENTS[4]);
   assert(myTestBuffer[1] == TEST_FILE_CONTENTS[5]);
   assert(myTestBuffer[2] == TEST_FILE_CONTENTS[2]);
   assert(myTestBuffer[3] == TEST_FILE_CONTENTS[3]);
   assert(numBytesRead == 2);
   totalBytesPreviouslyRead += numBytesRead;
   // This read should not have affected the internal buffer.
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   // Should not be at eof
   assert(myFileTypePtr->eof() == false);
   assert(ifeof() == false);
  
   // Read the rest of the file.
   // Determine expected results for reading the rest of the file by
   // taking the substring starting after what had been previously read.
   numBytesRead = readFromFile(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   // Read the rest of the file, so the number of bytes read is
   // what was left in the file.
   assert(numBytesRead == (TEST_FILE_SIZE - totalBytesPreviouslyRead));
   assert(numBytesRead != MAX_TEST_BUFFER_SIZE);
   for(int i = 0; i < numBytesRead; i++)
   {
      assert(myTestBuffer[i] ==
	     TEST_FILE_CONTENTS[totalBytesPreviouslyRead+i]);
   }
   totalBytesPreviouslyRead += numBytesRead;
   assert(myFileTypePtr->eof() != 0);
   assert(ifeof() != 0);

   // Try to read one more time, making sure it doesn't read anything.
    numBytesRead = readFromFile(myTestBuffer, MAX_TEST_BUFFER_SIZE);
    assert(numBytesRead == 0);
   // Should be at eof
   assert(myFileTypePtr->eof() != 0);
   assert(ifeof() != 0);

   ifclose();

   std::cout << "  Passed test_readFromFile" << std::endl;
}




void IFILE_Test::test_readTilChar(const char* extension)
{
   // First open the test file.
   openFile(extension);

   // Verify the file successfully opened.
   assert(myFileTypePtr != NULL);
   assert(isOpen());
   assert(myFileTypePtr->isOpen());

   // Track position of ending char found.
   int pos = 0;

   // Test readTilChar.
   std::string output = "";
   std::string endChars = "a5d";
   pos = readTilChar(endChars, output);
   assert(pos == 0);  // read til a
   assert(output == "ABCD");
   output.clear();
   pos = readTilChar(endChars, output);
   assert(pos == 2);  // read til d
   assert(output == "bc");
   pos = readTilChar(endChars, output);
   assert(pos == 1);  // read til 5
   assert(output == "bc1234\nEFGefg");
   output.clear();
   pos = readTilChar(endChars, output);
   assert(pos == -1);  // read til 5
   assert(output == "67\nhijklHIJKL8910");

   ifrewind();
   // Test readTilChar.
   pos = readTilChar(endChars);
   assert(pos == 0);  // read til a
   pos = readTilChar(endChars);
   assert(pos == 2);  // read til d
   pos = readTilChar(endChars);
   assert(pos == 1);  // read til 5
   pos = readTilChar(endChars);
   assert(pos == -1);  // read til 5

   ifclose();

   std::cout << "  Passed test_readTilChar" << std::endl;
}


void IFILE_Test::test_ifeof_ifrewind(const char* extension)
{
   // First open the test file.
   openFile(extension);
   
   // Verify the file successfully opened.
   assert(myFileTypePtr != NULL);
   assert(isOpen());
   assert(myFileTypePtr->isOpen());
   
   // Not at eof - verify that it reports not eof.
   assert(ifeof() == false);
   
   // Track the total number of the bytes that have been read from the file
   // at any given point.
   int totalBytesPreviouslyRead = 0;
   int numBytesRead = 0;

   //////////////////////////////////////////////////////////////
   // Test doing reads from file without IFILE internal buffering.
   disableBuffering();

   // Verify position in file.
   assert(iftell() == 0);

   // Read a character from the file.
   numBytesRead = readFromFile(myTestBuffer, 1);
   assert(numBytesRead == 1);
   assert(myTestBuffer[0] == TEST_FILE_CONTENTS[totalBytesPreviouslyRead]);
   // Now that we have tested based on the previous total bytes read, 
   // increment the count.
   totalBytesPreviouslyRead += numBytesRead;
   // Not at eof
   assert(ifeof() == false);

   // Perform char read.
   char readChar = ifgetc();
   assert(readChar == TEST_FILE_CONTENTS[totalBytesPreviouslyRead]);
   // Now that we have tested based on the previous total bytes read, 
   // increment the count.
   ++totalBytesPreviouslyRead;
   // Not at eof
   assert(ifeof() == false);
   assert(iftell() == totalBytesPreviouslyRead);

   // Now read the rest.
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == (TEST_FILE_SIZE - totalBytesPreviouslyRead));
   // Hit the end of the file before reading the entire requested size.
   assert(numBytesRead != MAX_TEST_BUFFER_SIZE);
   // Now that we have tested based on the previous total bytes read, 
   // increment the count.
   totalBytesPreviouslyRead += numBytesRead;

   assert(myFileTypePtr->eof() != 0);
   assert(ifeof() != 0);
   
   numBytesRead = readFromFile(myTestBuffer, 1);
   assert(numBytesRead == 0);
   // Now it registers eof
   assert(ifeof() != 0);

   // bgzf files use a specialized return value for iftell that
   // is not just straight file offset.
   if((strcmp(extension, "bam") == 0) || (strcmp(extension, "glf") == 0))
   {
       assert(iftell() == (BGZF_TEST_FILE_SIZE << 16));
   }
   else
   {
      assert(iftell() == TEST_FILE_SIZE);
   }

   ///////////////////////////////////
   // Test doing IFILE buffered reads.
   // rewind the file and verify that it no longer registers eof.
   ifrewind();
   totalBytesPreviouslyRead = 0;
   // No longer at eof
   assert(ifeof() == false);
   // Verify position in file.
   assert(iftell() == 0);
   
   // Buffer reads - may have been disabled for iftell to work for bgzf.
   bufferReads();

   // Read a character from the file.
   numBytesRead = readFromFile(myTestBuffer, 1);
   assert(numBytesRead == 1);
   assert(myTestBuffer[0] == TEST_FILE_CONTENTS[totalBytesPreviouslyRead]);
   // Now that we have tested based on the previous total bytes read, 
   // increment the count.
   totalBytesPreviouslyRead += numBytesRead;
   // Not at eof
   assert(ifeof() == false);

   // Perform char read.
   readChar = ifgetc();
   assert(readChar == TEST_FILE_CONTENTS[totalBytesPreviouslyRead]);
   // Now that we have tested based on the previous total bytes read, 
   // increment the count.
   ++totalBytesPreviouslyRead;
   // Not at eof
   assert(ifeof() == false);
   
   // bgzf files use a specialized return value for iftell that
   // is not just straight file offset.
   if((strcmp(extension, "bam") == 0) || (strcmp(extension, "glf") == 0))
   {
       bool caught = false;
       try
       {
           assert(iftell() == totalBytesPreviouslyRead);
       }
       catch (std::exception& e)
       {
           caught = true;
           assert(strcmp(e.what(), "IFILE: CANNOT use buffered reads and tell for BGZF files") == 0);
       }
       assert(caught);
   }
   else
   {
      assert(iftell() == totalBytesPreviouslyRead);
   }

   // Now read the rest.
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == (TEST_FILE_SIZE - totalBytesPreviouslyRead));
   // Now that we have tested based on the previous total bytes read, 
   // increment the count.
   totalBytesPreviouslyRead += numBytesRead;
   // Registers eof.
   assert(ifeof() != 0);

   // Read past eof.
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == 0);
   // Eof.
   assert(ifeof() != 0);

   // bgzf files use a specialized return value for iftell that
   // is not just straight file offset.
   if((strcmp(extension, "bam") == 0) || (strcmp(extension, "glf") == 0))
   {
       bool caught = false;
       try
       {
           assert(iftell() == (BGZF_TEST_FILE_SIZE << 16));
       }
       catch (std::exception& e)
       {
           caught = true;
           assert(strcmp(e.what(), "IFILE: CANNOT use buffered reads and tell for BGZF files") == 0);
       }
       assert(caught);
       disableBuffering();
       assert(iftell() == (BGZF_TEST_FILE_SIZE << 16));
   }
   else
   {
      assert(iftell() == TEST_FILE_SIZE);
   }

   // Verify that after rewind, eof is no longer registered.
   ifrewind();
  // reset since we are back to the beginning of the file.
   totalBytesPreviouslyRead = 0;
   // No longer at eof
   assert(ifeof() == false);
   // Verify position in file.
   assert(iftell() == 0);

   // Verify properly works even if already at the beginning.   
   ifrewind();
  // reset since we are back to the beginning of the file.
   totalBytesPreviouslyRead = 0;
   // Not eof
   assert(ifeof() == false);
   // Verify position in file.
   assert(iftell() == 0);

   // Buffer reads - may have been disabled for iftell to work for bgzf.
   bufferReads();

   //////////////////////
   // Close the test file.
   ifclose();
   
   std::cout << "  Passed test_ifeof_ifrewind" << std::endl;
}


void IFILE_Test::test_ifread_ifgetc(const char* extension)
{
   // First open the test file.
   openFile(extension);

   // Verify the file successfully opened.
   assert(myFileTypePtr != NULL);
   assert(isOpen());
   assert(myFileTypePtr->isOpen());

   int numBytesRead = 0;
   int totalBytesPreviouslyRead = 0;

   ////////////////////////////////////
   // Test reading entire file at once.
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == TEST_FILE_SIZE);
   
   for(int i = 0; i < TEST_FILE_SIZE; i++)
   {
      assert(myTestBuffer[i] == TEST_FILE_CONTENTS[i]);
   }
   totalBytesPreviouslyRead += numBytesRead;
  
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == TEST_FILE_SIZE);
   assert(myBufferIndex == TEST_FILE_SIZE);
   
   assert(myFileTypePtr->eof() != 0);
   assert(ifeof() != 0);

   // Try reading at end of file twice.
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() != 0);

   // 2nd read attempt at eof.   
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() != 0);
   

   // RESET
   ifrewind();
   totalBytesPreviouslyRead = 0;

   //////////////////////////////////////////
   // Test reading entire file using getc.
   // Loop through reading the file.
   char readChar;
   for(int index = 0; index < TEST_FILE_SIZE; index++)
   {
      // Read a character.
      readChar = ifgetc();
      assert(readChar == TEST_FILE_CONTENTS[index]);
      // Should affect the IFILE buffer
      assert(myCurrentBufferSize == TEST_FILE_SIZE);
      assert(myBufferIndex == index+1);
   }
   
   // Now that we have read the file, try reading again at eof.
   readChar = ifgetc();
   assert(readChar == EOF);
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);

   // Try again at eof.
   // Now that we have read the file, try reading again at eof.
   readChar = ifgetc();
   assert(readChar == EOF);
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);

   // RESET
   ifrewind();
   totalBytesPreviouslyRead = 0;

   ////////////////////////////////////////////////
   // Test reading just the beginning of the file.
   numBytesRead = ifread(myTestBuffer, 4);
   assert(myTestBuffer[0] == TEST_FILE_CONTENTS[0]);
   assert(myTestBuffer[1] == TEST_FILE_CONTENTS[1]);
   assert(myTestBuffer[2] == TEST_FILE_CONTENTS[2]);
   assert(myTestBuffer[3] == TEST_FILE_CONTENTS[3]);
   assert(numBytesRead == 4);
   totalBytesPreviouslyRead += numBytesRead;
   // This read should have affected the internal buffer.
   assert(myCurrentBufferSize == TEST_FILE_SIZE);
   assert(myBufferIndex == 4);
   // Should not be at eof
   assert(ifeof() == false);

   // Test reading rest of file.
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == (TEST_FILE_SIZE - (int)totalBytesPreviouslyRead));
   // Verify contents of what read.   
    for(int i = 0; i < numBytesRead; i++)
    {
       assert(myTestBuffer[i] == 
	      TEST_FILE_CONTENTS[i + totalBytesPreviouslyRead]);
    }
    totalBytesPreviouslyRead += numBytesRead;

   // Try at end of file twice.
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() != 0);

   // 2nd read attempt at eof.   
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() != 0);
   
    // RESET
   ifrewind();
   totalBytesPreviouslyRead = 0;

   //////////////////////////////////////
   // Test reading just the beginning.
   numBytesRead = ifread(myTestBuffer, 4);
   assert(myTestBuffer[0] == TEST_FILE_CONTENTS[0]);
   assert(myTestBuffer[1] == TEST_FILE_CONTENTS[1]);
   assert(myTestBuffer[2] == TEST_FILE_CONTENTS[2]);
   assert(myTestBuffer[3] == TEST_FILE_CONTENTS[3]);
   assert(numBytesRead == 4);
   totalBytesPreviouslyRead += numBytesRead;
   // This read should have affected the internal buffer.
   assert(myCurrentBufferSize == TEST_FILE_SIZE);
   assert(myBufferIndex == 4);
   // Should not be at eof
   assert(ifeof() == false);

   // Test doing 2 getc.
   readChar = ifgetc();
   assert(readChar == TEST_FILE_CONTENTS[totalBytesPreviouslyRead]);
   int bufferSize = TEST_FILE_SIZE;
   assert(myCurrentBufferSize == bufferSize);
   assert(myBufferIndex == 5);
   totalBytesPreviouslyRead++;

   readChar = ifgetc();
   assert(readChar == TEST_FILE_CONTENTS[totalBytesPreviouslyRead]);
   assert(myCurrentBufferSize == bufferSize);
   assert(myBufferIndex == 6);
   totalBytesPreviouslyRead++;

   // Test reading rest of file.
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == (TEST_FILE_SIZE - (int)totalBytesPreviouslyRead));
   // Verify contents of what read.   
    for(int i = 0; i < numBytesRead; i++)
    {
       assert(myTestBuffer[i] == 
	      TEST_FILE_CONTENTS[i + totalBytesPreviouslyRead]);
    }
    totalBytesPreviouslyRead += numBytesRead;

   // Try at end of file twice.
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() != 0);

   // 2nd read attempt at eof.   
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() != 0);
   
    // RESET
   ifrewind();
   totalBytesPreviouslyRead = 0;
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);

   //////////////////////////////////   
   // Start with 2 getc.
   readChar = ifgetc();
   assert(readChar == TEST_FILE_CONTENTS[totalBytesPreviouslyRead]);
   bufferSize = TEST_FILE_SIZE;
   assert(myCurrentBufferSize == bufferSize);
   assert(myBufferIndex == 1);
   totalBytesPreviouslyRead++;

   readChar = ifgetc();
   assert(readChar == TEST_FILE_CONTENTS[totalBytesPreviouslyRead]);
   assert(myCurrentBufferSize == bufferSize);
   assert(myBufferIndex == 2);
   totalBytesPreviouslyRead++;

   // Test reading part of the rest of the file.
   numBytesRead = ifread(myTestBuffer, 4);
   assert(myTestBuffer[0] == TEST_FILE_CONTENTS[totalBytesPreviouslyRead]);
   assert(myTestBuffer[1] == TEST_FILE_CONTENTS[totalBytesPreviouslyRead + 1]);
   assert(myTestBuffer[2] == TEST_FILE_CONTENTS[totalBytesPreviouslyRead + 2]);
   assert(myTestBuffer[3] == TEST_FILE_CONTENTS[totalBytesPreviouslyRead + 3]);
   assert(numBytesRead == 4);
   totalBytesPreviouslyRead += numBytesRead;
   // This read should have affected the internal buffer.
   assert(myCurrentBufferSize == bufferSize);
   assert(myBufferIndex == totalBytesPreviouslyRead);
   // Should not be at eof
   assert(ifeof() == false);

   // Test reading 2 char with getc.
   readChar = ifgetc();
   assert(readChar == TEST_FILE_CONTENTS[totalBytesPreviouslyRead]);
   assert(myCurrentBufferSize == bufferSize);
   totalBytesPreviouslyRead++;
   assert(myBufferIndex == totalBytesPreviouslyRead);

   readChar = ifgetc();
   assert(readChar == TEST_FILE_CONTENTS[totalBytesPreviouslyRead]);
   assert(myCurrentBufferSize == bufferSize);
   totalBytesPreviouslyRead++;
   assert(myBufferIndex == totalBytesPreviouslyRead);

   // Test reading rest of file.
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == (TEST_FILE_SIZE - (int)totalBytesPreviouslyRead));
   // Verify contents of what read.   
    for(int i = 0; i < numBytesRead; i++)
    {
       assert(myTestBuffer[i] == 
	      TEST_FILE_CONTENTS[i + totalBytesPreviouslyRead]);
    }
    totalBytesPreviouslyRead += numBytesRead;
   assert(myBufferIndex == 0);
   assert(myCurrentBufferSize == 0);

   // Try at end of file twice.
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() != 0);

   // 2nd read attempt at eof.   
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() != 0);
   
    // RESET
   ifrewind();
   totalBytesPreviouslyRead = 0;
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);

   //////////////
   // Close the test file.
   ifclose();

   ////////////////////////////////////////////////////////////////////////
   // Repeat the test on a test file that is larger than the IFILE
   // buffer size.

   // First open the test file.
   openLargeFile(extension);

   // This file contains DEFAULT_BUFFER_SIZE of '0's followed by "12345"
   // The size of the file is DEFAULT_BUFFER_SIZE + 5.
   int largeTestFileSize = DEFAULT_BUFFER_SIZE + 5;
   char largeBuffer[largeTestFileSize + 5];

   // Verify the file successfully opened.
   assert(myFileTypePtr != NULL);
   assert(isOpen());
   assert(myFileTypePtr->isOpen());

   numBytesRead = 0;
   totalBytesPreviouslyRead = 0;

   ////////////////////////////////////
   // Test reading part of the file, then more then the buffer size,
   // then the rest of the file (test buffer handling when read
   // available and directly into the file, then read more).
   numBytesRead = ifread(largeBuffer, 2);
   assert(numBytesRead == 2);
   numBytesRead = ifread(largeBuffer + 2, DEFAULT_BUFFER_SIZE * 3);
   assert(numBytesRead == DEFAULT_BUFFER_SIZE + 3);
   // Should be at the end of the file.
   assert(myFileTypePtr->eof() != 0);
   assert(ifeof() != 0);
   numBytesRead = ifread(largeBuffer + DEFAULT_BUFFER_SIZE + 3, 2);
   assert(numBytesRead == 0);
   
   // Validate all the 0s
   for(unsigned int i = 0; i < DEFAULT_BUFFER_SIZE; i++)
   {
      assert(largeBuffer[i] == '0');
   }
   // Now validate the "12345"
   assert(largeBuffer[DEFAULT_BUFFER_SIZE] == '1');
   assert(largeBuffer[DEFAULT_BUFFER_SIZE+1] == '2');
   assert(largeBuffer[DEFAULT_BUFFER_SIZE+2] == '3');
   assert(largeBuffer[DEFAULT_BUFFER_SIZE+3] == '4');
   assert(largeBuffer[DEFAULT_BUFFER_SIZE+4] == '5');

   totalBytesPreviouslyRead += numBytesRead;
  
   // Should affect the IFILE buffer - 0 because read
   // is bigger than the buffer, so just read directly
   // into the largeBuffer.
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   
   assert(myFileTypePtr->eof() != 0);
   assert(ifeof() != 0);

   // Try reading at end of file twice.
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() != 0);

   // 2nd read attempt at eof.   
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() != 0);
   

   // RESET
   ifrewind();
   totalBytesPreviouslyRead = 0;

   ////////////////////////////////////
   // Test reading entire file at once.
   numBytesRead = ifread(largeBuffer, largeTestFileSize + 4);
   assert(numBytesRead == largeTestFileSize);
   
   // Validate all the 0s
   for(unsigned int i = 0; i < DEFAULT_BUFFER_SIZE; i++)
   {
      assert(largeBuffer[i] == '0');
   }
   // Now validate the "12345"
   assert(largeBuffer[DEFAULT_BUFFER_SIZE] == '1');
   assert(largeBuffer[DEFAULT_BUFFER_SIZE+1] == '2');
   assert(largeBuffer[DEFAULT_BUFFER_SIZE+2] == '3');
   assert(largeBuffer[DEFAULT_BUFFER_SIZE+3] == '4');
   assert(largeBuffer[DEFAULT_BUFFER_SIZE+4] == '5');

   totalBytesPreviouslyRead += numBytesRead;
  
   // Should affect the IFILE buffer - 0 because read
   // is bigger than the buffer, so just read directly
   // into the largeBuffer.
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   
   assert(myFileTypePtr->eof() != 0);
   assert(ifeof() != 0);

   // Try reading at end of file twice.
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() != 0);

   // 2nd read attempt at eof.   
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() != 0);
   

   // RESET
   ifrewind();
   totalBytesPreviouslyRead = 0;

   //////////////////////////////////////////
   // Test reading entire file using getc.
   // Loop through reading the file.
   // First loop through verifying the 0's
   for(int index = 0; index < (int)DEFAULT_BUFFER_SIZE; index++)
   {
      // Read a character.
      readChar = ifgetc();
      assert(readChar == '0');
      // Should affect the IFILE buffer
      assert(myCurrentBufferSize == (int)DEFAULT_BUFFER_SIZE);
      assert(myBufferIndex == index+1);
   }
   // Now read the 12345.
   readChar = ifgetc();
   assert(readChar == '1');
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 5);
   assert(myBufferIndex == 1);
   readChar = ifgetc();
   assert(readChar == '2');
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 5);
   assert(myBufferIndex == 2);
   readChar = ifgetc();
   assert(readChar == '3');
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 5);
   assert(myBufferIndex == 3);
   readChar = ifgetc();
   assert(readChar == '4');
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 5);
   assert(myBufferIndex == 4);
   readChar = ifgetc();
   assert(readChar == '5');
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 5);
   assert(myBufferIndex == 5);

   // Now that we have read the file, try reading again at eof.
   readChar = ifgetc();
   assert(readChar == EOF);
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);

   // Try again at eof.
   // Now that we have read the file, try reading again at eof.
   readChar = ifgetc();
   assert(readChar == EOF);
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);

   // RESET
   ifrewind();
   totalBytesPreviouslyRead = 0;

   ////////////////////////////////////////////////
   // Test reading just the beginning of the file.
   numBytesRead = ifread(largeBuffer, 4);
   assert(largeBuffer[0] == '0');
   assert(largeBuffer[1] == '0');
   assert(largeBuffer[2] == '0');
   assert(largeBuffer[3] == '0');
   assert(numBytesRead == 4);
   totalBytesPreviouslyRead += numBytesRead;
   // This read should have affected the internal buffer.
   assert(myCurrentBufferSize == (int)DEFAULT_BUFFER_SIZE);
   assert(myBufferIndex == 4);
   // Should not be at eof
   assert(ifeof() == false);

   // Test reading rest of file.
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == (largeTestFileSize - (int)totalBytesPreviouslyRead));
   // Verify contents of what read.   First check the 0's
   for(int i = 0; i < (numBytesRead-5); i++)
   {
      assert(largeBuffer[i] == '0');
   }
   // Check the 12345
   assert(largeBuffer[numBytesRead - 5] == '1');
   assert(largeBuffer[numBytesRead - 5 + 1] == '2');
   assert(largeBuffer[numBytesRead - 5 + 2] == '3');
   assert(largeBuffer[numBytesRead - 5 + 3] == '4');
   assert(largeBuffer[numBytesRead - 5 + 4] == '5');
   totalBytesPreviouslyRead += numBytesRead;
   
   // Try at end of file twice.
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == 0);
   // Trying to read at the end cleared the buffer..
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() != 0);

   // 2nd read attempt at eof.   
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() != 0);
   
    // RESET
   ifrewind();
   totalBytesPreviouslyRead = 0;

   //////////////////////////////////////
   // Test reading just the beginning.
   numBytesRead = ifread(largeBuffer, 2);
   assert(largeBuffer[0] == '0');
   assert(largeBuffer[1] == '0');
   assert(numBytesRead == 2);
   totalBytesPreviouslyRead += numBytesRead;
   // This read should have affected the internal buffer.
   assert(myCurrentBufferSize == (int)DEFAULT_BUFFER_SIZE);
   assert(myBufferIndex == 2);
   // Should not be at eof
   assert(ifeof() == false);

   // Test doing 2 getc.
   readChar = ifgetc();
   assert(readChar == '0');
   bufferSize = DEFAULT_BUFFER_SIZE;
   assert(myCurrentBufferSize == bufferSize);
   assert(myBufferIndex == 3);
   totalBytesPreviouslyRead++;

   readChar = ifgetc();
   assert(readChar == '0');
   assert(myCurrentBufferSize == bufferSize);
   assert(myBufferIndex == 4);
   totalBytesPreviouslyRead++;

   // Test reading rest of file.
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == (largeTestFileSize - (int)totalBytesPreviouslyRead));
   // Verify contents of what read.   
   // All except the last 5 should be '0'
    for(int i = 0; i < numBytesRead - 5; i++)
    {
       assert(largeBuffer[i] == '0');
    }
    assert(largeBuffer[numBytesRead - 5] == '1');
    assert(largeBuffer[numBytesRead - 4] == '2');
    assert(largeBuffer[numBytesRead - 3] == '3');
    assert(largeBuffer[numBytesRead - 2] == '4');
    assert(largeBuffer[numBytesRead - 1] == '5');

    totalBytesPreviouslyRead += numBytesRead;

   // Try at end of file twice.
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == 0);
   // Reading at the end clears the buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() != 0);

   // 2nd read attempt at eof.   
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == 0);
   // Reading at the end clears the buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() != 0);
   
    // RESET
   ifrewind();
   totalBytesPreviouslyRead = 0;
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);

   //////////////////////////////////   
   // Start with 2 getc.
   readChar = ifgetc();
   assert(readChar == '0');
   bufferSize = DEFAULT_BUFFER_SIZE;
   assert(myCurrentBufferSize == bufferSize);
   assert(myBufferIndex == 1);
   totalBytesPreviouslyRead++;

   readChar = ifgetc();
   assert(readChar == '0');
   assert(myCurrentBufferSize == bufferSize);
   assert(myBufferIndex == 2);
   totalBytesPreviouslyRead++;

   // Test reading part of the rest of the file.
   numBytesRead = ifread(myTestBuffer, 2);
   assert(myTestBuffer[0] == '0');
   assert(myTestBuffer[1] == '0');
   assert(numBytesRead == 2);
   totalBytesPreviouslyRead += numBytesRead;
   // This read should have affected the internal buffer.
   assert(myCurrentBufferSize == bufferSize);
   assert(myBufferIndex == totalBytesPreviouslyRead);
   // Should not be at eof
   assert(ifeof() == false);

   // Test reading 2 char with getc.
   readChar = ifgetc();
   assert(readChar == '0');
   assert(myCurrentBufferSize == bufferSize);
   totalBytesPreviouslyRead++;
   assert(myBufferIndex == totalBytesPreviouslyRead);

   readChar = ifgetc();
   assert(readChar == '0');
   assert(myCurrentBufferSize == bufferSize);
   totalBytesPreviouslyRead++;
   assert(myBufferIndex == totalBytesPreviouslyRead);

   // Test reading rest of file.
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == (largeTestFileSize - (int)totalBytesPreviouslyRead));
   // Verify contents of what read.   
   for(int i = 0; i < numBytesRead - 5; i++)
   {
      assert(largeBuffer[i] == '0');
   }
   // Verify the 12345
   assert(largeBuffer[numBytesRead - 5] == '1');
   assert(largeBuffer[numBytesRead - 5 + 1] == '2');
   assert(largeBuffer[numBytesRead - 5 + 2] == '3');
   assert(largeBuffer[numBytesRead - 5 + 3] == '4');
   assert(largeBuffer[numBytesRead - 5 + 4] == '5');
   totalBytesPreviouslyRead += numBytesRead;
   bufferSize = 5;
   assert(myBufferIndex == bufferSize);
   assert(myCurrentBufferSize == bufferSize);

   // Try at end of file twice.
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == 0);
   // Reading at the end clears the buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() != 0);

   // 2nd read attempt at eof.   
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == 0);
   // Reading at the end clears the buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() != 0);
   
    // RESET
   ifrewind();
   totalBytesPreviouslyRead = 0;
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);

   ifclose();

   std::cout << "  Passed test_ifread_ifgetc" << std::endl;
}


// Test closing a file.
void IFILE_Test::test_ifclose(const char* extension)
{
   // First open the test file.
   openFile(extension);

   // Verify the file successfully opened.
   assert(myFileTypePtr != NULL);
   assert(isOpen());
   assert(myFileTypePtr->isOpen());

   ifclose();

   assert(myFileTypePtr == NULL);
   assert(isOpen() == false);

   std::cout << "  Passed test_ifclose" << std::endl;
}


void IFILE_Test::test_ifseek(const char* extension)
{
    disableBuffering();
   // First open the test file.
   openFile(extension);

   // Read a character from the file.
   int numBytesRead = readFromFile(myTestBuffer, 1);
   assert(numBytesRead == 1);
   assert(myTestBuffer[0] == TEST_FILE_CONTENTS[0]);

   // Get the current position.
   long int currentPos = iftell();
   
   // Read the next character from the file.
   numBytesRead = readFromFile(myTestBuffer, 1);
   assert(numBytesRead == 1);
   assert(myTestBuffer[0] == TEST_FILE_CONTENTS[1]);

   // Seek to just before the character that was just read and read again
   // Should be the same character.
   assert(ifseek(currentPos, SEEK_SET) == true);
   numBytesRead = readFromFile(myTestBuffer, 1);
   assert(numBytesRead == 1);
   assert(myTestBuffer[0] == TEST_FILE_CONTENTS[1]);
   
   ifclose();

   assert(myFileTypePtr == NULL);
   assert(isOpen() == false);

   // Buffer reads - may have been disabled for iftell to work for bgzf.
   bufferReads();

   std::cout << "  Passed test_ifseek" << std::endl;
}

void IFILE_Test::test_noExistRead(const char* extension)
{
   openNoExistFile(extension);

}


// Open a file for testing.
void IFILE_Test::openFile(const char* extension)
{
   std::string filename = "data/InputFileTest.";
   filename += extension;
   assert(InputFile::openFile(filename.c_str(), "rb", InputFile::DEFAULT) == true);
}

// Open a file for testing.
void IFILE_Test::openLargeFile(const char* extension)
{
   std::string filename = "data/InputFileTestLarge.";
   filename += extension;
   assert(InputFile::openFile(filename.data(), "rb", InputFile::DEFAULT) == true);
}


void IFILE_Test::openNoExistFile(const char* extension)
{
   std::string filename = "data/noExist.";
   filename += extension;
   assert(InputFile::openFile(filename.data(), "rb", InputFile::DEFAULT) == false);
}


void testWrite()
{
    std::string filenameNoExt = "results/InputFileTest.";
    std::string filename = filenameNoExt + "glf";
    
    IFILE filePtr = ifopen(filename.c_str(), "wt");
    assert(filePtr != NULL);
    
    assert(ifwrite(filePtr, 
                   IFILE_Test::TEST_FILE_CONTENTS.c_str(), 
                   IFILE_Test::TEST_FILE_CONTENTS.length()) 
           == IFILE_Test::TEST_FILE_CONTENTS.length());
    
    assert(ifclose(filePtr) == 0);

    filename = "results/uncompressedFile.glf";
    
    filePtr = ifopen(filename.c_str(), "wt", InputFile::UNCOMPRESSED);
    assert(filePtr != NULL);
    
    assert(ifwrite(filePtr,
                   IFILE_Test::TEST_FILE_CONTENTS.c_str(), 
                   IFILE_Test::TEST_FILE_CONTENTS.length()) 
           == IFILE_Test::TEST_FILE_CONTENTS.length());
    
    assert(ifclose(filePtr) == 0);

    filename = "results/bgzfFile.glf";
    
    filePtr = ifopen(filename.c_str(), "wt", InputFile::BGZF);
    assert(filePtr != NULL);
    
    assert(ifwrite(filePtr,
                   IFILE_Test::TEST_FILE_CONTENTS.c_str(), 
                   IFILE_Test::TEST_FILE_CONTENTS.length()) 
           == IFILE_Test::TEST_FILE_CONTENTS.length());
    
    assert(ifclose(filePtr) == 0);

    filename = "results/gzipFile.glf";
    
    filePtr = ifopen(filename.c_str(), "wt", InputFile::GZIP);
    assert(filePtr != NULL);
    
    assert(ifwrite(filePtr,
                   IFILE_Test::TEST_FILE_CONTENTS.c_str(), 
                   IFILE_Test::TEST_FILE_CONTENTS.length()) 
           ==IFILE_Test:: TEST_FILE_CONTENTS.length());
    
    assert(ifclose(filePtr) == 0);

    filename = "results/defaultFile.glf";
    
    filePtr = ifopen(filename.c_str(), "wt");
    assert(filePtr != NULL);
    
    assert(ifwrite(filePtr,
                   IFILE_Test::TEST_FILE_CONTENTS.c_str(), 
                   IFILE_Test::TEST_FILE_CONTENTS.length()) 
           == IFILE_Test::TEST_FILE_CONTENTS.length());
    
    assert(ifclose(filePtr) == 0);

    filename = "results/defaultFile.gz";
    
    filePtr = ifopen(filename.c_str(), "wt");
    assert(filePtr != NULL);
    
    assert(ifwrite(filePtr,
                   IFILE_Test::TEST_FILE_CONTENTS.c_str(), 
                   IFILE_Test::TEST_FILE_CONTENTS.length()) 
           == IFILE_Test::TEST_FILE_CONTENTS.length());
    
    assert(ifclose(filePtr) == 0);


    filename = "results/textFile.gz";
    
    unsigned int myuint = 99;
    int myint = -99;
    char mychar = 'z';

    filePtr = ifopen(filename.c_str(), "wt");
    (*filePtr) << "Hello\n";
    (*filePtr) << "Hello." << 3 << ' ' << -2 << "How are you";
    (*filePtr) << "?" << "\n";
    std::string mytext = "Bye\n";
    (*filePtr) << mytext;
    (*filePtr) << 3.125 << mychar;
    (*filePtr) << myuint;
    (*filePtr) << mychar;
    (*filePtr) << myint;
    String myString = "Good Bye!\n";
    (*filePtr) << myString;
    assert(ifclose(filePtr) == 0);

    filename = "results/textFile1.gz";
    InputFile& fileRef = *(ifopen(filename.c_str(), "wt"));
    fileRef << "Hello\n";
    fileRef << "Hello." << 3 << ' ' << -2 << "How are you";
    fileRef << "?" << "\n";
    fileRef << mytext;
    fileRef << 3.125 << mychar;
    fileRef << myuint;
    fileRef << mychar;
    fileRef << myint;
    fileRef << myString;
    InputFile* fileRefPtr = &fileRef;
    assert(ifclose(fileRefPtr) == 0);
    assert(fileRefPtr == NULL);

    // TODO - automatically verify that the files were written in the
    // correct format - rather than hand checking.
}



void testAdditional(const char* extension)
{
    std::string fileName = "data/InputFileTest2.";
    fileName += extension;
    IFILE testFile = ifopen(fileName.c_str(), "r");
    assert(testFile != NULL);

    std::string buffer = "989";
    std::string stopChars = "C5F2";

    // Test readTilChar that stores the string.
    assert(testFile->readTilChar(stopChars, buffer) == 0);
    assert(buffer == "989AB");
    buffer.clear();
    assert(testFile->readTilChar(stopChars, buffer) == 2);
    assert(buffer == "DE");
    assert(testFile->readTilChar(stopChars, buffer) == 3);
    assert(buffer == "DEG\tabcdefg\n1");

    // Test readTilChar that discards the string.
    assert(testFile->readTilChar(stopChars) == 1);
    buffer.clear();
    buffer = "t";
    assert(testFile->readTilTab(buffer) == 1);
    assert(buffer == "t6");
    assert(testFile->readTilTab(buffer) == 0);
    assert(buffer == "t6hijklm");
    assert(testFile->readTilTab(buffer) == 0);
    assert(buffer == "t6hijklm1");
    assert(testFile->readTilTab(buffer) == 1);
    assert(buffer == "t6hijklm1NOP");
    assert(testFile->readLine(buffer) == 0);
    assert(buffer == "t6hijklm1NOPQRST\tUVW");
    assert(testFile->readTilTab(buffer) == 0);
    assert(buffer == "t6hijklm1NOPQRST\tUVW");
    buffer.clear();
    assert(testFile->discardLine() == 0);
    assert(testFile->readLine(buffer) == -1);
    assert(buffer == "@#$");
    assert(testFile->discardLine() == -1);
    assert(testFile->readTilTab(buffer) == -1);
    assert(testFile->readTilChar(stopChars, buffer) == -1);
    assert(testFile->readTilChar(stopChars) == -1);
    assert(buffer == "@#$");

    ifclose(testFile);

}
