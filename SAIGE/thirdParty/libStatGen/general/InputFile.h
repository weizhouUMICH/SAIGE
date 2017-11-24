/*
 *  Copyright (C) 2010-2012  Regents of the University of Michigan
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
/*! \file */ 
#ifndef __INPUTFILE_H__
#define __INPUTFILE_H__

#include <stdio.h>
#include <iostream>
#include <cstring>
#include <stdint.h>

#include "FileType.h"

/// Class for easily reading/writing files without having to worry about
/// file type (uncompressed, gzip, bgzf) when reading.
/// It hides the low level file operations/structure from the user, allowing
/// them to generically open and operate on a file using the same
/// interface without knowing the file format (standard uncompressed,
/// gzip, or bgzf).  For writing, the user must specify the file type.
/// There is a typedef IFILE which is InputFile* and setup to mimic FILE
/// including global methods that take IFILE as a parameter.
class InputFile
{
    bool    myAttemptRecovery;  // use recovery techniques if possible
public:

    /// Compression to use when writing a file & decompression used when
    /// reading a file from stdin.  Any other read checks the file to determine
    ///  how to uncompress it.
    enum ifileCompression {
        DEFAULT,  ///< Check the extension, if it is ".gz", treat as gzip, otherwise treat it as UNCOMPRESSED.
        UNCOMPRESSED,  ///< uncompressed file.
        GZIP,  ///< gzip file.
        BGZF ///< bgzf file.
    };

    /// Default constructor
    InputFile()
    {
        myAttemptRecovery = false;
        myFileTypePtr = NULL;
        myBufferIndex = 0;
        myCurrentBufferSize = 0;
        // Default to buffer.
        myAllocatedBufferSize = DEFAULT_BUFFER_SIZE;
        myFileBuffer = new char[myAllocatedBufferSize];
        myFileName.clear();
    }

    /// Destructor
    ~InputFile();

    /// Constructor for opening a file.
    /// \param filename file to open
    /// \param mode same format as fopen: "r" for read & "w" for write.
    /// \param compressionMode set the type of file to open for writing or
    /// for reading from stdin (when reading files, the compression type is 
    /// determined by reading the file).
    InputFile(const char * filename, const char * mode,
              InputFile::ifileCompression compressionMode = InputFile::DEFAULT);

    /// Set the buffer size for reading from files so that bufferSize bytes
    /// are read at a time and stored until accessed by another read call.
    /// This improves performance over reading the file small bits at a time.
    /// Buffering reads disables the tell call for bgzf files.
    /// Any previous values in the buffer will be deleted.
    /// \param bufferSize number of bytes to read/buffer at a time,
    /// turn off read buffering by setting bufferSize = 1;
    inline void bufferReads(unsigned int bufferSize = DEFAULT_BUFFER_SIZE)
    {
        // If the buffer size is the same, do nothing.
        if(bufferSize == myAllocatedBufferSize)
        {
            return;
        }
        // Delete the previous buffer.
        if(myFileBuffer != NULL)
        {
            delete[] myFileBuffer;
        }
        myBufferIndex = 0;
        myCurrentBufferSize = 0;
        // The buffer size must be at least 1 so one character can be
        // read and ifgetc can just assume reading into the buffer.
        if(bufferSize < 1)
        {
            bufferSize = 1;
        }
        myFileBuffer = new char[bufferSize];
        myAllocatedBufferSize = bufferSize;

        if(myFileTypePtr != NULL)
        {
            if(bufferSize == 1)
            {
                myFileTypePtr->setBuffered(false);
            }
            else
            {
                myFileTypePtr->setBuffered(true);
            }
        }
    }


    /// Disable read buffering.
    inline void disableBuffering()
    {
        bufferReads(1);
        if(myFileTypePtr != NULL)
        {
            myFileTypePtr->setBuffered(false);
        }
    }

    
    /// Close the file.
    /// \return status of the close (0 is success).
    inline int ifclose()
    {
        if (myFileTypePtr == NULL)
        {
            return EOF;
        }
        int result = myFileTypePtr->close();
        delete myFileTypePtr;
        myFileTypePtr = NULL;
        myFileName.clear();
        return result;
    }

    /// Read size bytes from the file into the buffer.
    /// \param buffer pointer to memory at least size bytes big to write the
    /// data into.
    /// \param size number of bytes to be read
    /// \return number of bytes read, if it is not equal to size,
    /// there was either an error or the end of the file was reached, use
    /// ifeof to determine which case it was.
    inline int ifread(void * buffer, unsigned int size)
    {
        // There are 2 cases:
        //  1) There are already size available bytes in buffer.
        //  2) There are not size bytes in buffer.

        // Determine the number of available bytes in the buffer.
        unsigned int availableBytes = myCurrentBufferSize - myBufferIndex;
        int returnSize = 0;

        // Case 1: There are already size available bytes in buffer.
        if (size <= availableBytes)
        {
            //   Just copy from the buffer, increment the index and return.
            memcpy(buffer, myFileBuffer+myBufferIndex, size);
            // Increment the buffer index.
            myBufferIndex += size;
            returnSize = size;
        }
        // Case 2: There are not size bytes in buffer.
        else
        {
            // Check to see if there are some bytes in the buffer.
            if (availableBytes > 0)
            {
                // Size > availableBytes > 0
                // Copy the available bytes into the buffer.
                memcpy(buffer, myFileBuffer+myBufferIndex, availableBytes);
            }
            // So far availableBytes have been copied into the read buffer.
            returnSize = availableBytes;
            // Increment myBufferIndex  by what was read.
            myBufferIndex += availableBytes;

            unsigned int remainingSize = size - availableBytes;

            // Check if the remaining size is more or less than the
            // max buffer size.
            if(remainingSize < myAllocatedBufferSize)
            {
                // the remaining size is not the full buffer, but read
                //  a full buffer worth of data anyway.
                myCurrentBufferSize =
                    readFromFile(myFileBuffer, myAllocatedBufferSize);

                // Check for an error.
                if(myCurrentBufferSize <= 0)
                {
                    // No more data was successfully read, so check to see
                    // if any data was copied to the return buffer at all.
                    if( returnSize == 0)
                    {
                        // No data has been copied at all into the
                        // return read buffer, so just return the value
                        // returned from readFromFile.
                        returnSize = myCurrentBufferSize;
                        // Otherwise, returnSize is already set to the
                        // available bytes that was already copied (so no
                        // else statement is needed).
                    }
                    // Set myBufferIndex & myCurrentBufferSize to 0.
                    myCurrentBufferSize = 0;
                    myBufferIndex = 0;
                }
                else
                {
                    // Successfully read more data.
                    // Check to see how much was copied.
                    int copySize = remainingSize;
                    if(copySize > myCurrentBufferSize)
                    {
                        // Not the entire requested amount was read
                        // (either from EOF or there was a partial read due to
                        // an error), so set the copySize to what was read.
                        copySize = myCurrentBufferSize;
                    }

                    // Now copy the rest of the bytes into the buffer.
                    memcpy((char*)buffer+availableBytes, 
                           myFileBuffer, copySize);

                    // set the buffer index to the location after what we are
                    // returning as read.
                    myBufferIndex = copySize;
                
                    returnSize += copySize;
                }
            }
            else
            {
                // More remaining to be read than the max buffer size, so just
                // read directly into the output buffer.
                int readSize = readFromFile((char*)buffer + availableBytes,
                                            remainingSize);

                // Already used the buffer, so "clear" it.
                myCurrentBufferSize = 0;
                myBufferIndex = 0;
                if(readSize <= 0)
                {
                    // No more data was successfully read, so check to see
                    // if any data was copied to the return buffer at all.
                    if(returnSize == 0)
                    {
                        // No data has been copied at all into the
                        // return read buffer, so just return the value
                        // returned from readFromFile.
                        returnSize = readSize;
                        // Otherwise, returnSize is already set to the
                        // available bytes that was already copied (so no
                        // else statement is needed).
                    }
                }
                else
                {
                    // More data was read, so increment the return count.
                    returnSize += readSize;
                }
            }
        }
        return(returnSize);
    }

    /// Read until the specified characters, returning which character was
    /// found causing the stop, -1 returned for EOF, storing the other read
    /// characters into the specified string.
    /// Note: If stopChars is just '\n', readLine is faster and if
    /// stopChars is just '\n' and '\t', readTilTab is faster.
    /// \param stopChars characters to stop reading when they are hit.
    /// \param stringRef reference to a string that the read characters should
    /// be apppended to (does not include the stopchar).
    /// \return index of the character in stopChars that caused it to stop
    /// reading or -1 for EOF.
    int readTilChar(const std::string& stopChars, std::string& stringRef);

    /// Read until the specified characters, returning which character was
    /// found causing the stop, -1 returned for EOF, dropping all read chars.
    /// Note: If stopChars is just '\n', discardLine is faster.
    /// \param stopChars characters to stop reading when they are hit.
    /// \return index of the character in stopChars that caused it to stop
    /// reading or -1 for EOF.
    int readTilChar(const std::string& stopChars);

    /// Read until the end of the line, discarding the characters,
    /// returning -1 returned for EOF and returning 0 if the end of the line
    /// was found.
    /// \return 0 if the end of the line was found before EOF or -1 for EOF.
    int discardLine();

    /// Read, appending the characters into the specified string until new
    /// line or EOF is found, returning -1 if EOF is found first and 0 if new
    /// line is found first.  The new line and EOF are not written into the 
    /// specified string.
    /// \param line reference to a string that the read characters should
    /// be apppended to (does not include the new line or eof).
    /// \return 0 if new line and -1 for EOF.
    int readLine(std::string& line);

    /// Read, appending the characters into the specified string until tab, new
    /// line, or EOF is found, returning -1 if EOF is found first, 0 if new
    /// line is found first, or 1 if a tab is found first.  The tab, new line,
    /// and EOF are not written into the specified string.
    /// \param field reference to a string that the read characters should
    /// be apppended to (does not include the tab, new line, or eof).
    /// \return 1 if tab is found, 0 if new line, and -1 for EOF.
    int readTilTab(std::string& field);

    /// Get a character from the file.  Read a character from the internal
    /// buffer, or if the end of the buffer has been reached, read from the
    /// file into the buffer and return index 0.
    /// \return character that was read or EOF.
    inline int ifgetc()
    {
        if (myBufferIndex >= myCurrentBufferSize)
        {
            // at the last index, read a new buffer.
            myCurrentBufferSize = readFromFile(myFileBuffer, myAllocatedBufferSize);
            myBufferIndex = 0;
            // If the buffer index is still greater than or equal to the
            // myCurrentBufferSize, then we failed to read the file - return EOF.
            // NB: This only needs to be checked when myCurrentBufferSize
            // is changed.  Simplify check - readFromFile returns zero on EOF
            if (myCurrentBufferSize == 0)
            {
                return(EOF);
            }
        }
        return(myFileBuffer[myBufferIndex++]);
    }

    /// Get a line from the file.
    /// \param buffer the buffer into which data is to be placed
    /// \param max the maximum size of the buffer, in bytes
    /// \return true if the last character read was an EOF
    inline bool ifgetline(void *voidBuffer, size_t max)
    {
        int ch;
        char *buffer = (char *) voidBuffer;

        while( (ch=ifgetc()) != '\n' && ch != EOF) {
            *buffer++ = ch;
            if((--max)<2)
            {
                // truncate the line, so drop remainder
                while( (ch=ifgetc()) && ch != '\n' && ch != EOF)
                {
                }
                break;
            }
        }
        *buffer++ = '\0';
        return ch==EOF;
    }

    /// Reset to the beginning of the file.
    inline void ifrewind()
    {
        // Just set the myBufferIndex and the myCurrentBufferSize to 0 to simulate
        // clearing the buffer and call rewind to move to the beginning of the
        // file.
        if (myFileTypePtr == NULL)
        {
            // No pointer, so nothing to rewind.
            return;
        }
        myCurrentBufferSize = 0;
        myBufferIndex = 0;
        myFileTypePtr->rewind();
    }


    /// Check to see if we have reached the EOF.
    /// \return 0 if not EOF, any other value means EOF.
    inline int ifeof() const
    {
        // Not EOF if we are not at the end of the buffer.
        if (myBufferIndex < myCurrentBufferSize)
        {
            // There are still available bytes in the buffer, so NOT EOF.
            return false;
        }
        else
        {
            if (myFileTypePtr == NULL)
            {
                // No myFileTypePtr, so not eof (return 0).
                return 0;
            }
            // exhausted our buffer, so check the file for eof.
            return myFileTypePtr->eof();
        }
    }

    /// Write the specified buffer into the file.
    /// \param buffer buffer containing size bytes to write to the file.
    /// \param size number of bytes to write
    /// \return number of bytes written
    /// We do not buffer the write call, so just leave this as normal.
    inline unsigned int ifwrite(const void * buffer, unsigned int size)
    {
        if (myFileTypePtr == NULL)
        {
            // No myFileTypePtr, so return 0 - nothing written.
            return 0;
        }
        return myFileTypePtr->write(buffer, size);
    }

    /// Returns whether or not the file was successfully opened.
    /// \return true if the file is open, false if not.
    inline bool isOpen() const
    {
        // It is open if the myFileTypePtr is set and says it is open.
        if ((myFileTypePtr != NULL) && myFileTypePtr->isOpen())
        {
            return true;
        }
        // File was not successfully opened.
        return false;
    }

    /// Get current position in the file.
    /// \return current position in the file, -1 indicates an error.
    inline int64_t iftell()
    {
        if (myFileTypePtr == NULL)
        {
            // No myFileTypePtr, so return false - could not seek.
            return -1;
        }
        int64_t pos = myFileTypePtr->tell();
        pos -= (myCurrentBufferSize - myBufferIndex);
        return(pos);
    }


    /// Seek to the specified offset from the origin.
    /// \param offset offset into the file to move to (must be from a tell call)
    /// \param origin can be any of the following:
    /// Note: not all are valid for all filetypes.
    ///   SEEK_SET - Beginning of file
    ///   SEEK_CUR - Current position of the file pointer
    ///   SEEK_END - End of file
    /// \return true on successful seek and false on a failed seek.
    inline bool ifseek(int64_t offset, int origin)
    {
        if (myFileTypePtr == NULL)
        {
            // No myFileTypePtr, so return false - could not seek.
            return false;
        }
        // TODO - may be able to seek within the buffer if applicable.
        // Reset buffering since a seek is being done.
        myBufferIndex = 0;
        myCurrentBufferSize = 0;
        return myFileTypePtr->seek(offset, origin);
    }

    /// Get the filename that is currently opened.
    /// \return filename associated with this class
    const char* getFileName() const
    {
        return(myFileName.c_str());
    }

    /// Enable (default) or disable recovery.
    /// 
    /// When true, we can attach a myFileTypePtr
    /// that implements a recovery capable decompressor.
    /// This requires that the caller be able to catch
    /// the exception XXX "blah blah blah".
    ///
    void setAttemptRecovery(bool flag = false)
    {
        myAttemptRecovery = flag;
    }

    bool attemptRecoverySync(bool (*checkSignature)(void *data) , int length)
    {
        if(myFileTypePtr==NULL) return false; 
        return myFileTypePtr->attemptRecoverySync(checkSignature, length);
    }

    // Open a file. Called by the constructor.
    // Returns true if the file was successfully opened, false otherwise.
    bool openFile(const char * filename, const char * mode,
                  InputFile::ifileCompression compressionMode);

protected:
    // Read into a buffer from the file.  Since the buffer is passed in and
    // this would bypass the myFileBuffer used by this class, this method must
    // be protected.
    inline int readFromFile(void * buffer, unsigned int size)
    {
        // If no myFileTypePtr, return 0 - nothing read.
        if (myFileTypePtr == NULL)
        {
            return 0;
        }
        return myFileTypePtr->read(buffer, size);
    }

#ifdef __ZLIB_AVAILABLE__
    // Only necessary with zlib to determine what file type on a new
    // file.  Without zlib, there are only uncompressed files, so a special
    // method is not needed to determine the type of file to open.
    // Open a file.  This method will open a file with the specified name and
    // mode with the fileTypePtr associated with the specified compressionMode.
    void openFileUsingMode(const char* filename, const char* mode,
                           InputFile::ifileCompression compressionMode);
#endif

    // The size of the buffer used by this class.
    static const unsigned int DEFAULT_BUFFER_SIZE = 65536;

    // Pointer to a class that interfaces with different file types.
    FileType* myFileTypePtr;

    unsigned int myAllocatedBufferSize;

    // Buffer used to do large reads rather than 1 by 1 character reads
    // from the file.  The class is then managed to iterate through the buffer.
    char* myFileBuffer;

    // Current index into the buffer.  Used to track where we are in reading the
    // file from the buffer.
    int myBufferIndex;

    // Current number of entries in the buffer.  Used to ensure that
    // if a read did not fill the buffer, we stop before hitting the
    // end of what was read.
    int myCurrentBufferSize;

    std::string myFileName;
};


/// Define IFILE as a pointer to an InputFile object.
typedef InputFile* IFILE;


/// Open a file with the specified name and mode, using a filename of "-" to 
/// indicate stdin/stdout.
/// \param filename file to open ("-" meands stdin/stdout)
/// \param mode same format as fopen: "r" for read & "w" for write.
/// \param compressionMode set the type of file to open for writing or
/// for reading from stdin (when reading files not from stdin, the compression
/// type is determined by reading the file).
/// \return IFILE - pointer to the InputFile object that has been opened.
inline IFILE ifopen(const char * filename, const char * mode,
                    InputFile::ifileCompression compressionMode = InputFile::DEFAULT)
{
    IFILE file = new InputFile(filename, mode, compressionMode);
    if (!file->isOpen())
    {

        // Not open, so delete the file, and return null.
        delete file;
        file = NULL;
    }
    return file;
}


/// Close the file.
/// \param file file to be closed - IFILE is a pointer to an InputFile object
/// \return status of the close (0 is success or if NULL is passed in).
inline int ifclose(IFILE &file)
{
    if(file == NULL)
    {
        // NULL Pointer passed in, so return 0, since no file is open, so
        // does not need to be closed.
        return(0);
    }
    int result = file->ifclose();
    delete file;
    file = NULL;
    return(result);
}

/// Read up to size bytes from the file into the buffer.
/// \param file file to be read - IFILE is a pointer to an InputFile object
/// \param buffer pointer to memory at least size bytes big to write the
/// data into.
/// \param size number of bytes to be read
/// \return number of bytes read
inline unsigned int ifread(IFILE file, void * buffer, unsigned int size)
{
    if(file == NULL)
    {
        // No file was passed in, so 0 bytes were read.
        return(0);
    }
    return(file->ifread(buffer, size));
}

/// Get a character from the file.  Read a character from the internal
/// buffer, or if the end of the buffer has been reached, read from the
/// file into the buffer and return index 0.
/// \param file file to be read - IFILE is a pointer to an InputFile object
/// \return character that was read or EOF.
inline int ifgetc(IFILE file)
{
    if(file == NULL)
    {
        // return eof since there is no file.
        return(EOF);
    }
    return(file->ifgetc());
}

/// Get a line from the file.
/// \param file file to be read - IFILE is a pointer to an InputFile object
/// \param buffer the buffer into which data is to be placed
/// \param max the maximum size of the buffer, in bytes
/// \return true if the last character read was an EOF
inline bool ifgetline(IFILE file, void *buffer, size_t max)
{
    if(file == NULL)
    {
        // return eof since there is no file.
        return(true);
    }
    return(file->ifgetline(buffer, max));
}

/// Reset to the beginning of the file (cannot be done for stdin/stdout).
/// \param file file to be rewound - IFILE is a pointer to an InputFile object
inline void ifrewind(IFILE file)
{
    if(file == NULL)
    {
        return;
    }
    file->ifrewind();
}

/// Check to see if we have reached the EOF (returns 0 if not EOF).
/// \param file file to be checked - IFILE is a pointer to an InputFile object
/// \return 0 if not EOF, any other value means EOF.
inline int ifeof(IFILE file)
{
    if(file == NULL)
    {
        // No file, so that is considered to be EOF, so return 1.
        return(1);
    }
    return(file->ifeof());
}

/// Write the specified number of bytes from the specified buffer into the file.
/// \param file file to write to - IFILE is a pointer to an InputFile object
/// \param buffer buffer containing size bytes to write to the file.
/// \param size number of bytes to write
/// \return number of bytes written
inline unsigned int ifwrite(IFILE file, const void * buffer, unsigned int size)
{
    if(file == NULL)
    {
        // No file specified, so retun 0 bytes written.
        return(0);
    }
    return(file->ifwrite(buffer, size));
}

/// Get current position in the file.  Can be fed back into ifseek.
/// \param file file to perform tell on - IFILE is a pointer to an InputFile object
/// \return current position in the file, -1 indicates an error.
inline int64_t iftell(IFILE file)
{
    if(file == NULL)
    {
        return(-1);
    }
    return (file->iftell());
}

/// Seek to the specified position (result from an iftell), but cannot
/// be done for stdin/stdout.
/// \param file file to perform seek on - IFILE is a pointer to an InputFile object
/// \param offset offset into the file to move to (must be from a tell call)
/// \param origin can be any of the following:
/// Note: not all are valid for all filetypes.
///   SEEK_SET - Beginning of file
///   SEEK_CUR - Current position of the file pointer
///   SEEK_END - End of file
/// \return true on successful seek and false on a failed seek.
inline bool ifseek(IFILE file, int64_t offset, int origin)
{
    if(file == NULL)
    {
        // Could not see since no file was specified.
        return(false);
    }
    return (file->ifseek(offset, origin));
}

/// Write to a file using fprintf format.
/// \param file file to write to - IFILE is a pointer to an InputFile object
/// \param format printf format for writing, followed by parameters.
/// \return number of bytes written
int ifprintf(IFILE output, const char * format, ...);

/// Read a line from a file using streaming.  
/// Will not fail when the file hits EOF, so do not do: while(iFile >> iStr)
/// unless within your loop you check for ifeof and break.
/// Instead, do something like:
///    while(!iFile->ifeof() && iFile >> iStr)
/// \param stream file to read from - IFILE is a pointer to an InputFile object
/// \param str output string containing the line read from the file.
inline IFILE operator >> (IFILE stream, std::string &str)
{
    str.clear();
    int ch;
    // not safe... newline handling?
    while ((ch = stream->ifgetc())!=EOF && (ch != '\n')) str.push_back(ch);
    return stream;
}

/// Write to a file using streaming.
/// \param stream file to write to - IFILE is a pointer to an InputFile object
/// \param str string containing what should be written to the file.
inline InputFile& operator << (InputFile& stream, const std::string& str)
{
    unsigned int numExpected = str.length();
    unsigned int numWritten = 
        stream.ifwrite(str.c_str(), numExpected);
    if(numExpected != numWritten)
    {
        std::cerr << "Failed to stream to IFILE, expected " 
                  << numExpected << " but only wrote "
                  << numWritten << std::endl;
    }
    return(stream);
}

/// Write to a file using streaming.
/// \param stream file to write to - IFILE is a pointer to an InputFile object
/// \param str string containing what should be written to the file.
inline InputFile& operator << (InputFile& stream, const char* str)
{
    unsigned int numExpected = strlen(str);
    unsigned int numWritten = 
        stream.ifwrite(str, numExpected);
    if(numExpected != numWritten)
    {
        std::cerr << "Failed to stream to IFILE, expected " 
                  << numExpected << " but only wrote "
                  << numWritten << std::endl;
    }
    return(stream);
}


/// Write to a file using streaming.
/// \param stream file to write to - IFILE is a pointer to an InputFile object
/// \param num number that should be written to the file.
InputFile& operator << (InputFile& stream, double num);

/// Write to a file using streaming.
/// \param stream file to write to - IFILE is a pointer to an InputFile object
/// \param num number that should be written to the file.
InputFile& operator << (InputFile& stream, int num);

/// Write to a file using streaming.
/// \param stream file to write to - IFILE is a pointer to an InputFile object
/// \param num number that should be written to the file.
InputFile& operator << (InputFile& stream, unsigned int num);

/// Write to a file using streaming.
/// \param stream file to write to - IFILE is a pointer to an InputFile object
/// \param ch character that should be written to the file.
inline InputFile& operator << (InputFile& stream, char ch)
{
    unsigned int numWritten = 
        stream.ifwrite(&ch, 1);
    if(1 != numWritten)
    {
        std::cerr << "Failed to stream to IFILE, expected 1, but only wrote " 
                  << numWritten << std::endl;
    }
    return(stream);
}

#endif

