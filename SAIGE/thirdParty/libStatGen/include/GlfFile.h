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

#ifndef __GLF_FILE_H__
#define __GLF_FILE_H__

#include "InputFile.h"
#include "GlfHeader.h"
#include "GlfRefSection.h"
#include "GlfRecord.h"
#include "GlfStatus.h"

/// This class allows a user to easily read/write a GLF file.
class GlfFile
{
public:
    /// Enum for indicating whether to open the file for read or write.
    enum OpenType 
        {
            READ, ///< open for reading.
            WRITE ///< open for writing.
        };

    /// Default Constructor.
    GlfFile();

    /// Constructor that opens the specified file based on the specified mode
    /// (READ/WRITE).  Default is READ.
    /// \param filename name of the file to open.
    /// \param mode mode to use for opening the file (defaults to READ).
    GlfFile(const char* filename, OpenType mode = READ);

    /// Closes the file if there is one open, adding an end marker record
    /// if there is a previous section and one has not already been written.
    virtual ~GlfFile();
   
    /// Open a glf file for reading with the specified filename.
    /// \param  filename glf file to open for reading.
    /// \return true = success; false = failure.   
    bool openForRead(const char * filename);

    /// Open a glf file for reading with the specified filename and read the
    /// header into the specified header.
    /// \param  filename glf file to open for reading.
    /// \param  header header object to populate with the file's glf header.
    /// \return true = success; false = failure.   
    bool openForRead(const char * filename, GlfHeader& header);

    /// Open a glf file for writing with the specified filename.
    /// \param  filename glf file to open for writing.
    /// \param  compressed whether or not to compress the file, defaults to true
    /// \return true = success; false = failure.
    bool openForWrite(const char * filename, bool compressed = true);

    /// Close the file if there is one open, adding an end marker record
    /// if there is a previous section and one has not already been written.
    void close();

    /// Returns whether or not the end of the file has been reached.
    /// \return true = EOF; false = not eof.
    /// If the file is not open, true is returned.
    bool isEOF();
   
    /// Reads the header section from the file and stores it in
    /// the passed in header.
    /// \param  header header object to populate with the file's glf header.
    /// \return true = success; false = failure.
    bool readHeader(GlfHeader& header);
   
    /// Writes the specified header into the file.
    /// \param  header header object to write into the file.
    /// \return true = success; false = failure.
    bool writeHeader(GlfHeader& header);

    /// Gets the next reference section from the file & stores it in the
    /// passed in section, consuming records until a new section is found.
    /// \param  refSection object to populate with the file's next reference 
    ///                    section.
    /// \return true  = section was successfully set.
    ///         false = section was not successfully set.
    bool getNextRefSection(GlfRefSection& refSection);
   
    /// Write the reference section to the file, adding an end marker record
    /// if there is a previous section and one has not already been written.
    /// \param  refSection reference section to write to the file.
    /// \return true = succes; false = failure.
    bool writeRefSection(const GlfRefSection& refSection);

    /// Gets the nextrecord from the file & stores it in the
    /// passed in record.
    /// \param  record object to populate with the file's next record. 
    /// \return true  = record was successfully set.
    ///         false = record not successfully set or for the endMarker record.
    bool getNextRecord(GlfRecord& record);
   
    /// Writes the specified record into the file.
    /// \param record record to write to the file.
    /// \return true = success; false = failure.
    bool writeRecord(const GlfRecord& record);
   
    /// Return the number of records that have been read/written so far.
    /// \return number of records that have been read/written so far.
    uint32_t getCurrentRecordCount();

    /// Get the Status of the last call that sets status.
    /// To remain backwards compatable - will be removed later.
    inline GlfStatus::Status getFailure()
    {
        return(getStatus());
    }

    /// Get the Status of the last call that sets status.
    /// \return status of the last method that sets a status.
    inline GlfStatus::Status getStatus()
    {
        return(myStatus.getStatus());
    }

    /// Get the Status of the last call that sets status.
    /// \return status message of the last method that sets a status.
    inline const char* getStatusMessage()
    {
        return(myStatus.getStatusMessage());
    }

private:
    /// reset this file including all its attributes.
    void resetFile();

    /// Pointer to the file
    IFILE  myFilePtr;

    /// Flag to indicate if a file is open for reading.
    bool myIsOpenForRead;
    /// Flag to indicate if a file is open for writing.
    bool myIsOpenForWrite;

    /// End marker that is inserted when writing files if a new section
    /// is specified without one or if the file is closed without writing
    /// an endMarker.
    GlfRecord myEndMarker;

    /// Track the state of this file as to what it is expecting to read next.
    enum EXPECTED_SECTION
    {
        HEADER,
        REF_SECTION,
        RECORD
    } myNextSection;
    
    /// Keep count of the number of records that have been read/written so far.
    uint32_t myRecordCount;

    /// The status of the last GlfFile command.
    GlfStatus myStatus;
};


class GlfFileReader : public GlfFile
{
public:

    /// Default Constructor.
    GlfFileReader();

    /// Constructor that opens the specified file for read.
     /// \param filename file to open for reading.
   GlfFileReader(const char* filename);

    virtual ~GlfFileReader();
};


class GlfFileWriter : public GlfFile
{
public:
    /// Default Constructor.
    GlfFileWriter();

    /// Constructor that opens the specified file for write.
    /// \param filename file to open for writing.
    GlfFileWriter(const char* filename);

    virtual ~GlfFileWriter();
};

#endif
