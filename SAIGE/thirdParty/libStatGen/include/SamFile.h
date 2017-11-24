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

#ifndef __SAM_FILE_H__
#define __SAM_FILE_H__

#include "SamStatus.h"
#include "InputFile.h"
#include "SamFileHeader.h"
#include "SamRecord.h"
#include "GenericSamInterface.h"
#include "BamIndex.h"
#include "SamStatistics.h"

/// Allows the user to easily read/write a SAM/BAM file.
/// The SamFile class contains additional functionality that allows a user
/// to read specific sections of sorted & indexed BAM files.  In order to
/// take advantage of this capability, the index file must be read prior to
/// setting the read section.  This logic saves the time of having to read
/// the entire file and takes advantage of the seeking capability of BGZF.
class SamFile
{
public:
    /// Enum for indicating whether to open the file for read or write.
    enum OpenType {
        READ, ///< open for reading.
        WRITE ///< open for writing.
    };
    
    
    /// Enum for indicating the type of sort expected in the file.
    enum SortedType {
        UNSORTED = 0, ///< file is not sorted.
        FLAG,         ///< SO flag from the header indicates the sort type.
        COORDINATE,   ///< file is sorted by coordinate.
        QUERY_NAME    ///< file is sorted by queryname.
    };
    
    /// Default Constructor, initializes the variables, but does not open
    /// any files.
    SamFile();

    /// Constructor that sets the error handling type.
    /// \param errorHandlingType how to handle errors.
    SamFile(ErrorHandler::HandlingType errorHandlingType);

    /// Constructor that opens the specified file based on the specified mode
    /// (READ/WRITE), aborts if the file could not be opened.
    /// \param filename name of the file to open.
    /// \param mode mode to use for opening the file.
    SamFile(const char* filename, OpenType mode);

    /// Constructor that opens the specified file based on the specified mode
    /// (READ/WRITE) and handles errors per the specified handleType.
    /// \param filename name of the file to open.
    /// \param mode mode to use for opening the file.
    /// \param errorHandlingType how to handle errors.
    SamFile(const char* filename, OpenType mode,
            ErrorHandler::HandlingType errorHandlingType);

    /// Constructor that opens the specified file based on the specified mode
    /// (READ/WRITE) and reads the header, aborts if the file could not be
    /// opened or the header not read.
    /// \param filename name of the file to open.
    /// \param mode mode to use for opening the file.
    /// \param header to read into or write from
    SamFile(const char* filename, OpenType mode, SamFileHeader* header);

    /// Constructor that opens the specified file based on the specified mode
    /// (READ/WRITE) and reads the header, handling errors per the specified
    /// handleType.
    /// \param filename name of the file to open.
    /// \param mode mode to use for opening the file.
    /// \param errorHandlingType how to handle errors.
    /// \param header to read into or write from
    SamFile(const char* filename, OpenType mode,
            ErrorHandler::HandlingType errorHandlingType,
            SamFileHeader* header);

    /// Destructor
    virtual ~SamFile();
   
    /// Open a sam/bam file for reading with the specified filename,
    /// determing the type of file and SAM/BAM  by reading the file
    /// (if not stdin).
    /// \param  filename the sam/bam file to open for reading.
    /// \param header to read into or write from (optional)
    /// \return true = success; false = failure.
    bool OpenForRead(const char * filename, SamFileHeader* header = NULL);

    /// Open a sam/bam file for writing with the specified filename,
    /// determining SAM/BAM from the extension (.bam = BAM).
    /// \param  filename the sam/bam file to open for writing.
    /// \param header to read into or write from (optional)
    /// \return true = success; false = failure.
    bool OpenForWrite(const char * filename, SamFileHeader* header = NULL);

    /// Read the specified bam index file.  It must be read prior to setting a
    /// read section, for seeking and reading portions of a bam file.
    /// \param filename the name of the bam index file to be read.
    /// \return true = success; false = failure.
    bool ReadBamIndex(const char * filename);

    /// Read the bam index file using the BAM filename as a base. 
    /// It must be read prior to setting a read section, for seeking
    /// and reading portions of a bam file.
    /// Must be read after opening the BAM file since it uses the
    /// BAM filename as a base name for the index file.
    /// First it tries filename.bam.bai. If that fails, it tries
    /// it without the .bam extension, filename.bai.
    /// \return true = success; false = failure.
    bool ReadBamIndex();

    /// Sets the reference to the specified genome sequence object.
    /// \param reference pointer to the GenomeSequence object.
    void SetReference(GenomeSequence* reference);

    /// Set the type of sequence translation to use when reading
    /// the sequence.  Passed down to the SamRecord when it is read.  
    /// The default type (if this method is never called) is
    /// NONE (the sequence is left as-is).
    /// \param translation type of sequence translation to use.
    void SetReadSequenceTranslation(SamRecord::SequenceTranslation translation);

    /// Set the type of sequence translation to use when writing
    /// the sequence.  Passed down to the SamRecord when it is written.
    /// The default type (if this method is never called) is
    /// NONE (the sequence is left as-is).
    /// \param translation type of sequence translation to use.
    void SetWriteSequenceTranslation(SamRecord::SequenceTranslation translation);

    /// Close the file if there is one open.
    void Close();

    /// Returns whether or not the file has been opened successfully.
    /// \return true = open; false = not open.
    bool IsOpen();
   
    /// Returns whether or not the end of the file has been reached.
    /// \return true = EOF; false = not eof.
    /// If the file is not open, true is returned.
    bool IsEOF();
   
    /// Returns whether or not the file has been opened for streaming 
    /// input/output.
    /// \return true = stream; false = not a stream.
    bool IsStream();

    /// Reads the header section from the file and stores it in
    /// the passed in header.
    /// \return true = success; false = failure.
    bool ReadHeader(SamFileHeader& header);
   
    /// Writes the specified header into the file.
    /// \return true = success; false = failure.
    bool WriteHeader(SamFileHeader& header);

    /// Reads the next record from the file & stores it in the passed in record.
    ///
    /// If it is an indexed BAM file and SetReadSection was called, 
    /// only alignments in the section specified by SetReadSection are read.
    /// If they all have already been read, this method returns false.
    ///
    /// Validates that the record is sorted according to the value set by
    /// setSortedValidation. No sorting validation is done if specified to be
    /// unsorted, or setSortedValidation was never called.
    /// \return true  = record was successfully set (and sorted if applicable),
    ///                false = record was not successfully set 
    ///                (or not sorted as expected).
    bool ReadRecord(SamFileHeader& header, SamRecord& record);
   
    /// Writes the specified record into the file.
    /// Validates that the record is sorted according to the value set by
    /// setSortedValidation. No sorting validation is done if specified to
    /// be unsorted, or setSortedValidation was never called.  Returns false
    /// and does not write the record if the record was not properly sorted.
    /// \return true = success; false = failure.
    bool WriteRecord(SamFileHeader& header, SamRecord& record);
   
    /// Set the flag to validate that the file is sorted as it is read/written.
    /// Must be called after the file has been opened.
    /// Sorting validation is reset everytime SetReadPosition is called since
    /// it can jump around in the file.
    /// \param sortType specifies the type of sort to be checked for.
    void setSortedValidation(SortedType sortType);

    /// Return the number of records that have been read/written so far.
    uint32_t GetCurrentRecordCount();

    /// Deprecated, get the Status of the last call that sets status.
    /// To remain backwards compatable - will be removed later.
    inline SamStatus::Status GetFailure()
    {
        return(GetStatus());
    }

    /// Get the Status of the last call that sets status.
    inline SamStatus::Status GetStatus()
    {
        return(myStatus.getStatus());
    }

    /// Get the Status Message of the last call that sets status.
    inline const char* GetStatusMessage()
    {
        return(myStatus.getStatusMessage());
    }

    /// Sets which reference id (index into the BAM list of reference
    /// information) of the BAM file should be read.  The records
    /// for that reference id will be retrieved on each ReadRecord call.  
    /// Reference ids start at 0, and -1 indicates reads with no reference.
    /// When all records have been retrieved for the specified reference id, 
    /// ReadRecord will return failure until a new read section is set.
    /// Must be called only after the file has been opened for reading.
    /// Sorting validation is reset everytime SetReadPosition is called since
    /// it can jump around in the file.
    /// \param  refID the reference ID of the records to read from the file.
    /// \return true = success; false = failure.
    bool SetReadSection(int32_t refID);

    /// Sets which reference name of the BAM file should be read.  The records
    /// for that reference name will be retrieved on each ReadRecord call.
    /// Specify "" or "*" to read records not associated with a reference.
    /// When all records have been retrieved for the specified reference name,
    /// ReadRecord will return failure until a new read section is set.
    /// Must be called only after the file has been opened for reading.
    /// Sorting validation is reset everytime SetReadPosition is called since
    /// it can jump around in the file.
    /// \param  refName the reference name of the records to read from the file.
    /// \return true = success; false = failure.
    bool SetReadSection(const char* refName);

    /// Sets which reference id (index into the BAM list of reference
    /// information) & start/end positions of the BAM file should be read.
    /// The records for that reference id and positions will be retrieved on
    /// each ReadRecord call.  Reference ids start at 0, and -1 indicates
    /// reads with no reference.  When all records have been retrieved for the
    /// specified reference id, ReadRecord will return failure until a new read
    /// section is set.  Must be called only after the file has been opened
    /// for reading.  Sorting validation is reset everytime SetReadPosition is
    /// called since it can jump around in the file.
    /// \param  refID the reference ID of the records to read from the file.
    /// \param  start inclusive 0-based start position of records that should be read for this refID.
    /// \param  end exclusive 0-based end position of records that should be read for this refID.
    /// \param overlap When true (default), return reads that just overlap the region; when false, only return reads that fall completely within the region
    /// \return true = success; false = failure.   
    bool SetReadSection(int32_t refID, int32_t start, int32_t end, 
                        bool overlap = true);

    /// Sets which reference name & start/end positions of the BAM file should
    /// be read.  The records for this reference name & positions will be
    /// retrieved on each ReadRecord call.  Specify "" or "*" to indicate
    /// reads with no reference.  When all records have been retrieved for
    /// the specified section, ReadRecord will return failure until a new read
    /// section is set.  Must be called only after the file has been opened for
    /// reading.  Sorting validation is reset everytime SetReadSection is
    /// called since it can jump around in the file.
    /// \param  refName the reference name of the records to read from the file.
    /// \param  start inclusive 0-based start position of records that should be read for this refID.
    /// \param  end exclusive 0-based end position of records that should be read for this refID.
    /// \param overlap When true (default), return reads that just overlap the region; when false, only return reads that fall completely within the region
    /// \return true = success; false = failure.   
    bool SetReadSection(const char* refName, int32_t start, int32_t end, 
                        bool overlap = true);

    /// Specify which reads should be returned by ReadRecord.
    /// Reads will only be returned by ReadRecord that contain the specified
    /// required flags and that do not contain any of the specified excluded
    /// flags.  ReadRecord will continue to read from the file until a record
    /// that complies with these flag settings is found or until the end of the
    /// file/region.
    /// \param requiredFlags flags that are required to be in records
    /// returned by ReadRecord (set to 0x0 if there are no required flags).
    /// \param excludedFlags flags that are required to not be in records
    /// returned by ReadRecord (set to 0x0 if there are no excluded flags).
    void SetReadFlags(uint16_t requiredFlags, uint16_t excludedFlags);

    /// Get the number of mapped reads in the specified reference id.  
    /// Returns -1 for out of range refIDs.
    /// \param refID reference ID for which to extract the number of mapped reads.
    /// \return number of mapped reads for the specified reference id.
    int32_t getNumMappedReadsFromIndex(int32_t refID);

    /// Get the number of unmapped reads in the specified reference id.  
    /// Returns -1 for out of range refIDs.
    /// \param refID reference ID for which to extract the number of unmapped reads.
    /// \return number of unmapped reads for the specified reference id.
    int32_t getNumUnMappedReadsFromIndex(int32_t refID);

    /// Get the number of mapped reads in the specified reference name.
    /// Returns -1 for unknown reference names.
    /// \param refName reference name for which to extract the number of mapped reads.
    /// \param header header object containing the map from refName to refID
    /// \return number of mapped reads for the specified reference name.
    int32_t getNumMappedReadsFromIndex(const char* refName,
                                       SamFileHeader& header);

    /// Get the number of unmapped reads in the specified reference name.
    /// Returns -1 for unknown reference names.
    /// \param refName reference name for which to extract the number of unmapped reads.
    /// \param header header object containing the map from refName to refID
    /// \return number of unmapped reads for the specified reference name.
    int32_t getNumUnMappedReadsFromIndex(const char* refName,
                                         SamFileHeader& header);

    /// Returns the number of bases in the passed in read that overlap the
    /// region that is currently set.  Overlapping means that the bases occur
    /// in both the read and the reference as either matches or mismatches.
    /// This does not count insertions, deletions, clips, pads, or skips.
    /// \param samRecord to check for overlapping bases.
    /// \return number of bases that overlap region that is currently set.
    uint32_t GetNumOverlaps(SamRecord& samRecord);

    /// Whether or not statistics should be generated for this file.
    /// The value is carried over between files and is not reset, but
    /// the statistics themselves are reset between files.
    /// \param genStats set to true if statistics should be generated, false if not.
    void GenerateStatistics(bool genStats);

    /// Return the bam index if one has been opened.
    /// \return const pointer to the bam index, or null if one has not been opened.
    const BamIndex* GetBamIndex();

    /// Get the current file position.
    /// \return current position in the file.
    inline int64_t GetCurrentPosition()
    {
        return(iftell(myFilePtr));
    }
    
    /// Turn off file read buffering.
    inline void DisableBuffering()
    {
        if(myFilePtr != NULL)
        {
            myFilePtr->disableBuffering();
        }
    }

    /// Print the statistics that have been recorded due to a call to
    /// GenerateStatistics.
    inline void PrintStatistics() {if(myStatistics != NULL) myStatistics->print();}

protected:
    void init();
    void init(const char* filename, OpenType mode, SamFileHeader* header);

    /// Resets the file prepping for a new file.
    void resetFile();

    /// Validate that the record is sorted compared to the previously read
    /// record if there is one, according to the specified sort order.
    /// If the sort order is UNSORTED, true is returned.
    /// Sorting validation is reset everytime SetReadPosition is called since
    /// it can jump around in the file.
    bool validateSortOrder(SamRecord& record, SamFileHeader& header);
   
    // Return the sort order as defined by the header.  If it is undefined
    // or set to an unknown value, UNSORTED is returned.
    SortedType getSortOrderFromHeader(SamFileHeader& header);

    bool processNewSection(SamFileHeader &header);

    // Check if there is more to read in the current chunk, if not,
    // move to the next chunk.
    // If no sections are specified or it successfully found a chunk to read,
    // return true.
    // Sets the status and returns false if it was unable to move to a new chunk
    // or there are no more chunks to read, otherwise returns true.
    bool ensureIndexedReadPosition();

    // Check whether or not the record falls within the specified section. 
    // If no sections are specified or this read falls within the
    // specified sections, return true.
    // If it does not, return false.
    // If the record position indicates there will be no more records within the
    // region, return false AND set the sam status to indicate NO_MORE_RECS.
    bool checkRecordInSection(SamRecord& record);

    IFILE  myFilePtr;
    GenericSamInterface* myInterfacePtr;

    /// Flag to indicate if a file is open for reading.
    bool myIsOpenForRead;
    /// Flag to indicate if a file is open for writing.
    bool myIsOpenForWrite;
    /// Flag to indicate if a header has been read/written - required before
    /// being able to read/write a record.
    bool myHasHeader;

    SortedType mySortedType;

    /// Previous values used for checking if the file is sorted.
    int32_t myPrevCoord;
    int32_t myPrevRefID;
    String myPrevReadName;

    /// Keep a count of the number of records that have been read/written so far.
    uint32_t myRecordCount;

    /// Pointer to the statistics for this file.
    SamStatistics* myStatistics;
   
    /// The status of the last SamFile command.
    SamStatus myStatus;

    /// Values for reading Sorted BAM files via the index.
    bool myIsBamOpenForRead;
    bool myNewSection;
    // whether to return reads that overlap (true) the section or
    // are fully enclosed (false) in the section.
    bool myOverlapSection;
    int32_t myRefID;
    int32_t myStartPos;
    int32_t myEndPos;
    uint64_t myCurrentChunkEnd;
    SortedChunkList myChunksToRead;
    BamIndex* myBamIndex;

    GenomeSequence* myRefPtr;
    SamRecord::SequenceTranslation myReadTranslation;
    SamRecord::SequenceTranslation myWriteTranslation;
    
    std::string myRefName;

private:
    bool    myAttemptRecovery;

    uint16_t myRequiredFlags;
    uint16_t myExcludedFlags;

public:

    bool attemptRecoverySync(bool (*checkSignature)(void *data) , int length);

    void setAttemptRecovery(bool flag = false)
    {
        myAttemptRecovery = flag;
    }

};


/// Child class of SamFile for reading files.
class SamFileReader : public SamFile
{
public:

    /// Default Constructor.
    SamFileReader();

    /// Constructor that opens the specified file for read.
    SamFileReader(const char* filename);

    /// Constructor that opens the specified file for read.
    SamFileReader(const char* filename,
                  ErrorHandler::HandlingType errorHandlingType);

    /// Constructor that opens the specified file for read and reads
    /// the header from the file.
    SamFileReader(const char* filename,
                  SamFileHeader* header);

    /// Constructor that opens the specified file for read and reads
    /// the header from the file.
    SamFileReader(const char* filename,
                  ErrorHandler::HandlingType errorHandlingType,
                  SamFileHeader* header);

    virtual ~SamFileReader();
};


/// Child class of SamFile for writing files.
class SamFileWriter : public SamFile
{
public:
    /// Default Constructor.
    SamFileWriter();

    /// Constructor that opens the specified file for write.
    SamFileWriter(const char* filename);

    /// Constructor that opens the specified file for write.
    SamFileWriter(const char* filename,
                  ErrorHandler::HandlingType errorHandlingType);

    /// Constructor that opens the specified file for write and write
    /// the specified header into the file.
    SamFileWriter(const char* filename,
                  SamFileHeader* header);

    /// Constructor that opens the specified file for write and write
    /// the specified header into the file.
    SamFileWriter(const char* filename,
                  ErrorHandler::HandlingType errorHandlingType,
                  SamFileHeader* header);

    virtual ~SamFileWriter();
};

#endif
