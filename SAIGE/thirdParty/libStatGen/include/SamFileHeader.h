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

#ifndef __SAM_FILE_HEADER_H__
#define __SAM_FILE_HEADER_H__

#include <map>
#include <stdint.h>

#include "SamReferenceInfo.h"
#include "SamHeaderHD.h"
#include "SamHeaderSQ.h"
#include "SamHeaderRG.h"
#include "SamHeaderPG.h"

/// This class allows a user to get/set the fields in a SAM/BAM Header.
/// Sam/Bam headers contain comments and multiple SamHeaderRecords 
/// (HD, SQs, RGs, PGs) comprised of tag/value pairs with each tag only
/// appearing once within a specific record.
class SamFileHeader
{
public:
    SamFileHeader();
    ~SamFileHeader();

    /////////////////////////////
    /// @name  Copying a Header
    /// These methods are ways of copying the contents of one header into
    /// another one.
    //@{

    /// Copy Constructor copies the specified header into this one.
    SamFileHeader(const SamFileHeader& header);

    /// Overload operator = to copy the passed in header into this header.
    SamFileHeader & operator = (const SamFileHeader& header);

    /// Copy method copies the passed in header into this header.
    /// Returns true if at least one header line was successfully copied.
    bool copy(const SamFileHeader& header);
    //@}

    /// Initialize the header.
    void resetHeader();

    /////////////////////////////
    /// @name  Get the Entire Header
    /// Get the entire header as a single string.
    //@{

    /// Set the passed in string to the entire header string, clearing its
    /// current contents.
    /// \return true if successfully set (even if set to "")
    bool getHeaderString(std::string& header) const;

    //@}

    /// Get the reference ID for the specified reference name (chromosome).
    /// If addID is set to true, a reference id will be created for the
    /// referenceName if one does not already exist.  If addID is set to
    /// false (default), it will return SamReferenceInfo::NO_REF_ID.
    int getReferenceID(const String & referenceName, bool addID = false);

    /// Get the reference ID for the specified reference name (chromosome).
    /// If addID is set to true, a reference id will be created for the
    /// referenceName if one does not already exist.  If addID is set to
    /// false (default), it will return SamReferenceInfo::NO_REF_ID.
    int getReferenceID(const char* referenceName, bool addID = false);

    /// Return the reference name (chromosome) for the specified reference id.
    const String & getReferenceLabel(int id) const;

    /// Get the Reference Information
    const SamReferenceInfo& getReferenceInfo() const;

    // Get the Reference Information for updating separately when reading
    // BAMs...should only be called by BamInterface.
    SamReferenceInfo& getReferenceInfoForBamInterface();

    ////////////////////////////////////////////////////////////////////////
    // Set Values in the header
    ////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////
    /// @name  Adding an entire header/comment line.
    /// These methods are ways of adding an entire header line at once.
    //@{

    /// Add a header line that is just one tag with a const char* value.
    /// Note: This method will only do one tag per type on a line, so if a
    /// type has multiple tags, the whole line needs to be added at once,
    /// and a different method should be used.
    bool addHeaderLine(const char* type, const char* tag, const char* value); 

    /// Add a header line that is already preformatted in a const char*.
    /// Returns true if at least one header line was successfully added.
    bool addHeaderLine(const char* headerLine);

    /// Add a header that is already preformatted in a const char*.
    /// Returns true if at least one header line was successfully added.
    bool addHeader(const char* header);

    /// Add the specified comment to the header (do not include "@CO" or "\n").
    /// \return true if successfully added, false if not.
    bool addComment(const char* comment);

    //@}

    
    /////////////////////////////////////////
    /// @name  Set/Add/Remove a Single Tag
    /// The passed in tag should be the two character SAM tag as defined
    /// in the SAM spec.  A tag is removed from the header record by setting
    /// it to "". For the SQ and RG header types, the key tags (SN for SQ
    /// and ID for RG) may not be modified or removed once set. This is
    /// because these values are used as a lookup key for the header record, 
    /// so the entire record must be removed.
    //@{

//     // Set the specified header type tag to the specified value in the 
//     // header with the specified keyID.  keyID must be specified when
//     // type = SQ, RG, or PG.
//     bool setTag(SamHeaderRecord::SamHeaderRecordType type, const char* tag,
//                 const char* value, const char* keyID = NULL);

    /// Set the specified tag to the specified value in the HD header, remove
    /// the tag by specifying value="".
    /// \return true if the tag was successfully set, false if not.
    bool setHDTag(const char* tag, const char* value);

    /// Set the specified tag to the specified value in the SQ header with
    /// the specified name, remove the tag by specifying value="".  If the
    /// header does not yet exist, the tag must be "LN" and the header is added
    /// with the specified LN value and the SN value passed in name.  
    /// The SN & LN tags may not be modified or removed after they are
    /// set unless the entire record is deleted.
    /// \return true if the tag was successfully set, false if not.
    bool setSQTag(const char* tag, const char* value, const char* name);

    /// Set the specified tag to the specified value in the RG header with
    /// the specified id, remove the tag by specifying value="".  If the
    /// header does not yet exist, the header is added and so is the ID tag
    /// with the value set to the passed in id.  The ID tag may not be 
    /// modified or removed after it is set unless the entire record is deleted.
    /// \return true if the tag was successfully set, false if not.
    bool setRGTag(const char* tag, const char* value, const char* id);

    /// Set the specified tag to the specified value in the PG header with
    /// the specified id, remove the tag by specifying value="".  If the
    /// header does not yet exist, the header is added and so is the ID tag
    /// with the value set to the passed in id.  The ID tag may not be 
    /// modified or removed after it is set unless the entire record is deleted.
    /// \return true if the tag was successfully set, false if not.
    bool setPGTag(const char* tag, const char* value, const char* id);

    //@}

    /////////////////////////////////////////
    /// @name  Add an Already Setup SamHeaderRecord
    /// NOTE: These methods add a pointer to the passed in record.
    /// The header record will be deleted when it's cleaned up from this header.
    /// NOTE: Do NOT delete the passed in record, the SamFileHeader class
    /// takes care of that itself.
    //@{

    /// Add the HD record to the header.
    /// Note: it adds a pointer to the passed in header record.  The header
    /// record will be deleted when it is cleaned up from this header.
    /// \return true if the record was successfully added, false otherwise.
    bool addHD(SamHeaderHD* hd);

    /// Add the SQ record to the header.
    /// Note: it adds a pointer to the passed in header record.  The header
    /// record will be deleted when it is cleaned up from this header.
    /// \return true if the record was successfully added, false otherwise.
    bool addSQ(SamHeaderSQ* sq);

    /// Add the RG record to the header.
    /// Note: it adds a pointer to the passed in header record.  The header
    /// record will be deleted when it is cleaned up from this header.
    /// \return true if the record was successfully added, false otherwise.
    bool addRG(SamHeaderRG* rg);

    /// Add the PG record to the header.
    /// Note: it adds a pointer to the passed in header record.  The header
    /// record will be deleted when it is cleaned up from this header.
    /// \return true if the record was successfully added, false otherwise.
    bool addPG(SamHeaderPG* pg);

    /// Add a copy of the specified header record to the header.
    /// Note: it creates a new header record that is identical to the specified
    /// one and adds it to the header.  The passed in pointer will not be
    /// deleted due to this.
    /// \return true if the record was successfully added, false otherwise.
    bool addRecordCopy(const SamHeaderRecord& hdrRec);

    //@}

    ////////////////////////////////////////////////////////////////////////
    /// @name  Remove an Entire Header Record
    //@{

    /// Remove the HD record.
    /// \return true if successfully removed or did not exist, false if
    /// the record still exists.
    bool removeHD();

    /// Remove SQ record with the specified key.
    /// NOTE: Does not remove it from the BAM index.
    /// \return true if successfully removed or did not exist, false if
    /// the record still exists.
    bool removeSQ(const char* name);

    /// Remove RG record with the specified key.
    /// \return true if successfully removed or did not exist, false if
    /// the record still exists.
    bool removeRG(const char* id);

    /// Remove PG record with the specified key.
    /// \return true if successfully removed or did not exist, false if
    /// the record still exists.
    bool removePG(const char* id);

    //@}

    ////////////////////////////////////////////////////////////////////////
    /// @name  Get a Specific Tag
    /// These methods return the value associated with the specified tag.
    /// If the tag does not exist in the record "" is returned.
    ///
    /// For SQ, RG, and PG the value returned is for the tag associated with
    /// the specified key (name/id). If a record with that key does not exist
    /// or if the tag does not exist for the record with that key, "" is 
    /// returned.
    //@{

    /// Returns the value associated with the specified HD tag, returning "" if
    /// the tag does not exist in the header.
    const char* getHDTagValue(const char* tag);

    /// Get the value associated with the specified tag on the SQ line with
    /// the specified sequence name, returning "" if the tag or key does
    /// not exist.
    const char* getSQTagValue(const char* tag, const char* name);

    /// Get the value associated with the specified tag on the RG line with
    /// the specified read group identifier, returning "" if the tag or key does
    /// not exist.
    const char* getRGTagValue(const char* tag, const char* id);

    /// Get the value associated with the specified tag on the RG line with
    /// the specified id, returning "" if the tag or key does
    /// not exist.
    const char* getPGTagValue(const char* tag, const char* id);

    //@}

    /// Get the number of SQ objects.
    int getNumSQs();

    /// Get the number of RG objects.
    int getNumRGs();

    /// Get the number of PG objects.
    int getNumPGs();

    ////////////////////////////////////////////////////////////////////////
    /// @name  Get a Specific Header Record
    /// These methods return a reference to the specific record that was
    /// requested, returning NULL if that record does not exist in the header.
    ///
    /// The returned record can be modified to add/remove some tags.
    /// Since a reference is returned, the SamHeaderFile automatically 
    /// reflects these changes.
    //@{

    /// Get the HD object, returning NULL if there is no HD record.
    SamHeaderHD* getHD();

    /// Get the SQ object with the specified sequence name, returning NULL
    /// if there is no SQ object with that key.
    SamHeaderSQ* getSQ(const char* name);

    /// Get the RG object with the specified read group identifier, returning
    /// NULL if there is no RG object with that key..
    SamHeaderRG* getRG(const char* id);

    /// Get the PG object with the specified id, returning NULL
    /// if there is no PG object with that key..
    SamHeaderPG* getPG(const char* id);

    //@}

//     //////////////////////////////////
//     // Set methods for header fields.
//     bool setVersion(const char* version);
//     bool setSortOrder(const char* sortOrder);
//     bool addSequenceName(const char* sequenceName);
//     bool setSequenceLength(const char* keyID, int sequenceLength);
//     bool setGenomeAssemblyId(const char* keyID, const char* genomeAssemblyId);
//     bool setMD5Checksum(const char* keyID, const char* md5sum);
//     bool setURI(const char* keyID, const char* uri);
//     bool setSpecies(const char* keyID, const char* species);
//     bool addReadGroupID(const char* readGroupID);
//     bool setSample(const char* keyID, const char* sample);
//     bool setLibrary(const char* keyID, const char* library);
//     bool setDescription(const char* keyID, const char* description);
//     bool setPlatformUnit(const char* keyID, const char* platform);
//     bool setPredictedMedianInsertSize(const char* keyID, const char* isize);
//     bool setSequencingCenter(const char* keyID, const char* center);
//     bool setRunDate(const char* keyID, const char* runDate);
//     bool setTechnology(const char* keyID, const char* technology);
//     bool addProgram(const char* programID);
//     bool setProgramVersion(const char* keyID, const char* version);
//     bool setCommandLine(const char* keyID, const char* commandLine);
    
//     ///////////////////////////////////
//     // Get methods for header fields.
//     // Returns the number of SQ entries in the header.
//     int32_t getSequenceDictionaryCount();

    /// Return the Sort Order value that is set in the Header, returning ""
    /// if this field does not exist.
    const char* getSortOrder();


    /// DEPRECATED
    const char* getTagSO();

    /////////////////////////////
    /// @name  Get the Header Record/Comment/Line by Record/Comment/Line
    /// These methods iterate through the header.
    /// NOTE: both getNextHeaderRecord and getNextHeaderLine increment the
    /// same iterator.  getNextHeaderRecord that takes a header type
    /// uses the same iterator as the getNextXXRecord with that type.
    /// Otherwise the iterators are independent.
    //@{

    /// Get the next SQ header record.  After all SQ headers have been
    /// retrieved, NULL is returned until a reset is called.
    /// Independent from getNextHeaderRecord, getNextHeaderLine and the
    /// other getNextXXRecord methods and the associated reset methods.
    SamHeaderRecord* getNextSQRecord();

    /// Get the next RG header record.  After all RG headers have been
    /// retrieved, NULL is returned until a reset is called.
    /// Independent from getNextHeaderRecord, getNextHeaderLine and the
    /// other getNextXXRecord methods and the associated reset methods.
    SamHeaderRecord* getNextRGRecord();

    /// Get the next PG header record.  After all PG headers have been
    /// retrieved, NULL is returned until a reset is called.
    /// Independent from getNextHeaderRecord, getNextHeaderLine and the
    /// other getNextXXRecord methods and the associated reset methods.
    SamHeaderRecord* getNextPGRecord();

    /// Reset to the beginning of the header records so the next call
    /// to getNextSQRecord returns the first SQ header record.
    void resetSQRecordIter();

    /// Reset to the beginning of the header records so the next call
    /// to getNextRGRecord returns the first RG header record.
    void resetRGRecordIter();

    /// Reset to the beginning of the header records so the next call
    /// to getNextPGRecord returns the first PG header record.
    void resetPGRecordIter();

    /// Get the next header record of the specified type starting from the
    /// specified index and update the index.
    /// After all headers of that type have been retrieved,
    /// NULL is returned until a reset is called for that type.
    SamHeaderRecord* getNextHeaderRecord(uint32_t& index, 
                                         SamHeaderRecord::SamHeaderRecordType headerType);

    /// Get the next header record, but not comment line.  After all headers
    /// have been retrieved, NULL is returned until a reset is called.
    /// NOTE: both getNextHeaderRecord and getNextHeaderLine increment the
    /// same iterator.
    SamHeaderRecord* getNextHeaderRecord();

    /// Set the passed in string to the next header line, overwritting
    /// the passed in string.  If there are no more header lines or there
    /// is an error, false is returned and the passed in string is set to ""
    /// until a rest is called.
    /// NOTE: both getNextHeaderRecord and getNextHeaderLine increment the
    /// same iterator.
    bool getNextHeaderLine(std::string &headerLine);

    /// Reset to the beginning of the header records so the next call
    /// to getNextHeaderRecord returns the first header line.
    void resetHeaderRecordIter();
   
    /// Append all of the comment lines to the specified string.
    void appendCommentLines(std::string &commentLines);

    /// Returns the comment on the next comment line.  Returns "" if all comment
    /// lines have been returned, until resetCommentIter is called.
    const char* getNextComment();

    /// Resets to the beginning of the comments so getNextComment returns
    /// the first comment.
    void resetCommentIter();

    //@}


    /// Get the failure message if a method returned failure.
    const char* getErrorMessage()  { return(myErrorMessage.c_str()); }

    static const std::string EMPTY_RETURN;

private:
    // Parse the header string. 
    bool parseHeader(String& header);

    // Parse the specified line of the header.
    bool parseHeaderLine(const String& headerLine);

    // Set the passed in string to the header line at the specified index.
    // It does NOT clear the current contents of header.
    bool getHeaderLine(unsigned int index, std::string& header) const;

    int16_t makeKey(char ch1, char ch2)
    {
        return((ch1 << 8) + ch2);
    }

    // Only one HD type is allowed per file.
    SamHeaderHD* myHD;

    // There can be multiple SQ Types, indexed by SN.
    StringHash mySQs;

    // There can be multiple RG Types, indexed by ID.
    StringHash myRGs;

    // There can be multiple PG types, indexed by ID.
    StringHash myPGs;

    // Reference Name information
    SamReferenceInfo myReferenceInfo;

    // Vector of comments
    std::vector<std::string> myComments;

    std::vector<SamHeaderRecord*> myHeaderRecords;

    std::string myErrorMessage;

    uint32_t myCurrentSQIndex;

    uint32_t myCurrentRGIndex;

    uint32_t myCurrentPGIndex;

    uint32_t myCurrentHeaderIndex;

    uint32_t myCurrentCommentIndex;
};

#endif

