/*
 *  Copyright (C) 2010-2011  Regents of the University of Michigan
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

#ifndef __SAM_RECORD_H__
#define __SAM_RECORD_H__

#include <stdint.h>

#include "GenomeSequence.h"
#include "SamStatus.h"
#include "LongHash.h"
#include "MathVector.h"
#include "StringArray.h"
#include "IntArray.h"
#include "SamFileHeader.h"
#include "CigarRoller.h"

/// Structure of a BAM record.
struct bamRecordStruct
{
public:
    int32_t      myBlockSize;
    int32_t      myReferenceID;
    int32_t      myPosition;
    uint32_t     myReadNameLength : 8, myMapQuality : 8, myBin : 16;
    uint32_t     myCigarLength : 16, myFlag : 16;
    int32_t      myReadLength;
    int32_t      myMateReferenceID;
    int32_t      myMatePosition;
    int32_t      myInsertSize;             // Outer fragment length
    char  myData[1];
};


/// Class providing an easy to use interface to get/set/operate on the
/// fields in a SAM/BAM record. 
class SamRecord
{
public:
    /// Enum containing the settings on how to translate the sequence if a
    /// reference is available.  If no reference is available, no translation
    /// is done.
    enum SequenceTranslation { 
        NONE,   ///< Leave the sequence as is.
        EQUAL,  ///< Translate bases that match the reference to '='
        BASES,  ///< Translate '=' to the actual base.
    };

    /// Default Constructor.
    SamRecord();

    /// Constructor that sets the error handling type.
    /// \param errorHandlingType how to handle errors.
    SamRecord(ErrorHandler::HandlingType errorHandlingType);

    /// Destructor
    ~SamRecord();

    /// Reset the fields of the record to a default value.
    /// This is not necessary when you are reading a SAM/BAM file, 
    /// but if you are setting fields, it is a good idea to clean
    /// out a record before reusing it. Clearing it allows you to 
    /// not have to set any empty fields. 
    void resetRecord();

    /// Returns whether or not the record is valid, setting the status to
    /// indicate success or failure.
    /// \param header SAM Header associated with the record.  Used to perform
    /// some validation against the header.
    /// \return true if the record is valid, false if not.
    bool isValid(SamFileHeader& header);

    /// Set the reference to the specified genome sequence object.
    /// \param reference pointer to the GenomeSequence object.
    void setReference(GenomeSequence* reference);

    /// Set the type of sequence translation to use when getting
    /// the sequence.  The default type (if this method is never called) is
    /// NONE (the sequence is left as-is).  Can be over-ridden by using 
    /// the accessors that take a SequenceTranslation parameter.
    /// \param translation type of sequence translation to use.
    void setSequenceTranslation(SequenceTranslation translation);

    ///////////////////////
    /// @name  Set Alignment Data
    /// Set methods for record fields.  All of the "set" methods set the
    /// status to indicate success or the failure reason.
    //@{

    /// Set QNAME to the passed in name.
    /// \param readName the readname to set the QNAME to.
    /// \return true if successfully set, false if not.
    bool setReadName(const char* readName);

    /// Set the bitwise FLAG to the specified value.
    /// \param flag integer flag to use.
    /// \return true if successfully set, false if not.
    bool setFlag(uint16_t flag);
    
    /// Set the reference sequence name (RNAME) to the specified name, using
    /// the header to determine the reference id.
    /// \param header SAM/BAM header to use to determine the reference id.
    /// \param referenceName reference name to use.
    /// \return true if successfully set, false if not
    bool setReferenceName(SamFileHeader& header, 
                          const char* referenceName);

    /// Set the leftmost position (POS) using the specified 1-based (SAM format)
    /// value.
    /// Internal processing handles the switching between SAM/BAM formats 
    /// when read/written.
    /// \param position 1-based start position
    /// \return true if successfully set, false if not.
    bool set1BasedPosition(int32_t position);

    /// Set the leftmost position using the specified 0-based (BAM format)
    /// value.
    /// Internal processing handles the switching between SAM/BAM formats 
    /// when read/written.
    /// \param position 0-based start position
    /// \return true if successfully set, false if not.
    bool set0BasedPosition(int32_t position);

    /// Set the mapping quality (MAPQ).
    /// \param mapQuality map quality to set in the record.
    /// \return true if successfully set, false if not.
    bool setMapQuality(uint8_t mapQuality);

    /// Set the CIGAR to the specified SAM formatted cigar string.
    /// Internal processing handles the switching between SAM/BAM formats 
    /// when read/written.
    /// \param cigar string containing the SAM formatted cigar.
    /// \return true if successfully set, false if not.
    bool setCigar(const char* cigar);

    /// Set the CIGAR to the specified Cigar object.
    /// Internal processing handles the switching between SAM/BAM formats 
    /// when read/written.
    /// \param cigar object to set this record's cigar to have.
    /// \return true if successfully set, false if not.
    bool setCigar(const Cigar& cigar);


    /// Set the mate/next fragment's reference sequence name (RNEXT) to the
    /// specified name, using the header to determine the mate reference id.
    /// \param header SAM/BAM header to use to determine the mate reference id.
    /// \param referenceName mate reference name to use.
    /// \return true if successfully set, false if not
    bool setMateReferenceName(SamFileHeader& header,
                              const char* mateReferenceName);

    /// Set the mate/next fragment's leftmost position (PNEXT) using the
    /// specified 1-based (SAM format) value.
    /// Internal processing handles the switching between SAM/BAM formats 
    /// when read/written.
    /// \param position 1-based start position
    /// \return true if successfully set, false if not.
    bool set1BasedMatePosition(int32_t matePosition);

    /// Set the mate/next fragment's leftmost position using the specified
    /// 0-based (BAM format) value.
    /// Internal processing handles the switching between SAM/BAM formats 
    /// when read/written.
    /// \param position 0-based start position
    /// \return true if successfully set, false if not.
    bool set0BasedMatePosition(int32_t matePosition);

    /// Sets the inferred insert size (ISIZE)/observed template length (TLEN).
    /// \param insertSize inferred insert size/observed template length.
    /// \return true if successfully set, false if not.
    bool setInsertSize(int32_t insertSize);

    /// Sets the sequence (SEQ) to the specified SAM formatted sequence string.
    ///  Internal processing handles switching between SAM/BAM formats when
    /// read/written.
    /// \param seq SAM sequence string.  May contain '='.
    /// \return true if successfully set, false if not.
    bool setSequence(const char* seq);

    /// Sets the quality (QUAL) to the specified SAM formatted quality string.
    /// Internal processing handles switching between SAM/BAM formats when
    /// read/written.
    /// \param quality SAM quality string.
    /// \return true if successfully set, false if not.
    bool setQuality(const char* quality);

    /// Shift the indels (if any) to the left by updating the CIGAR.
    /// \return true if the cigar was shifted, false if not.
    bool shiftIndelsLeft();

    /// Sets the SamRecord to contain the information in the BAM formatted
    /// fromBuffer.
    /// \param fromBuffer buffer to read the BAM record from.
    /// \param fromBufferSize size of the buffer containing the BAM record.
    /// \param header BAM header for the record.
    /// \return status of reading the BAM record from the buffer.
    SamStatus::Status setBuffer(const char* fromBuffer, uint32_t fromBufferSize,
                                SamFileHeader& header);

    /// Read the BAM record from a file.
    /// \param filePtr file to read the buffer from.
    /// \param header BAM header for the record.
    /// \return status of the reading the BAM record from the file.
    SamStatus::Status setBufferFromFile(IFILE filePtr, SamFileHeader& header);

    //@}

    ///////////////////////
    /// @name  Set Tag Data
    /// Set methods for tags.
    //@{

    /// Add the specified integer tag to the record.  Internal processing
    /// handles switching between SAM/BAM formats when read/written and 
    /// determining the type for BAM format.  If the tag is already there
    /// this code will replace it if the specified value is different.
    /// \param tag two character tag to be added to the SAM/BAM record.
    /// \param value value for the specified tag.
    /// \return true if the tag was successfully added, false otherwise.
    bool addIntTag(const char* tag, int32_t value);

    /// Add the specified tag,vtype,value to the record.  Vtype can be SAM/BAM
    /// format.  Internal processing handles switching between SAM/BAM formats
    /// when read/written.  If the tag is already there this code will replace
    /// it if the specified value is different.
    /// \param tag two character tag to be added to the SAM/BAM record.
    /// \param vtype vtype of the specified value - either SAM/BAM vtypes.
    /// \param value value as a string for the specified tag.
    /// \return true if the tag was successfully added, false otherwise.
    bool addTag(const char* tag, char vtype, const char* value);

    /// Clear the tags in this record.
    /// Does not set SamStatus.
    void clearTags();
   
    /// Remove a tag.
    /// \param tag tag to remove.
    /// \param type of the tag to be removed.
    /// \return true if the tag no longer exists in the record, false if it could not be removed (Returns true if the tag was not found in the record).
    bool rmTag(const char* tag, char type);

    /// Remove tags.
    /// The delimiter between the tags is ',' or ';'.  ',' was added since 
    /// the original delimiter, ';', requires the string to be quoted on the
    /// command-line.
    /// \param tags tags to remove, formatted as  Tag:Type,Tag:Type,Tag:Type...
    /// \return true if all tags no longer exist in the record, false if any could not be removed
    /// (Returns true if the tags were not found in the record).
    /// SamStatus is set to INVALID if the tags are incorrectly formatted.
    bool rmTags(const char* tags);

    //@}

    ///////////////////////
    /// @name  Get Alignment Data
    /// Get methods for record fields.  All of the "get" methods set the
    /// status to indicate success or the failure reason.
    //@{

    /// Get a const pointer to the buffer that contains the BAM representation
    /// of the record.
    /// \return const pointer to the buffer that contains the BAM representation
    /// of the record.
    const void* getRecordBuffer();

    /// Get a const pointer to the buffer that contains the BAM representation
    /// of the record using the specified translation on the sequence.
    /// \param translation type of sequence translation to use.
    /// \return const pointer to the buffer that contains the BAM representation
    /// of the record.
    const void* getRecordBuffer(SequenceTranslation translation);

    /// Write the record as a BAM into the specified already opened file.
    /// \param filePtr file to write the BAM record into.
    /// \return status of the write.
    SamStatus::Status writeRecordBuffer(IFILE filePtr);

    /// Write the record as a BAM into the specified already opened file using
    /// the specified translation on the sequence.
    /// \param filePtr file to write the BAM record into.
    /// \param translation type of sequence translation to use.
    /// \return status of the write.
    SamStatus::Status writeRecordBuffer(IFILE filePtr, 
                                        SequenceTranslation translation);

    /// Get the block size of the record (BAM format).
    /// \return BAM block size of the record.
    int32_t getBlockSize();

    /// Get the reference sequence name (RNAME) of the record.
    /// \return reference sequence name
    const char* getReferenceName();

    /// Get the reference sequence id of the record (BAM format rid).
    /// \return reference sequence id
    int32_t getReferenceID();

    /// Get the 1-based(SAM) leftmost position (POS) of the record.
    /// \return 1-based leftmost position.
    int32_t get1BasedPosition();
 
    /// Get the 0-based(BAM) leftmost position of the record.
    /// \return 0-based leftmost position.
   int32_t get0BasedPosition();

    /// Get the length of the readname (QNAME) including the null.
    /// \return length of the read name (including null).
    uint8_t getReadNameLength();

    /// Get the mapping quality (MAPQ) of the record.
    /// \return map quality.
    uint8_t getMapQuality();

    /// Get the BAM bin for the record.
    /// \return BAM bin
    uint16_t getBin();

    /// Get the length of the BAM formatted CIGAR.
    /// \return length of BAM formatted cigar.
    uint16_t getCigarLength();

    /// Get the flag (FLAG).
    /// \return flag.
    uint16_t getFlag();

    /// Get the length of the read.
    /// \return read length.
    int32_t getReadLength();

    /// Get the mate/next fragment's reference sequence name (RNEXT).  If it
    /// is equal to the reference name, it still returns the reference name.
    /// \return reference sequence name
    const char* getMateReferenceName();

    /// Get the mate/next fragment's reference sequence name (RNEXT),
    /// returning "=" if it is the same as the reference name, unless 
    /// they are both "*" in which case "*" is returned.
    /// \return reference sequence name or '='
    const char* getMateReferenceNameOrEqual();

    /// Get the mate reference id of the record
    /// (BAM format: mate_rid/next_refID).
    /// \return reference id
    int32_t getMateReferenceID();

    /// Get the 1-based(SAM) leftmost mate/next fragment's position (PNEXT).
    /// \return 1-based leftmost position.
    int32_t get1BasedMatePosition();

    /// Get the 0-based(BAM) leftmost mate/next fragment's position.
    /// \return 0-based leftmost position.
    int32_t get0BasedMatePosition();

    /// Get the inferred insert size of the read pair (ISIZE) or
    /// observed template length (TLEN).
    /// \return inferred insert size or observed template length.
    int32_t getInsertSize();

    /// Returns the 0-based inclusive rightmost position of the
    /// clipped sequence.
    /// \return 0-based inclusive rightmost position
    int32_t get0BasedAlignmentEnd();

    /// Returns the 1-based inclusive rightmost position of the
    /// clipped sequence.
    /// \return 1-based inclusive rightmost position
    int32_t get1BasedAlignmentEnd();
   
    /// Returns the length of the clipped sequence, returning 0 if the cigar
    /// is '*'.
    /// \return length of the clipped sequence.
    int32_t getAlignmentLength();

    /// Returns the 0-based inclusive left-most position adjusted for
    /// clipped bases.
    /// \return 0-based inclusive leftmost position including clips.
    int32_t get0BasedUnclippedStart();

    /// Returns the 1-based inclusive left-most position adjusted for
    /// clipped bases.
    /// \return 1-based inclusive leftmost position including clips.
    int32_t get1BasedUnclippedStart();

    /// Returns the 0-based inclusive right-most position adjusted for
    /// clipped bases.
    /// \return 0-based inclusive rightmost position including clips.
    int32_t get0BasedUnclippedEnd();
 
    /// Returns the 1-based inclusive right-most position adjusted for
    /// clipped bases.
    /// \return 1-based inclusive rightmost position including clips.
    int32_t get1BasedUnclippedEnd();

    /// Returns the SAM formatted Read Name (QNAME).
    /// \return read name.
    const char* getReadName();

    /// Returns the SAM formatted CIGAR string.
    /// \return cigar string.
    const char* getCigar();

    /// Returns the SAM formatted sequence string (SEQ), translating the base as
    /// specified by setSequenceTranslation.
    /// \return sequence string.
    const char* getSequence();

    /// Returns the SAM formatted sequence string (SEQ) performing the specified
    /// sequence translation.
    /// \param translation type of sequence translation to use.
    /// \return sequence string.
    const char* getSequence(SequenceTranslation translation);

    /// Returns the SAM formatted quality string (QUAL).
    /// \return quality string.
    const char* getQuality();

    /// Get the sequence base at the specified index into this sequence 0 to
    /// readLength - 1, translating the base as specified by
    /// setSequenceTranslation.  Throws an exception if index is out of range.
    /// \param index index into the sequence string (0 to readLength-1).
    /// \return the sequence base at the specified index into the sequence.
    char getSequence(int index);
    
    /// Get the sequence base at the specified index into this sequence 0 to
    /// readLength - 1 performing the specified sequence translation. 
    /// Throws an exception if index is out of range.
    /// \param index index into the sequence string (0 to readLength-1).
    /// \param translation type of sequence translation to use.
    /// \return the sequence base at the specified index into the sequence.
    char getSequence(int index, SequenceTranslation translation);
    
    /// Get the quality character at the specified index into the quality 0 to
    /// readLength - 1.  Throws an exception if index is out of range.
    /// \param index index into the quality string (0 to readLength-1).
    /// \return the quality character at the specified index into the quality.
    char getQuality(int index);
   
    /// Returns a pointer to the Cigar object associated with this record.  
    /// The object is essentially read-only, only allowing modifications 
    /// due to lazy evaluations.
    /// \return pointer to the Cigar object.
    Cigar* getCigarInfo();

    /// Return the number of bases in this read that overlap the passed in
    /// region.  Matches & mismatches between the read and the reference
    /// are counted as overlaps, but insertions, deletions, skips, clips, and
    /// pads are not counted.
    /// \param start inclusive 0-based start position (reference position) of
    ///              the region to check for overlaps in.
    ///              (-1 indicates to start at the beginning of the reference.)
    /// \param end   exclusive 0-based end position (reference position) of the
    ///              region to check for overlaps in.
    ///              (-1 indicates to go to the end of the reference.)
    /// \return number of overlapping bases
    uint32_t getNumOverlaps(int32_t start, int32_t end);

    /// Returns the values of all fields except the tags.
    /// \param recStruct structure containing the contents of all 
    /// non-variable length fields.
    /// \param readName read name from the record (return param)
    /// \param cigar cigar string from the record (return param)
    /// \param sequence sequence string from the record (return param)
    /// \param quality quality string from the record (return param)
    /// \return true if all fields were successfully set, false otherwise.
    bool getFields(bamRecordStruct& recStruct, String& readName, 
                   String& cigar, String& sequence, String& quality);

    /// Returns the values of all fields except the tags using the specified
    /// sequence translation.
    /// \param recStruct structure containing the contents of all 
    /// non-variable length fields.
    /// \param readName read name from the record (return param)
    /// \param cigar cigar string from the record (return param)
    /// \param sequence sequence string from the record (return param)
    /// \param quality quality string from the record (return param)
    /// \param translation type of sequence translation to use.
    /// \return true if all fields were successfully set, false otherwise.
    bool getFields(bamRecordStruct& recStruct, String& readName, 
                   String& cigar, String& sequence, String& quality,
                   SequenceTranslation translation);

    /// Returns a pointer to the genome sequence object associated with this
    /// record if it was set (NULL if it was not set).
    /// \return pointer to the GenomeSequence object or NULL if there isn't one.
    GenomeSequence* getReference();

    //@}

    ///////////////////////
    /// @name  Get Tag Methods
    /// Get methods for obtaining information on tags.
    //@{

    /// Returns the length of the BAM formatted tags.
    /// \return length of the BAM formatted tags.
    uint32_t getTagLength();

    /// Get the next tag from the record.
    /// Sets the Status to SUCCESS when a tag is successfully returned or
    /// when there are no more tags.  Otherwise the status is set to describe
    /// why it failed (parsing, etc).
    /// \param tag set to the tag when a tag is read.
    /// \param vtype set to the vtype when a tag is read.
    /// \param value pointer to the value of the tag (will need to cast
    /// to int, float, char, or string based on vtype).
    /// \return true if a tag was read, false if there are no more tags.
    bool getNextSamTag(char* tag, char& vtype, void** value);

    /// Reset the tag iterator to the beginning of the tags.
    void resetTagIter();
 
    /// Returns whether or not the specified vtype is an integer type.
    /// Does not set SamStatus.
    /// \param vtype value type to check.
    /// \return true if the passed in vtype is an integer ('c', 'C', 's',
    /// 'S', 'i', 'I'), false otherwise.
    static bool isIntegerType(char vtype);

    /// Returns whether or not the specified vtype is a float type.
    /// Does not set SamStatus.
    /// \param vtype value type to check.
    /// \return true if the passed in vtype is a float ('f'), false otherwise.
    static bool isFloatType(char vtype);

    /// Returns whether or not the specified vtype is a char type.
    /// Does not set SamStatus.
    /// \param vtype value type to check.
    /// \return true if the passed in vtype is a char ('A'), false otherwise.
    static bool isCharType(char vtype);

    /// Returns whether or not the specified vtype is a string type.
    /// Does not set SamStatus.
    /// \param vtype value type to check.
    /// \return true if the passed in vtype is a string ('Z'/'B'), false othwerise.
    static bool isStringType(char vtype);

    /// Get the string representation of the tags from the record, formatted
    /// as TAG:TYPE:VALUE<delim>TAG:TYPE:VALUE...
    /// Sets the Status to SUCCESS when the tags are successfully returned or
    /// the tags were not found.  If a different error occured, the status is
    /// set appropriately.
    /// The delimiter between the tags to retrieve is ',' or ';'.  ',' was added
    /// since the original delimiter, ';', requires the string to be quoted on
    /// the command-line.
    /// \param tags the tags to retrieve, formatted as TAG:TYPE,TAG:TYPE...
    /// \param returnString the String to set (this method first clears returnString)
    ///                     to TAG:TYPE:VALUE<delim>TAG:TYPE:VALUE...
    /// \param delim delimiter to use to separate two tags, default is a tab.
    /// \return true if there were not any errors even if no tags were found.
    bool getTagsString(const char* tags, String& returnString, char delim = '\t');

    /// Get the string value for the specified tag.
    /// \param tag tag to retrieve
    /// \param pointer to the tag's string value if found, NULL if not found.
    const String* getStringTag(const char * tag);

    /// Get the integer value for the specified tag, DEPRECATED, use one that returns a bool (success/failure).
    /// \param tag tag to retrieve
    /// \retun pointer to the tag's integer value if found, NULL if not found.
    int* getIntegerTag(const char * tag);

    /// Get the integer value for the specified tag.
    /// \param tag tag to retrieve
    /// \param tagVal return parameter with integer value for the tag
    /// \retun bool true if Integer tag was found and tagVal was set, 
    ///             false if not.
    bool getIntegerTag(const char * tag, int& tagVal);

    /// Get the float value for the specified tag.
    /// \param tag tag to retrieve
    /// \param tagVal return parameter with integer value for the tag
    /// \return bool true if Float tag was found and tagVal was set,
    ///         false if not.
    bool getFloatTag(const char * tag, float& tagVal);

    /// Get the string value for the specified tag.
    const String & getString(const char * tag);

    /// Get the integer value for the specified tag, DEPRECATED, use getIntegerTag that returns a bool.
    int &    getInteger(const char * tag);

    /// Check if the specified tag contains a string.
    /// Does not set SamStatus.
    /// \param tag SAM tag to check contents of.
    /// \return true if the value associated with the tag is a string.
    bool checkString(const char * tag)
    { return(checkTag(tag, 'Z') || checkTag(tag, 'B')); }
    
    /// Check if the specified tag contains an integer.
    /// Does not set SamStatus.
    /// \param tag SAM tag to check contents of.
    /// \return true if the value associated with the tag is a string.
    bool checkInteger(const char * tag)   { return checkTag(tag, 'i'); }
    
    /// Check if the specified tag contains a string.
    /// Does not set SamStatus.
    /// \param tag SAM tag to check contents of.
    /// \return true if the value associated with the tag is a string.
    bool checkFloat(const char * tag)    { return checkTag(tag, 'f'); }
     
    /// Check if the specified tag contains a value of the specified vtype.
    /// Does not set SamStatus.
    /// \param tag SAM tag to check contents of.
    /// \param type value type to check if the SAM tag matches.
    /// \return true if the value associated with the tag is a string.
   bool checkTag(const char * tag, char type);
    //@}

    /// Returns the status associated with the last method that sets the status.
    /// \return SamStatus of the last command that sets status.
    const SamStatus& getStatus();


private:
    static int MAKEKEY(char ch1, char ch2, char type)
    { return (getKeyType(type) << 16) + (ch2 << 8) + ch1; }

    static char getKeyType(char type)
    {
        switch(type)
        {
            // For any char/integer type, return 'i'
            case 'A' :
            case 'c' :
            case 'C' :
            case 's' :
            case 'S' :
            case 'i' :
            case 'I' :
                return('i');
                break;
            default:
                // For all other types, return the actual type.
                return(type);
        };
    }

    static inline int getNumericTagTypeSize(char type)
    {
        switch(type)
        {
            case 'A':
            case 'c':
            case 'C':
                return(1);
                break;
            case 's':
            case 'S':
                return(2);
                break;
            case 'i':
            case 'I':
            case 'f':
                return(4);
            default:
                // Not a numeric type.
                return(0);
        }
    }

    // Allocate space for the record - does a realloc.  
    // The passed in size is the size of the entire record including the
    // block size field.
    // Adds any errors to myStatus.
    bool allocateRecordStructure(int size);

    void* getStringPtr(int offset);
    void* getIntegerPtr(int offset, char& vtype);
    void* getFloatPtr(int offset);

    // Fixes the buffer to match the variable length fields.
    // Adds any errors to myStatus.
    bool fixBuffer(SequenceTranslation translation);

    // Sets the Sequence and Quality strings from the buffer.
    // They are done together in one method because they require the same
    // loop, so might as well be done at the same time.
    // Adds any errors to myStatus.
    void setSequenceAndQualityFromBuffer();

    // Parse the cigar to calculate the alignment/unclipped ends and convert
    // to SAM/BAM format.
    // Adds any errors to myStatus.
    bool parseCigar();
    // Parse the cigar string to calculate the cigar length and alignment end
    // and convert to SAM format.
    // Adds any errors to myStatus.
    bool parseCigarBinary();
    // Parse the cigar string to calculate the cigar length and alignment end
    // and convert to BAM format.
    // Adds any errors to myStatus.
    bool parseCigarString();

    // Set the tags from the buffer.
    // Adds any errors to myStatus.
    bool setTagsFromBuffer();

    // Set the tags in the buffer.
    // Adds any errors to myStatus.
    bool setTagsInBuffer();

    void setVariablesForNewBuffer(SamFileHeader& header);

    void getTypeFromKey(int key, char& type) const;
    void getTag(int key, char* tag) const;

    String & getString(int offset);
    int &    getInteger(int offset);
    const char &   getIntegerType(int offset) const;
    float & getFloat(int offset);

    // Append the string representation of the value at the specified index
    // of the int array.
    inline void appendIntArrayValue(int index, String& strVal) const
    {
        appendIntArrayValue(intType[index], integers[index], strVal);
    }

    void appendIntArrayValue(char type, int value, String& strVal) const;

    int getBtagBufferSize(String& tagStr);
    int setBtagBuffer(String& tagStr, char* extraPtr);
    int getStringFromBtagBuffer(unsigned char* buffer, String& tagStr);

    static const int DEFAULT_BLOCK_SIZE = 40;
    static const int DEFAULT_BIN = 4680;
    static const int DEFAULT_READ_NAME_LENGTH = 8;
    static const char* DEFAULT_READ_NAME;
    static const char* FIELD_ABSENT_STRING;

    bamRecordStruct * myRecordPtr;
    int allocatedSize;

    // Pointer to a temporary cigar buffer that can be used during string
    // parsing before it is ready to be copied into the actual record.
    uint32_t* myCigarTempBuffer;

    // Size of the currently allocated temporary cigar buffer.
    int myCigarTempBufferAllocatedSize;

    // Length of the cigar currently contained in the temporary buffer.
    int myCigarTempBufferLength;

    // Track if the buffer is in sync with the Strings/Tags.
    // Set to false if any of the variable length fields are modified.
    // Set to true when the buffer is updated to match the variable length
    // fields.
    bool myIsBufferSynced;

    // Track if the tags need to be set from the buffer.
    bool myNeedToSetTagsFromBuffer;

    // Trag if the tags need to be set in the buffer.
    // Allows you to set just the tags if they are the only thing that changed
    // in the buffer.
    bool myNeedToSetTagsInBuffer;

    int myTagBufferSize;
    int myLastTagIndex;

    String myReadName;
    String myReferenceName;
    String myMateReferenceName;
    String myCigar;
    String mySequence;
    String myQuality;

    std::string mySeqWithEq;
    std::string mySeqWithoutEq;

    // The length of the alignment.
    int32_t myAlignmentLength;
    // Unclipped alignment positions.
    int32_t myUnclippedStartOffset;
    int32_t myUnclippedEndOffset;
    
    CigarRoller myCigarRoller;

    LongHash<int>  extras;
    // Note: not all values in strings, integers, and floats are always
    // in extras.  They will not be if the tags were removed.  Removed
    // tags are removed from extras, but not from strings, integers, or floats
    // since if one was removed from these arrays, all other entries would
    // need their indices updated in extras.
    StringArray    strings;
    IntArray       integers;
    std::vector<char> intType; // contains the type of int at same position in integers.
    std::vector<float> floats;


    // Track whether or not the buffer values are correct for
    // each setting.
    bool myIsReadNameBufferValid;
    bool myIsCigarBufferValid;
    bool myIsSequenceBufferValid;
    bool myIsQualityBufferValid;
    bool myIsTagsBufferValid;
    bool myIsBinValid;

    unsigned char* myPackedSequence;
    unsigned char* myPackedQuality;


    SamStatus myStatus;

    // The current translation of the sequence as it occurs in the buffer.
    // Only applicable if myIsSequenceBufferValid == true.
    SequenceTranslation myBufferSequenceTranslation;


    // Track the Reference.
    GenomeSequence* myRefPtr;

    // The type of translation to do when getting a sequence.
    SequenceTranslation mySequenceTranslation;

    String NOT_FOUND_TAG_STRING;
    int NOT_FOUND_TAG_INT;

    static const int myMaxWarns = 5;
    static int myNumWarns;
};

#endif
