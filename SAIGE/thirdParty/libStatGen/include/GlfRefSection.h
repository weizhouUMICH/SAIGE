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

#ifndef __GLF_REFSECTION_H__
#define __GLF_REFSECTION_H__

#include <stdint.h>

#include "InputFile.h" 
#include "CharBuffer.h"

///This class allows a user to easily get/set the fields in a
/// GLF section/chromosome header.
/// The GlfRefSection contains:
/// - Reference Sequence Name
/// - Reference Sequence Length 
class GlfRefSection
{
public:
    GlfRefSection();
    ~GlfRefSection();

    /// Copy Constructor
    /// \param refSection reference section to copy into this one.
    GlfRefSection(const GlfRefSection& refSection);

    /// Overload operator= to copy the passed in refSection into this one.
    /// \param refSection reference section to copy into this one.
    GlfRefSection & operator = (const GlfRefSection& refSection);

    /// Copy the passed in refSection into this refSection.
    /// \param refSection reference section to copy into this one.
    bool copy(const GlfRefSection& refSection);

    /// Clear this reference section back to the default setting.
    void resetRefSection();
   
    /// Read the refSection from the specified file (file MUST be in
    /// the correct position for reading a refSection).
    /// \param filePtr file to read from that is in the correct position.
    /// \return true if the reference section was successfully read from the 
    /// file, false if not.
    bool read(IFILE filePtr);

    /// Write the refSection to the specified file.
    /// \param filePtr file to write to that is in the correct position.
    /// \return true if the reference section was successfully written to the 
    /// file, false if not.
    bool write(IFILE filePtr) const;

    /////////////
    // Accessors.

    /// Get the reference name.
    /// \param name string to populate with the reference name.
    /// \return true if the name was successfully returned, false if not.
    bool getName(std::string& name) const;

    /// Get the length of the reference sequence.
    /// \return reference sequence length for this reference section.
    uint32_t getRefLen() const;

    /// Set the reference name.
    /// \param name reference name to set this section to.
    /// \return true if the name was successfully set, false if not.
    bool setName(const std::string& name);
    /// Set the length of the reference sequence.
    /// \param refLen reference sequence length to set this section to.
    /// \return true if the length was successfully set, false if not.
    bool setRefLen(uint32_t refLen);

    /// Print the reference section in a readable format.
    void print() const;

private:
    CharBuffer myRefName;
    uint32_t myRefLen;
};

#endif

