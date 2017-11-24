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

#ifndef __GLF_HEADER_H__
#define __GLF_HEADER_H__

#include <stdint.h>

#include "InputFile.h" 
#include "CharBuffer.h"

/// This class allows a user to easily get/set the fields in a GLF header.
/// The GlfHeader contains:
/// - Variable length text string 
class GlfHeader
{
public:
    GlfHeader();
    ~GlfHeader();

    /// Copy Constructor   
    /// \param header glfheader to copy into this one.
    GlfHeader(const GlfHeader& header);

    /// Overload operator= to copy the passed in header into this header.
    /// \param header glfheader to copy into this one.
    GlfHeader & operator = (const GlfHeader& header);

    /// Copy the passed in header into this header.
    /// \param header glfheader to copy into this one.
    bool copy(const GlfHeader& header);

    /// Clear this header back to the default setting.
    void resetHeader();
   
    /// Read the header from the specified file (file MUST be in 
    /// the correct position for reading the header).
    /// \param filePtr file to read from that is in the correct position.
    /// \return true if the header was successfully read from the 
    /// file, false if not.
    bool read(IFILE filePtr);

    /// Write the header to the specified file.
    /// \param filePtr file to write to that is in the correct position.
    /// \return true if the header was successfully written to the 
    /// file, false if not.
    bool write(IFILE filePtr) const;

    /// Set the passed in string to the text string stored in this header.
    /// \param text string to populate with the header text string.
    /// \return true if text was successfully returned, false if not.
    bool getHeaderTextString(std::string& text);

    /// Set the header to the passed in string.
    /// \param text header text to assign to this header.
    /// \return true if the text was successfully set, false if not.
    bool setHeaderTextString(const std::string& text);

private:
    int32_t myTextLen;
    CharBuffer myText;

    static const std::string GLF_MAGIC;
    static const int GLF_MAGIC_LEN = 4;
};

#endif

