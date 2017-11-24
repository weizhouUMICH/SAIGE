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

#ifndef __SAMHEADER_TAG_H__
#define __SAMHEADER_TAG_H__

#include <string>

class SamHeaderTag
{
public:
    SamHeaderTag(const char* tag, const char* value);
    SamHeaderTag(const SamHeaderTag&);

    ~SamHeaderTag();

    // Add this tag to the passed in tag string.
    // If the tag value is blank, the tag will not be added to the
    // passed in string.
    // NOTE: does not clear tagString.
    bool getTagString(std::string& tagString);

    // Set this tag to the passed in tag and value.
    bool setTag(const char* tag, const char* value);

    // Set the value associated with this tag to the passed in value.
    bool setValue(const char* value);

    // Return the tag for this tag.
    const char* getTag();

    // Return the value associated with this tag.
    const char* getValue();

    // Return true if there is a non-blank value associated with this tag.
    bool hasValue();

private:
    SamHeaderTag();
    SamHeaderTag& operator=(const SamHeaderTag&);

    std::string myTag;
    std::string myValue;
};

#endif
