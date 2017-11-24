/*
 *  Copyright (C) 2011-2012  Regents of the University of Michigan,
 *                           Hyun Min Kang, Matthew Flickenger, Matthew Snyder,
 *                           and Goncalo Abecasis
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
#include "VcfRecordInfo.h"

#include <string>

VcfRecordInfo::VcfRecordInfo()
{
    reset();
}


VcfRecordInfo::~VcfRecordInfo()
{
}


bool VcfRecordInfo::read(IFILE filePtr)
{
    // Clear out any previously set values.
    reset();
    
    if(ifeof(filePtr))
    {
        // End of file, just return false.
        return(false);
    }
    
    static const std::string keyStopChars   = "\n\t;=";
    static const std::string valueStopChars = "\n\t;";
    // The start of the first character in stopChars
    // that means there is more information for this object, so
    // continue reading the file.
    static const int contPos = 2;
    static const int tabPos = 1;

    // Keep reading.  Loop will be exited
    // when a \t, \n, or EOF is found.
    int stopPos = contPos;
    while(stopPos >= contPos)
    {
        // Get the next element to write the key into.
        InfoElement& nextElement = myInfo.getNextEmpty();
        // Read the next key.
        stopPos = filePtr->readTilChar(keyStopChars, nextElement.key);

        if(keyStopChars[stopPos] == '=')
        {
            // Stoped at the value part, so read the value
            // associated with the key.
            stopPos = filePtr->readTilChar(valueStopChars, nextElement.value);
        }
    }

    // Return whether or not a tab was found at the end of the field.
    return(stopPos == tabPos);
}


bool VcfRecordInfo::write(IFILE filePtr)
{
    // If there are no entries, write '.'.
    int infoSize = myInfo.size();
    if(infoSize <= 0)
    {
        return(ifprintf(filePtr, "%c", EMPTY_INFO) == 1);
    }

    int numWritten = 0;
    int numExpected = 0;

    for(int i = 0; i < infoSize; i++)
    {
        if(i != 0)
        {
            numWritten += ifprintf(filePtr, ";");
            ++numExpected;
        }
        InfoElement& info = myInfo.get(i);
        if(info.value.empty())
        {
            // No value, just a key.
            numWritten += ifprintf(filePtr, "%s", info.key.c_str());
            numExpected += info.key.size();
        }
        else
        {
            // write the key & the value.
            numWritten += ifprintf(filePtr, "%s=%s", info.key.c_str(),
                                  info.value.c_str());
            numExpected += info.key.size() + info.value.size() + 1;
        }
    }
    return(numWritten == numExpected);
}


void VcfRecordInfo::reset()
{
    myInfo.reset();
}


void VcfRecordInfo::setString(const char* key, const char* stringVal)
{
    // Check if the field is already there.
    int infoSize = myInfo.size();
    for(int i = 0; i < infoSize; i++)
    {
        InfoElement& info = myInfo.get(i);
        if(info.key.compare(key) == 0)
        {
            // Set the value and return.
            info.value = stringVal;
            return;
        }
    }

    // Not found, so add a new entry.
    InfoElement& newElement = myInfo.getNextEmpty();
    newElement.key = key;
    newElement.value = stringVal;
}


const std::string* VcfRecordInfo::getString(const char* key)
{
    // Check if the field is already there.
    int infoSize = myInfo.size();
    for(int i = 0; i < infoSize; i++)
    {
        InfoElement& info = myInfo.get(i);
        if(info.key.compare(key) == 0)
        {
            // Found, so return the value.
            return(&(info.value));
        }
    }

    // Not found, so return NULL..
    return(NULL);
}


const std::string* VcfRecordInfo::getString(int index)
{
    if(index >= myInfo.size())
    {
        // Out of range.
        return(NULL);
    }

    return(&(myInfo.get(index).value));
}

std::pair<std::string, std::string> VcfRecordInfo::getInfoPair(int index) const
{
    if (index < myInfo.size())
    {
        InfoElement& e = myInfo.get(index);
        return std::pair<std::string, std::string>(e.key, e.value);
    }

    return std::pair<std::string, std::string>();
}