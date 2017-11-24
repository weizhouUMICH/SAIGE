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

#include "GlfStatus.h"

const char* GlfStatus::enumStatusString[] = {
    "SUCCESS",
    "UNKNOWN",
    "FAIL_IO",
    "FAIL_ORDER",
    "FAIL_PARSE",
    "INVALID", 
    "FAIL_MEM"
};


const char* GlfStatus::getStatusString(GlfStatus::Status statusEnum)
{
    return(enumStatusString[statusEnum]);
}


// Returns whether or not it is "safe" to keep processing the file
// after the specified status return.
bool GlfStatus::isContinuableStatus(GlfStatus::Status status)
{
    if(status == GlfStatus::SUCCESS || status == GlfStatus::FAIL_PARSE || 
       status == GlfStatus::INVALID)
    {
        // The status is such that file processing can continue.
        return(true);
    }
    // UNKNOWN, FAIL_IO, FAIL_ORDER, FAIL_MEM
    return(false);
}


// Constructor
GlfStatus::GlfStatus()
{
    reset();
}

   
// Destructor
GlfStatus::~GlfStatus()
{
}


// Resets this status.
void GlfStatus::reset()
{
    myType = UNKNOWN;
    myMessage = "";
}


// Set the status with the specified values.
void GlfStatus::setStatus(Status newStatus, const char* newMessage)
{
    myType = newStatus;
    myMessage = getStatusString(newStatus);
    myMessage += ": ";
    myMessage += newMessage;
}


// Adds the specified error message to the status message.
// Sets the status to newStatus if the current status is SUCCESS.
void GlfStatus::addError(Status newStatus, const char* newMessage)
{
    if(myType == GlfStatus::SUCCESS)
    {
        myType = newStatus;
    }
    else
    {
        myMessage += "\n";
    }
    myMessage += getStatusString(newStatus);
    myMessage += ": ";
    myMessage += newMessage;
}


// Adds the specified status to the status message.
// Sets the status to newStatus if the current status is SUCCESS.
void GlfStatus::addError(GlfStatus newStatus)
{
    if(myType == GlfStatus::SUCCESS)
    {
        myType = newStatus.myType;
    }
    else
    {
        myMessage += "\n";
    }
    myMessage += newStatus.myMessage;
}


// Return the enum for this status.
GlfStatus::Status GlfStatus::getStatus() const
{
    return(myType);
}


// Return the status message.
const char* GlfStatus::getStatusMessage() const
{
    return(myMessage.c_str());
}


// Overload operator = to set the glf status type to the
// passed in status and to clear the message string.
GlfStatus & GlfStatus::operator = (GlfStatus::Status newStatus)
{
    reset();
    myType = newStatus;
    return(*this);
}

// Overload operator != to determine if the passed in type is not equal
// to this status's type.
bool GlfStatus::operator != (const GlfStatus::Status& compStatus) const
{
    return(compStatus != myType);
}
// Overload operator != to determine if the passed in type is equal
// to this status's type.
bool GlfStatus::operator == (const GlfStatus::Status& compStatus) const
{
    return(compStatus == myType);
}
