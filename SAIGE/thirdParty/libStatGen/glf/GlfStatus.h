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

#ifndef __GLF_STATUS_H__
#define __GLF_STATUS_H__

#include <iostream>

/// This class is used to track the status results of some methods in the
/// GLF classes using the status enum that is defined in this class to 
/// describe the return value of a method. 
class GlfStatus
{
public:

    /// Return value enum for the GlfFile class methods.
    enum Status {
        SUCCESS = 0, ///< method completed successfully.
        UNKNOWN,     ///< unknown result (default value should never be used)
        FAIL_IO,     ///< method failed due to an I/O issue.
        FAIL_ORDER,  ///< method failed because it was called out of order,
                     ///< like trying to read a file without opening it for
                     ///< read or trying to read a record before the header.
        FAIL_PARSE,  ///< failed to parse a record/header - invalid format.
        INVALID,     ///< invalid.
        FAIL_MEM     ///< fail a memory allocation.
    };

    /// Returns the string representation of the specified enum.
    /// \param statusEnum enum to convert to a string
    /// \return string representation of the enum
    static const char* getStatusString(GlfStatus::Status statusEnum);

    /// Returns whether or not it is "safe" to keep processing the file
    /// after the specified status return.
    /// \param status enum to check if it is "safe" to continue processing.
    /// \return whether or not it is "safe" to keep processing the file
    ///         after receiving the specified enum.
    static bool isContinuableStatus(GlfStatus::Status status);

    /// Constructor
    GlfStatus();
   
    /// Destructor
    ~GlfStatus();

    /// Resets this status.
    void reset();

    /// Set the status with the specified values.
    /// \param newStatus new status to set this object to.
    /// \param newMessage message associated with the new status
    void setStatus(Status newStatus, const char* newMessage);

    /// Adds the specified error message to the status message, setting
    /// the status to newStatus if the current status is SUCCESS.
    /// \param newStatus status to add to this object.
    /// \param newMessage message to add to this object
    void addError(Status newStatus, const char* newMessage);


    /// Adds the specified status to the status message, setting
    /// the status to newStatus if the current status is SUCCESS.
    /// \param newStatus status to add to this object.
    void addError(GlfStatus newStatus);

    /// Return the enum for this status.
    /// \return enum for this status object.
    Status getStatus() const;

    /// Return the status message.
    /// \return status message associate with this status object.
    const char* getStatusMessage() const;

    /// Overload operator = to set the glf status type to the
    /// passed in status and to clear the message string.
    /// \param newStatus new status to set this object to.
    /// \return this object.
    GlfStatus & operator = (Status newStatus);
   
    // Overload operator = to set the glf status.
    //    GlfStatus & operator = (GlfStatus newStatus);
   
    /// Overload operator != to determine if the passed in type is not equal
    /// to this status's type.
    /// \param compStatus status enum to compare this status object to.
    /// \return true if they are not equal, false if they are.
    bool operator != (const GlfStatus::Status& compStatus) const;

    /// Overload operator != to determine if the passed in type is equal
    /// to this status's type.
    /// \param compStatus status enum to compare this status object to.
    /// \return true if they are equal, false if they are not.
    bool operator == (const GlfStatus::Status& compStatus) const;
      
private:
    static const char* enumStatusString[];

    Status myType;
    std::string myMessage;
};

#endif
