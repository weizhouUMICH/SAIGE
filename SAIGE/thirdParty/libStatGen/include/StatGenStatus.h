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

#ifndef __STATGEN_STATUS_H__
#define __STATGEN_STATUS_H__

#include <iostream>
#include "ErrorHandler.h"

/// This class is used to track the status results of some methods in the BAM
/// classes. It contains a status enum that describing the status.
class StatGenStatus
{
public:

    /// Return value enum for StatGenFile methods.
    enum Status 
        { SUCCESS = 0, ///< method completed successfully.
          UNKNOWN, ///< unknown result (default value should never be used)
          /// NO_MORE_RECS: failed to read a record since there are no more to
          /// read either in the file or section if section based reading.
          NO_MORE_RECS,
          FAIL_IO, ///< method failed due to an I/O issue.
          /// FAIL_ORDER: method failed because it was called out of order,
          /// like trying to read a file without opening it for read or trying
          /// to read a record before the header.
          FAIL_ORDER,
          FAIL_PARSE, ///< failed to parse a record/header - invalid format.
          INVALID_SORT, ///< record is invalid due to it not being sorted.
          INVALID, ///< invalid other than for sorting.
          FAIL_MEM ///< fail a memory allocation.
        };

    /// Return a string representation of the passed in status enum.
    static const char* getStatusString(StatGenStatus::Status statusEnum);

    /// Returns whether or not it is "safe" to keep processing the file
    /// after the specified status return.
    static bool isContinuableStatus(StatGenStatus::Status status);

    /// Constructor that takes in the handling type, defaulting it to exception.
    StatGenStatus(ErrorHandler::HandlingType handleType = ErrorHandler::EXCEPTION);
   
    /// Destructor
    ~StatGenStatus();

    /// Reset this status to a default state.
    void reset();

    /// Set how to handle the errors when they are set.
    void setHandlingType(ErrorHandler::HandlingType handleType);

    /// Set the status with the specified status enum and message.
    void setStatus(Status newStatus, const char* newMessage);

    /// Add the specified error message to the status message, setting
    /// the status to newStatus if the current status is SUCCESS.
    void addError(Status newStatus, const char* newMessage);


    /// Add the specified status to the status message, setting
    /// the status to newStatus if the current status is SUCCESS.
    void addError(StatGenStatus newStatus);

    /// Return the enum for this status object.
    Status getStatus() const;

    /// Return the status message for this object.
    const char* getStatusMessage() const;

    /// Overload operator = to set the StatGen status type to the
    /// passed in status and to clear the message string.
    StatGenStatus & operator = (Status newStatus);

    /// Overload operator = to copy the specified status object to this one.
    StatGenStatus & operator = (StatGenStatus newStatus);

    /// Overload operator != to determine if the passed in type is not equal
    /// to this status's type.
    bool operator != (const StatGenStatus::Status& compStatus) const;

    /// Overload operator == to determine if the passed in type is equal
    /// to this status's type.
    bool operator == (const StatGenStatus::Status& compStatus) const;
      
private:
    // Handle an error based on the error handling type.
    void handleError(Status newType, const char* newMessage);


    static const char* enumStatusString[];

    Status myType;
    std::string myMessage;
    ErrorHandler::HandlingType myHandlingType;
};


#endif
