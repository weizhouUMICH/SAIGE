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

#ifndef __GLF_EXCEPTION_H__
#define __GLF_EXCEPTION_H__

#include <stdexcept> // stdexcept header file

#include "GlfStatus.h"

/// GlfException objects should be thrown by functions that operate on 
/// Glf files for exceptions.
class GlfException : public std::exception
{
public:
    /// Constructor that sets the exception to a default status
    /// and error message.
    GlfException();

    /// Constructor that sets the exception to a default status
    /// and the specified error message. 
    /// \param what_arg error message associated with this exception.
    GlfException(const std::string& what_arg);

    /// Constructor that sets the exception to the specified status
    /// and error message. 
    /// \param status glf status associated with this exception.
    /// \param errorMsg error message associated with this exception.
    GlfException(GlfStatus::Status status, const std::string& errorMsg);

    /// Constructor that sets the exception to the specified status. 
    /// \param status glf status associated with this exception.
    GlfException(const GlfStatus& status);

    virtual ~GlfException() throw();

    /// Returns the error message of this exception. 
    /// \return errror message
    virtual const char* what() const throw();

private:
    GlfStatus myStatus;
}; // end class GlfException


#endif
