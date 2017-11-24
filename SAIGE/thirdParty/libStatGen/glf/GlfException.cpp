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

#include "GlfException.h"

GlfException::GlfException()
    : myStatus()
{
    myStatus.setStatus(GlfStatus::UNKNOWN, "Failed operating on a GLF.");
}


GlfException::GlfException(const std::string& errorMsg)
    : myStatus()
{
    myStatus.setStatus(GlfStatus::UNKNOWN, errorMsg.c_str());
}
 
GlfException::GlfException(GlfStatus::Status status,
                           const std::string& errorMsg)
    : myStatus()
{
    myStatus.setStatus(status, errorMsg.c_str());
}
 
GlfException::GlfException(const GlfStatus& status)
    : myStatus()
{
    myStatus.addError(status);
}
 
GlfException::~GlfException() throw()
{
}

const char* GlfException::what() const throw()
{
    return(myStatus.getStatusMessage());
}
