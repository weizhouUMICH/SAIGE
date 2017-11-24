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

#include "ErrorHandler.h"
#include "PhoneHome.h"

#include <stdexcept>
#include <stdlib.h>

// Constructor
ErrorHandler::ErrorHandler()
{
}

   
// Destructor
ErrorHandler::~ErrorHandler()
{
}


void ErrorHandler::handleError(const char* message, 
                               HandlingType handlingType)
{
    // Check the handling type.
    switch(handlingType)
    {
        case(EXCEPTION):
            throw(std::runtime_error(message));
            break;
        case(ABORT):
            std::cerr << message << "\nExiting" << std::endl;
            PhoneHome::completionStatus("ErrorHandler: Exiting due to Error"); 
            exit(-1);
            break;
        case(RETURN):
            return;
            break;
        default:
            std::cerr << message << "\nUnknown Handle Type: Exiting" 
                      << std::endl;
            PhoneHome::completionStatus("Exiting, ErrorHandler::unknown handle type."); 
            exit(-1);
            break;
    }
}

