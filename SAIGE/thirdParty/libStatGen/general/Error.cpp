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

#include "Error.h"

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include "PhoneHome.h"

// Declare a dummy class to ensure that compilers recognize this as C++ code
class String;

void error(const char * msg, ...)
{
    va_list  ap;

    va_start(ap, msg);

    printf("\nFATAL ERROR - \n");
    vprintf(msg, ap);
    printf("\n\n");

    va_end(ap);

    PhoneHome::completionStatus("error: Exiting due to Fatal Error"); 
    exit(EXIT_FAILURE);
}

void warning(const char * msg, ...)
{
    va_list  ap;

    va_start(ap, msg);

    fprintf(stderr,"\n\aWARNING - \n");
    vfprintf(stderr,msg, ap);
    fprintf(stderr,"\n");

    va_end(ap);
}

void numerror(const char * msg , ...)
{
    va_list  ap;

    va_start(ap, msg);

    printf("\nFATAL NUMERIC ERROR - ");
    vprintf(msg, ap);
    printf("\n\n");

    va_end(ap);

    exit(EXIT_FAILURE);
}
