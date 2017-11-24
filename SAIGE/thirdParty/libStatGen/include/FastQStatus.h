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

#ifndef __FASTQ_STATUS_H__
#define __FASTQ_STATUS_H__

#include <string>

/// Status for FastQ operations.
class FastQStatus
{
 public:

    /// Return value enum for the FastQFile class methods, indicating
    /// success or error codes.
   enum Status 
       {
           FASTQ_SUCCESS = 0,      ///< indicates method finished successfully.
           FASTQ_INVALID,          ///< means that the sequence was invalid.
           FASTQ_ORDER_ERROR,      ///< means the methods are called out of order, like trying to read a file before opening it.
           FASTQ_OPEN_ERROR,       ///< means the file could not be opened.
           FASTQ_CLOSE_ERROR,      ///< means the file could not be closed.
           FASTQ_READ_ERROR,       ///< means that a problem occurred on a read.
           FASTQ_NO_SEQUENCE_ERROR ///< means there were no errors, but no sequences read.
       };

   /// Get the enum string for the status.
   static const char* getStatusString(Status status);

private:
   static const char* enumString[];
};


#endif
