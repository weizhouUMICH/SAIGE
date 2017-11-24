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

#include "FastQStatus.h"

const char* FastQStatus::enumString[] = {"FASTQ_SUCCESS", "FASTQ_INVALID", "FASTQ_ORDER_ERROR", "FASTQ_OPEN_ERROR", "FASTQ_CLOSE_ERROR", "FASTQ_READ_ERROR", "FASTQ_NO_SEQUENCE_ERROR"};


const char* FastQStatus::getStatusString(Status status)
{
   return(enumString[status]);
}
