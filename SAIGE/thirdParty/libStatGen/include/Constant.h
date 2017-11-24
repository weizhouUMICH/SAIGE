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

#ifndef _CONSTANT_H_
#define _CONSTANT_H_

#define  COMPAREFUNC (int (*)(const void *, const void *))

#define  BUFSIZE     1024
#define  FILENAMELEN 100
#define  IDLEN       20

#define  SEPARATORS  " \t\n\r\f/"
#define  WHITESPACE  " \t\n\r\f"

#define  SWTABLESKIP 9
#define  SWTABLEMAX  10000

#define  _NAN_       ((double) (6.66666e-66))

#define  QTDTDATA    "qtdt.dat"
#define  QTDTPED     "qtdt.ped"
#define  QTDTIBD     "qtdt.ibd"
#define  QTDTRAW     "regress.tbl"
#define  GENIHDATAIN "genih.dat"

#ifndef  _WIN32
#define  stricmp     strcasecmp
#endif

// Constants for older haplotype handling programs
// Constants for HAPLOXT
#define XT_MAX_ALLELES  50          // Maximum alleles for crosstabulation
#define XT_VECTORSIZE   10000       // Total haplotypes in population
#define XT_POOLTRESH    7           // Threshold for pooling rare alleles
// Simwalk Haplotype Vectors
#define HV_MAXSIZE      100         // Haplotypes in single SimWalk pedigree
#define HV_INFOTRESH    75          // Percentage of loci typed
#define HV_STATELENGTH  100         // Markers per haplotype
#define HV_SKIPLINES    4           // lines to skip at bottom of family tree
// Simwalk Summary Files
#define HT_TABLE_SIZE   1000
#define HT_SKIP_LINES   9

#endif

