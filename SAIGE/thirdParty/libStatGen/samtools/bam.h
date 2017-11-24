/* The MIT License

   Copyright (c) 2008-2010 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#ifndef BAM_BAM_H
#define BAM_BAM_H

/*!
  @header

  BAM library provides I/O and various operations on manipulating files
  in the BAM (Binary Alignment/Mapping) or SAM (Sequence Alignment/Map)
  format. It now supports importing from or exporting to SAM, sorting,
  merging, generating pileup, and quickly retrieval of reads overlapped
  with a specified region.

  @copyright Genome Research Ltd.
 */


#include <stdint.h>

/*
 * Only small section is pulled out that we use.
 * 4/29/2011 - Mary Kate Trost
 */


/*!
  @abstract    Calculate the minimum bin that contains a region [beg,end).
  @param  beg  start of the region, 0-based
  @param  end  end of the region, 0-based
  @return      bin
 */
static inline int bam_reg2bin(uint32_t beg, uint32_t end)
{
	--end;
	if (beg>>14 == end>>14) return 4681 + (beg>>14);
	if (beg>>17 == end>>17) return  585 + (beg>>17);
	if (beg>>20 == end>>20) return   73 + (beg>>20);
	if (beg>>23 == end>>23) return    9 + (beg>>23);
	if (beg>>26 == end>>26) return    1 + (beg>>26);
	return 0;
}

#endif
