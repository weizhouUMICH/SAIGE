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

#include "SamFile.h"
#include "WriteFiles.h"
#include "TestValidate.h"

#include <assert.h>
#include <stdio.h>

void testWrite()
{
    testHeaderWrite();
    testWriteCopiedHeader("testFiles/testSam.sam");
#ifdef __ZLIB_AVAILABLE__
    testWriteCopiedHeader("testFiles/testBam.bam");
#endif
}

void testHeaderWrite()
{
    SamFile samOut;

    samOut.OpenForWrite("results/MyTestOut.sam");

    // Create a sam header.
    SamFileHeader samHeader;
   
    std::string headerString = "";

    // Test getting HD & PG and the HD-SO tag when they do not exist.
    assert(samHeader.getHD() == NULL);
    assert(samHeader.getPG("1") == NULL);
    assert(strcmp(samHeader.getTagSO(), "") == 0);

    // Test removing the HD tag that does not exist.
    assert(samHeader.removeHD() == true);
    assert(samHeader.getHD() == NULL);
    assert(strcmp(samHeader.getHDTagValue("VN"), "") == 0);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "");

    char type[3] = "HD";
    char tag[3] = "VN";
    // Verify it has not yet been added to the parsed header.
    assert(strcmp(samHeader.getHDTagValue("VN"), "") == 0);
    assert(samHeader.addHeaderLine(type, tag, "1.0") == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tVN:1.0\n");

    // Verify it was added to the parsed header.
    assert(strcmp(samHeader.getHDTagValue("VN"), "1.0") == 0);

    type[0] = 'S';
    type[1] = 'Q';
    tag[0] = 'L';
    tag[1] = 'N';
    // Cannot add SQ LN tag without adding the SN tag also.
    assert(samHeader.addHeaderLine(type, tag, "123") == false);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tVN:1.0\n");

    // Has not yet been added, so returns blank.
    assert(strcmp(samHeader.getSQTagValue("LN", "chr20"), "") == 0);

    // Can't add the SQ type without a LN.
    std::string line = "@SQ\tSN:123";
    assert(samHeader.addHeaderLine(line.c_str()) == false);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tVN:1.0\n");

    // Successfully add a SQ line.
    line = "@SQ\tLN:123\tSN:chr20";
    assert(samHeader.addHeaderLine(line.c_str()) == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tVN:1.0\n@SQ\tLN:123\tSN:chr20\n");
    // Verify it was added to the parsed header.
    assert(strcmp(samHeader.getSQTagValue("LN", "chr20"), "123") == 0);

    // Test to make sure nothing changes if try to copy into self.
    samHeader = samHeader;
    assert(samHeader.addHeaderLine(line.c_str()) == false);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tVN:1.0\n@SQ\tLN:123\tSN:chr20\n");
    // Verify it was added to the parsed header.
    assert(strcmp(samHeader.getSQTagValue("LN", "chr20"), "123") == 0);

    samHeader.copy(samHeader);
    assert(samHeader.addHeaderLine(line.c_str()) == false);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tVN:1.0\n@SQ\tLN:123\tSN:chr20\n");
    // Verify it was added to the parsed header.
    assert(strcmp(samHeader.getSQTagValue("LN", "chr20"), "123") == 0);

    // Test adding an HD that is already there.
    assert(samHeader.addHeaderLine("@HD\tVN:1.1") == false);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tVN:1.0\n@SQ\tLN:123\tSN:chr20\n");
    // Verify it was added to the parsed header.
    assert(strcmp(samHeader.getHDTagValue("VN"), "1.0") == 0);

    // Test copying the header.
    SamFileHeader newHeader = samHeader;
    assert(newHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tVN:1.0\n@SQ\tLN:123\tSN:chr20\n");
    // Verify it was added to the parsed header.
    assert(strcmp(newHeader.getSQTagValue("LN", "chr20"), "123") == 0);
    

    // Modify one of the tags.
    assert(samHeader.setHDTag("VN", "1.1") == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tVN:1.1\n@SQ\tLN:123\tSN:chr20\n");
    // Verify it was modified.
    assert(strcmp(samHeader.getHDTagValue("VN"), "1.1") == 0);

    // Remove the version.
    assert(samHeader.setHDTag("VN", "") == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@SQ\tLN:123\tSN:chr20\n");
    // Verify it was removed.
    assert(strcmp(samHeader.getHDTagValue("VN"), "") == 0);
   
    // Remove the SN from the SQ type - fails because SN is the key.
    assert(samHeader.setSQTag("SN", "", "chr20") == false);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@SQ\tLN:123\tSN:chr20\n");
    // Verify it was not removed.
    assert(strcmp(samHeader.getSQTagValue("SN", "chr20"), "chr20") == 0);
    assert(strcmp(samHeader.getSQTagValue("LN", "chr20"), "123") == 0);
   
    // Can't remove the LN from the SQ type
    assert(samHeader.setSQTag("LN", "", "chr20") == false);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@SQ\tLN:123\tSN:chr20\n");
    // Verify it was not removed.
    assert(strcmp(samHeader.getSQTagValue("LN", "chr20"), "123") == 0);
    assert(strcmp(samHeader.getSQTagValue("SN", "chr20"), "chr20") == 0);

    // Delete the SQ line.
    assert(samHeader.removeSQ("chr20") == true);
    // There is no header string.
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "");
    // Verify it was removed.
    assert(strcmp(samHeader.getSQTagValue("SN", "chr20"), "") == 0);
    assert(strcmp(samHeader.getSQTagValue("LN", "chr20"), "") == 0);

    // Create an SQ record and add it back in.
    SamHeaderSQ* sq = new SamHeaderSQ();
    assert(sq->setTag("LN", "123") == true);
    assert(sq->setTag("SN", "chr20") == true);
    assert(samHeader.addSQ(sq) == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@SQ\tLN:123\tSN:chr20\n");
    // Verify it was added.
    assert(strcmp(samHeader.getSQTagValue("SN", "chr20"), "chr20") == 0);
    assert(strcmp(samHeader.getSQTagValue("LN", "chr20"), "123") == 0);

    // Modify a tag.
    assert(sq->setTag("LN", "222") == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@SQ\tLN:222\tSN:chr20\n");
    // Verify it was modified.
    assert(strcmp(samHeader.getSQTagValue("SN", "chr20"), "chr20") == 0);
    assert(strcmp(samHeader.getSQTagValue("LN", "chr20"), "222") == 0);
   

    // Test adding another SQ with the same key.
    SamHeaderSQ* sq2 = new SamHeaderSQ();
    assert(sq2->setTag("LN", "333") == true);
    assert(sq2->setTag("SN", "chr20") == true);
    assert(samHeader.addSQ(sq2) == false);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@SQ\tLN:222\tSN:chr20\n");
    // Verify it was not added.
    assert(strcmp(samHeader.getSQTagValue("SN", "chr20"), "chr20") == 0);
    assert(strcmp(samHeader.getSQTagValue("LN", "chr20"), "222") == 0);
    delete sq2;

    // Add a new tag to the SQ tag.
    assert(samHeader.setSQTag("AS", "HG18", "chr20") == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@SQ\tLN:222\tSN:chr20\tAS:HG18\n");
    // Verify it was added.
    assert(strcmp(samHeader.getSQTagValue("AS", "chr20"), "HG18") == 0);
   
    // Modify the AS tag.
    assert(samHeader.setSQTag("AS", "HG19", "chr20") == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@SQ\tLN:222\tSN:chr20\tAS:HG19\n");
    // Verify it was added.
    assert(strcmp(samHeader.getSQTagValue("AS", "chr20"), "HG19") == 0);

    // Add a new tag .
    sq2 = new SamHeaderSQ();
    assert(sq2->setTag("LN", "333") == true);
    assert(sq2->setTag("SN", "chr1") == true);
    assert(samHeader.addSQ(sq2) == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@SQ\tLN:222\tSN:chr20\tAS:HG19\n@SQ\tLN:333\tSN:chr1\n");
    // Verify it was added.
    assert(strcmp(samHeader.getSQTagValue("SN", "chr1"), "chr1") == 0);
    assert(strcmp(samHeader.getSQTagValue("LN", "chr1"), "333") == 0);

    // Test removing an SQ tag that does not exist.
    assert(samHeader.removeSQ("chr100") == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@SQ\tLN:222\tSN:chr20\tAS:HG19\n@SQ\tLN:333\tSN:chr1\n");

    // Remove the newly added sq2 by resetting it.
    sq2->reset();
    // Verify it was removed.
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@SQ\tLN:222\tSN:chr20\tAS:HG19\n");
    assert(strcmp(samHeader.getSQTagValue("SN", "chr1"), "") == 0);
    assert(strcmp(samHeader.getSQTagValue("LN", "chr1"), "") == 0);

    // Test getting HD which does exist since it was set before.  Even
    // though the VN was removed so it doesn't appear in the header string,
    // it was never actually removed.
    SamHeaderHD* hd = samHeader.getHD();
    assert(hd != NULL);
    // Blank since the sort order was never set.
    assert(strcmp(samHeader.getTagSO(), "") == 0);

    // Set the version number.
    assert(hd->setTag("VN", "2.1") == true);
    // Verify it was added.
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tVN:2.1\n@SQ\tLN:222\tSN:chr20\tAS:HG19\n");
    assert(strcmp(samHeader.getHDTagValue("VN"), "2.1") == 0);

    // Set the SO
    assert(hd->setTag("SO", "coordinate") == true);
    // Verify it was added.
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@HD\tVN:2.1\tSO:coordinate\n@SQ\tLN:222\tSN:chr20\tAS:HG19\n");
    assert(strcmp(samHeader.getHDTagValue("SO"), "coordinate") == 0);

    // Reset the header.
    samHeader.resetHeader();
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "");
   
    // Add a new HD tag.
    assert(samHeader.setHDTag("SO", "queryname") == true);
    assert(strcmp(samHeader.getHDTagValue("SO"), "queryname") == 0);
    assert(samHeader.getHeaderString(headerString) == true);
    // Blank since missing VN.
    assert(headerString == "");

    // Set the VN.
    assert(samHeader.setHDTag("VN", "3.1") == true);
    assert(strcmp(samHeader.getHDTagValue("SO"), "queryname") == 0);
    assert(strcmp(samHeader.getHDTagValue("VN"), "3.1") == 0);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tSO:queryname\tVN:3.1\n");
   
    //////////////////////////////////////////////////////////////
    // Test removing a non-existent PG.
    assert(samHeader.removePG("1") == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tSO:queryname\tVN:3.1\n");

    // Test adding a null PG.
    SamHeaderPG* pg = NULL;
    assert(samHeader.addPG(pg) == false);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tSO:queryname\tVN:3.1\n");
   
    // Add a PG tag.
    assert(strcmp(samHeader.getPGTagValue("ID", "pid"), "") == 0);
    assert(samHeader.setPGTag("ID", "pid", "pid") == true);
    assert(strcmp(samHeader.getPGTagValue("ID", "pid"), "pid") == 0);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tSO:queryname\tVN:3.1\n@PG\tID:pid\n");

    // Verify can't modify the key.
    assert(samHeader.setPGTag("ID", "pid1", "pid") == false);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tSO:queryname\tVN:3.1\n@PG\tID:pid\n");
    assert(strcmp(samHeader.getPGTagValue("ID", "pid"), "pid") == 0);

    // Test adding a null PG.
    pg = NULL;
    assert(samHeader.addPG(pg) == false);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tSO:queryname\tVN:3.1\n@PG\tID:pid\n");
   
    // Test adding a PG header when it already exists.
    pg = new SamHeaderPG();
    assert(pg->setTag("ID", "pid") == true);
    assert(samHeader.addPG(pg) == false);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tSO:queryname\tVN:3.1\n@PG\tID:pid\n");
    delete pg;

    // Get a PG that does not exist.
    pg = samHeader.getPG("pid1");
    assert(pg == NULL);

    // Get a PG tag that does not exist.
    assert(strcmp(samHeader.getPGTagValue("CL", "pid"), "") == 0);

    // Get the PG tag.
    pg = samHeader.getPG("pid");
    assert(pg != NULL);
    assert(strcmp(pg->getTagValue("ID"), "pid") == 0);
    // Add a tag to the PG.
    assert(pg->setTag("VN", "pg1") == true);
    assert(strcmp(samHeader.getPGTagValue("VN", "pid"), "pg1") == 0);
    assert(strcmp(samHeader.getPGTagValue("ID", "pid"), "pid") == 0);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tSO:queryname\tVN:3.1\n@PG\tID:pid\tVN:pg1\n");
   
    // Test modifying the key tag - fails.
    assert(pg->setTag("ID", "pid1") == false);
    assert(strcmp(samHeader.getPGTagValue("VN", "pid"), "pg1") == 0);
    assert(strcmp(samHeader.getPGTagValue("ID", "pid"), "pid") == 0);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tSO:queryname\tVN:3.1\n@PG\tID:pid\tVN:pg1\n");

    // Test modifying the VN tag.
    assert(samHeader.setPGTag("VN", "pg", "pid") == true);
    assert(strcmp(samHeader.getPGTagValue("VN", "pid"), "pg") == 0);
    assert(strcmp(samHeader.getPGTagValue("ID", "pid"), "pid") == 0);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tSO:queryname\tVN:3.1\n@PG\tID:pid\tVN:pg\n");
   
    // Test removing the VN tag.
    assert(pg->setTag("VN", "") == true);
    assert(strcmp(samHeader.getPGTagValue("VN", "pid"), "") == 0);
    assert(strcmp(samHeader.getPGTagValue("ID", "pid"), "pid") == 0);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tSO:queryname\tVN:3.1\n@PG\tID:pid\n");
   
    // Test removing a PG that does not exist.
    assert(samHeader.removePG("pid1") == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tSO:queryname\tVN:3.1\n@PG\tID:pid\n");
    assert(strcmp(samHeader.getPGTagValue("VN", "pid"), "") == 0);
    assert(strcmp(samHeader.getPGTagValue("ID", "pid"), "pid") == 0);

    // Test removing the PG.
    assert(samHeader.removePG("pid") == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tSO:queryname\tVN:3.1\n");
    assert(strcmp(samHeader.getPGTagValue("VN", "pid"), "") == 0);
    assert(strcmp(samHeader.getPGTagValue("ID", "pid"), "") == 0);

    // Test adding a PG header.
    pg = new SamHeaderPG();
    assert(pg->setTag("ID", "newID") == true);
    assert(samHeader.addPG(pg) == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n");

    // Test adding a PG that is already there.
    assert(samHeader.addHeaderLine("@PG\tID:newID") == false);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n");
    // Verify it was added to the parsed header.
    assert(strcmp(samHeader.getPGTagValue("ID", "newID"), "newID") == 0);

    // Test adding another PG header.
    pg = new SamHeaderPG();
    assert(pg->setTag("ID", "newID1") == true);
    assert(samHeader.addPG(pg) == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@PG\tID:newID1\n");

    // Test adding another PG header.
    pg = new SamHeaderPG();
    assert(pg->setTag("ID", "pid") == true);
    assert(samHeader.addPG(pg) == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@PG\tID:newID1\n@PG\tID:pid\n");

    // Test removing the new pg.
    assert(samHeader.removePG("newID1") == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@PG\tID:pid\n");

    // Test removing the other new pg.
    assert(samHeader.removePG("pid") == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n");

    // Test adding a tag
    assert(samHeader.setPGTag("VN", "1.0", "newID") == true);
    assert(strcmp(samHeader.getPGTagValue("VN", "newID"), "1.0") == 0);
    assert(strcmp(samHeader.getPGTagValue("ID", "newID"), "newID") == 0);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\tVN:1.0\n");
   
    // Test removing a tag
    assert(samHeader.setPGTag("VN", "", "newID") == true);
    assert(strcmp(samHeader.getPGTagValue("VN", "newID"), "") == 0);
    assert(strcmp(samHeader.getPGTagValue("ID", "newID"), "newID") == 0);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n");

    ////////////////////////////////////////////////////////////////////
    // Add an SQ, but fail since LN is not specified.
    assert(samHeader.setSQTag("AS", "HG18", "newName") == false);
    assert(samHeader.getHeaderString(headerString) == true);
    // SQ does not show up since it is missing the LN field.
    assert(headerString == 
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n");
    // Add the SQ's SN, but fail since LN is not specified.
    assert(samHeader.setSQTag("SN", "newName", "newName") == false);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n");
     sq = samHeader.getSQ("newName");
    assert(sq == NULL);
   // Add the SQ with the LN tag.
    assert(samHeader.setSQTag("LN", "111", "newName") == true);
    assert(strcmp(samHeader.getSQTagValue("SN", "newName"), "newName") == 0);
    assert(strcmp(samHeader.getSQTagValue("AS", "newName"), "") == 0);
    assert(strcmp(samHeader.getSQTagValue("LN", "newName"), "111") == 0);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\n");
    // Add the AS.
    assert(samHeader.setSQTag("AS", "HG18", "newName") == true);
    assert(strcmp(samHeader.getSQTagValue("AS", "newName"), "HG18") == 0);
    assert(strcmp(samHeader.getSQTagValue("SN", "newName"), "newName") == 0);
    assert(strcmp(samHeader.getSQTagValue("LN", "newName"), "111") == 0);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\n");

    // Get the SQ.
    sq = samHeader.getSQ("newName");
    assert(sq != NULL);
    // Modify the SQ
    assert(sq->setTag("SP", "species") == true);
    assert(strcmp(samHeader.getSQTagValue("SN", "newName"), "newName") == 0);
    assert(strcmp(samHeader.getSQTagValue("AS", "newName"), "HG18") == 0);
    assert(strcmp(samHeader.getSQTagValue("LN", "newName"), "111") == 0);
    assert(strcmp(samHeader.getSQTagValue("SP", "newName"), "species") == 0);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n");

    //////////////////////////////////////////////////////////////////////
    // Add a new RG Tag
// Remove this test since SM is no longer required, the RG would be added
//    assert(samHeader.setRGTag("ID", "rgID", "rgID") == true);
//    assert(samHeader.getHeaderString(headerString) == true);
//    // New RG does not show up since it is still missing a required field.
//    assert(headerString == 
//           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n");
//    assert(strcmp(samHeader.getRGTagValue("ID", "rgID"), "rgID") == 0);
   
    // Add the missing SM field.
    assert(samHeader.setRGTag("SM", "sm1", "rgID") == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n");
    assert(strcmp(samHeader.getRGTagValue("ID", "rgID"), "rgID") == 0);
   
    // Verify can't modify the key.
    assert(samHeader.setRGTag("ID", "rgID1", "rgID") == false);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n");
    assert(strcmp(samHeader.getRGTagValue("ID", "rgID"), "rgID") == 0);

    // Verify that the copied header did not change.
    assert(newHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tVN:1.0\n@SQ\tLN:123\tSN:chr20\n");
    // Verify it was added to the parsed header.
    assert(strcmp(newHeader.getSQTagValue("LN", "chr20"), "123") == 0);

    // Add a new RG Tag
    assert(samHeader.setRGTag("SM", "sample1", "rgID1") == true);
    assert(samHeader.getHeaderString(headerString) == true);
    // String does not show the tag until all required fields are there.
    assert(headerString ==
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n@RG\tID:rgID1\tSM:sample1\n");
    assert(strcmp(samHeader.getRGTagValue("ID", "rgID1"), "rgID1") == 0);
    assert(strcmp(samHeader.getRGTagValue("SM", "rgID1"), "sample1") == 0);
   
    // Modify an RG tag.
    assert(samHeader.setRGTag("SM", "sample", "rgID1") == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n@RG\tID:rgID1\tSM:sample\n");
    assert(strcmp(samHeader.getRGTagValue("ID", "rgID1"), "rgID1") == 0);
    assert(strcmp(samHeader.getRGTagValue("SM", "rgID1"), "sample") == 0);

    // Test removing an rg that does not exist.
    assert(samHeader.removeRG("rgID2") == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n@RG\tID:rgID1\tSM:sample\n");
   
    // Create a new RG, set some values and add it.
    SamHeaderRG* rg = new SamHeaderRG();
    // Try adding it without a key.
    assert(samHeader.addRG(rg) == false);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n@RG\tID:rgID1\tSM:sample\n");
    // Set some values in rg.
    assert(rg->setTag("ID", "rgID2") == true);
    assert(rg->setTag("SM", "sm2") == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n@RG\tID:rgID1\tSM:sample\n");
    // Add the new RG.
    assert(samHeader.addRG(rg) == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n@RG\tID:rgID1\tSM:sample\n@RG\tID:rgID2\tSM:sm2\n");
    assert(strcmp(samHeader.getRGTagValue("ID", "rgID2"), "rgID2") == 0);

    // Test trying to add another one with the same key.
    rg = new SamHeaderRG();
    assert(rg->setTag("ID", "rgID2") == true);
    assert(samHeader.addRG(rg) == false);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n@RG\tID:rgID1\tSM:sample\n@RG\tID:rgID2\tSM:sm2\n");


    // Test removing the rg again.
    assert(samHeader.removeRG("rgID2") == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString ==
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n@RG\tID:rgID1\tSM:sample\n");
   
    // Test getting an rg tag that doesn't exist
    assert(strcmp(samHeader.getRGTagValue("DS", "rgID"), "") == 0);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n@RG\tID:rgID1\tSM:sample\n");
   
    // Test getting an rg tag from a removed key
    assert(strcmp(samHeader.getRGTagValue("ID", "rgID2"), "") == 0);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n@RG\tID:rgID1\tSM:sample\n");

    // Test getting an rg tag from an key that doesn't exist
    assert(strcmp(samHeader.getRGTagValue("ID", "rgID22"), "") == 0);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n@RG\tID:rgID1\tSM:sample\n");

    // Test adding a null header.
    rg = NULL;
    assert(samHeader.addRG(rg) == false);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n@RG\tID:rgID1\tSM:sample\n");
   
    // Test adding the deleted header back in.
    rg = new SamHeaderRG();
    assert(rg->setTag("ID", "rgID2") == true);
    assert(rg->setTag("SM", "sm2") == true);
    assert(samHeader.addRG(rg) == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString ==
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n@RG\tID:rgID1\tSM:sample\n@RG\tID:rgID2\tSM:sm2\n");

    // Test adding an RG that is already there.
    assert(samHeader.addHeaderLine("@RG\tID:rgID\tSM:sm5") == false);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n@RG\tID:rgID1\tSM:sample\n@RG\tID:rgID2\tSM:sm2\n");
    // Verify it was added to the parsed header.
    assert(strcmp(samHeader.getRGTagValue("SM", "rgID"), "sm1") == 0);


    // Get an RG record then modify it.
    rg = samHeader.getRG("rgID1");
    assert(rg != NULL);
    assert(rg->setTag("SM", "sample1") == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == 
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n@RG\tID:rgID1\tSM:sample1\n@RG\tID:rgID2\tSM:sm2\n");
   
    // Try to modify the key.
    assert(rg->setTag("ID", "rgID111") == false);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString ==
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n@RG\tID:rgID1\tSM:sample1\n@RG\tID:rgID2\tSM:sm2\n");
   
    ////////////////////////////////////////////////////////////////////////////
    // Test getting a comment when there aren't any.
    assert(strcmp(samHeader.getNextComment(), "") == 0);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString ==
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n@RG\tID:rgID1\tSM:sample1\n@RG\tID:rgID2\tSM:sm2\n");

    // Test getting each headerline when there are no comments.
    const char* hdrlinechar;
    std::string hdrline;
    assert(samHeader.getNextHeaderLine(hdrline));
    hdrlinechar = hdrline.c_str();
    // Test to make sure there is not memory corruption.
    std::string tmpString = "@SQ\tSN:queryname\tVN:3.1\n";
    assert(hdrline == "@HD\tSO:queryname\tVN:3.1\n");
    assert(strcmp(hdrlinechar,
                  "@HD\tSO:queryname\tVN:3.1\n") == 0);

    assert(samHeader.getNextHeaderLine(hdrline));
    assert(hdrline == "@PG\tID:newID\n");
    assert(samHeader.getNextHeaderLine(hdrline));
    assert(hdrline == "@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n");
    assert(samHeader.getNextHeaderLine(hdrline)); 
    assert(hdrline == "@RG\tID:rgID\tSM:sm1\n");
    assert(samHeader.getNextHeaderLine(hdrline));
    assert(hdrline == "@RG\tID:rgID1\tSM:sample1\n");
    assert(samHeader.getNextHeaderLine(hdrline)); 
    assert(hdrline == "@RG\tID:rgID2\tSM:sm2\n");
    assert(!samHeader.getNextHeaderLine(hdrline));
    assert(hdrline == "");
    assert(!samHeader.getNextHeaderLine(hdrline));
    assert(hdrline == "");
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString ==
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n@RG\tID:rgID1\tSM:sample1\n@RG\tID:rgID2\tSM:sm2\n");

    // Verify that getHeaderRecord returns nothing.
    assert(samHeader.getNextHeaderRecord() == NULL);

    // Reset the header record iter.
    samHeader.resetHeaderRecordIter();

    // Test getting each headerrecord when there are no comments.
    SamHeaderRecord* hdrRec = samHeader.getNextHeaderRecord();
    assert(hdrRec != NULL);
    assert(strcmp(hdrRec->getTypeString(), "HD") == 0);
    assert(hdrRec->getType() == SamHeaderRecord::HD);
    assert(strcmp(hdrRec->getTagValue("SO"), "queryname") == 0);
    assert(strcmp(hdrRec->getTagValue("VN"), "3.1") == 0);
    hdrRec = samHeader.getNextHeaderRecord();
    assert(hdrRec != NULL);
    assert(strcmp(hdrRec->getTypeString(), "PG") == 0);
    assert(hdrRec->getType() == SamHeaderRecord::PG);
    assert(strcmp(hdrRec->getTagValue("ID"), "newID") == 0);
    hdrRec = samHeader.getNextHeaderRecord();
    assert(hdrRec != NULL);
    assert(strcmp(hdrRec->getTypeString(), "SQ") == 0);
    assert(hdrRec->getType() == SamHeaderRecord::SQ);
    assert(strcmp(hdrRec->getTagValue("SN"), "newName") == 0);
    assert(strcmp(hdrRec->getTagValue("AS"), "HG18") == 0);
    assert(strcmp(hdrRec->getTagValue("LN"), "111") == 0);
    hdrRec = samHeader.getNextHeaderRecord();
    assert(hdrRec != NULL);
    assert(strcmp(hdrRec->getTypeString(), "RG") == 0);
    assert(hdrRec->getType() == SamHeaderRecord::RG);
    assert(strcmp(hdrRec->getTagValue("ID"), "rgID") == 0);
    assert(strcmp(hdrRec->getTagValue("SM"), "sm1") == 0);

    // Get the SQ Header Record (should have no affect on the general
    // getNextHeaderRecord calls).
    hdrRec = samHeader.getNextSQRecord();
    assert(hdrRec != NULL);
    assert(strcmp(hdrRec->getTypeString(), "SQ") == 0);
    assert(hdrRec->getType() == SamHeaderRecord::SQ);
    assert(strcmp(hdrRec->getTagValue("SN"), "newName") == 0);
    assert(strcmp(hdrRec->getTagValue("AS"), "HG18") == 0);
    assert(strcmp(hdrRec->getTagValue("LN"), "111") == 0);
    // Only one SQ Header Record.
    hdrRec = samHeader.getNextSQRecord();
    assert(hdrRec == NULL);

    // Get the RG/PG Header Records (should have no affect on the general
    // getNextHeaderRecord calls).
    hdrRec = samHeader.getNextRGRecord();
    assert(hdrRec != NULL);
    assert(strcmp(hdrRec->getTypeString(), "RG") == 0);
    assert(hdrRec->getType() == SamHeaderRecord::RG);
    assert(strcmp(hdrRec->getTagValue("ID"), "rgID") == 0);
    assert(strcmp(hdrRec->getTagValue("SM"), "sm1") == 0);
    // Get the next RG record.
    hdrRec = samHeader.getNextRGRecord();
    assert(hdrRec != NULL);
    assert(strcmp(hdrRec->getTypeString(), "RG") == 0);
    assert(hdrRec->getType() == SamHeaderRecord::RG);
    assert(strcmp(hdrRec->getTagValue("ID"), "rgID1") == 0);
    assert(strcmp(hdrRec->getTagValue("SM"), "sample1") == 0);
    // Get the PG record.
    hdrRec = samHeader.getNextPGRecord();
    assert(hdrRec != NULL);
    assert(strcmp(hdrRec->getTypeString(), "PG") == 0);
    assert(hdrRec->getType() == SamHeaderRecord::PG);
    assert(strcmp(hdrRec->getTagValue("ID"), "newID") == 0);
    // Get the last RG record.
    hdrRec = samHeader.getNextRGRecord();
    assert(hdrRec != NULL);
    assert(strcmp(hdrRec->getTypeString(), "RG") == 0);
    assert(hdrRec->getType() == SamHeaderRecord::RG);
    assert(strcmp(hdrRec->getTagValue("ID"), "rgID2") == 0);
    assert(strcmp(hdrRec->getTagValue("SM"), "sm2") == 0);
    // Already got all RG Records.
    hdrRec = samHeader.getNextRGRecord();
    assert(hdrRec == NULL);
    // Reset the RG record.
    samHeader.resetRGRecordIter();
    // Get the RG record.
    hdrRec = samHeader.getNextRGRecord();
    assert(hdrRec != NULL);
    assert(strcmp(hdrRec->getTypeString(), "RG") == 0);
    assert(hdrRec->getType() == SamHeaderRecord::RG);
    assert(strcmp(hdrRec->getTagValue("ID"), "rgID") == 0);
    assert(strcmp(hdrRec->getTagValue("SM"), "sm1") == 0);
    // No more PG records.
    hdrRec = samHeader.getNextPGRecord();
    assert(hdrRec == NULL);
    // No more SQ records.
    hdrRec = samHeader.getNextSQRecord();
    assert(hdrRec == NULL);    
    // Reset the SQ record iterator.
    samHeader.resetSQRecordIter();
    // No more PG records.
    hdrRec = samHeader.getNextPGRecord();
    assert(hdrRec == NULL);
    // Get the now reset SQ record.
    hdrRec = samHeader.getNextSQRecord();
    assert(hdrRec != NULL);
    assert(strcmp(hdrRec->getTypeString(), "SQ") == 0);
    assert(hdrRec->getType() == SamHeaderRecord::SQ);
    assert(strcmp(hdrRec->getTagValue("SN"), "newName") == 0);
    assert(strcmp(hdrRec->getTagValue("AS"), "HG18") == 0);
    assert(strcmp(hdrRec->getTagValue("LN"), "111") == 0);
    // Only one SQ Header Record.
    hdrRec = samHeader.getNextSQRecord();
    assert(hdrRec == NULL);
    // Reset the PG record iterator.
    samHeader.resetPGRecordIter();
    // No more SQ records.
    hdrRec = samHeader.getNextSQRecord();
    assert(hdrRec == NULL);    
    // Get the next RG record.
    hdrRec = samHeader.getNextRGRecord();
    assert(hdrRec != NULL);
    assert(strcmp(hdrRec->getTypeString(), "RG") == 0);
    assert(hdrRec->getType() == SamHeaderRecord::RG);
    assert(strcmp(hdrRec->getTagValue("ID"), "rgID1") == 0);
    assert(strcmp(hdrRec->getTagValue("SM"), "sample1") == 0);
    // Get the PG record.
    hdrRec = samHeader.getNextPGRecord();
    assert(hdrRec != NULL);
    assert(strcmp(hdrRec->getTypeString(), "PG") == 0);
    assert(hdrRec->getType() == SamHeaderRecord::PG);
    assert(strcmp(hdrRec->getTagValue("ID"), "newID") == 0);


    hdrRec = samHeader.getNextHeaderRecord();
    assert(hdrRec != NULL);
    assert(strcmp(hdrRec->getTypeString(), "RG") == 0);
    assert(hdrRec->getType() == SamHeaderRecord::RG);
    assert(strcmp(hdrRec->getTagValue("ID"), "rgID1") == 0);
    assert(strcmp(hdrRec->getTagValue("SM"), "sample1") == 0);
    hdrRec = samHeader.getNextHeaderRecord();
    assert(hdrRec != NULL);
    assert(strcmp(hdrRec->getTypeString(), "RG") == 0);
    assert(hdrRec->getType() == SamHeaderRecord::RG);
    assert(strcmp(hdrRec->getTagValue("ID"), "rgID2") == 0);
    assert(strcmp(hdrRec->getTagValue("SM"), "sm2") == 0);
    hdrRec = samHeader.getNextHeaderRecord();
    assert(hdrRec == NULL);
    hdrRec = samHeader.getNextHeaderRecord();
    assert(hdrRec == NULL);
    assert(!samHeader.getNextHeaderLine(hdrline));
    assert(hdrline == "");
    assert(!samHeader.getNextHeaderLine(hdrline));
    assert(hdrline == "");

    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString ==
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n@RG\tID:rgID1\tSM:sample1\n@RG\tID:rgID2\tSM:sm2\n");

    // Add some comments.
    assert(samHeader.addComment("My Comment") == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString ==
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n@RG\tID:rgID1\tSM:sample1\n@RG\tID:rgID2\tSM:sm2\n@CO\tMy Comment\n");

    // Call getNextHeaderRecord - still nothing.
    hdrRec = samHeader.getNextHeaderRecord();
    assert(hdrRec == NULL);   

    // Call getNextHeaderLine - should return the comment.
    assert(samHeader.getNextHeaderLine(hdrline));
    assert(hdrline == "@CO\tMy Comment\n");
    assert(!samHeader.getNextHeaderLine(hdrline));
    assert(hdrline == "");
    assert(!samHeader.getNextHeaderLine(hdrline));
    assert(hdrline == "");

    // Call getNextCommentLine - should return the comment.
    assert(strcmp(samHeader.getNextComment(), "My Comment") == 0);
    assert(strcmp(samHeader.getNextComment(), "") == 0);
    assert(strcmp(samHeader.getNextComment(), "") == 0);

    // Add another comment.
    assert(samHeader.addComment("My Comment2") == true);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString ==
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n@RG\tID:rgID1\tSM:sample1\n@RG\tID:rgID2\tSM:sm2\n@CO\tMy Comment\n@CO\tMy Comment2\n");

    newHeader = samHeader;
    assert(newHeader.getHeaderString(headerString) == true);
    assert(headerString ==
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n@RG\tID:rgID1\tSM:sample1\n@RG\tID:rgID2\tSM:sm2\n@CO\tMy Comment\n@CO\tMy Comment2\n");

    // Call getNextHeaderLine - should return the comment.
    assert(samHeader.getNextHeaderLine(hdrline));
    assert(hdrline == "@CO\tMy Comment2\n");
    assert(!samHeader.getNextHeaderLine(hdrline));
    assert(hdrline == "");
    assert(!samHeader.getNextHeaderLine(hdrline));
    assert(hdrline == "");

    // Call getNextCommentLine - should return the comment.
    assert(strcmp(samHeader.getNextComment(), "My Comment2") == 0);
    assert(strcmp(samHeader.getNextComment(), "") == 0);
    assert(strcmp(samHeader.getNextComment(), "") == 0);

    // Reset the header record iter.
    samHeader.resetHeaderRecordIter();
   
    // Recall getNextCommentLine - should not return anything.
    assert(strcmp(samHeader.getNextComment(), "") == 0);
    assert(strcmp(samHeader.getNextComment(), "") == 0);   

    // Reset the next comment iter.
    samHeader.resetCommentIter();

    // Call the get next headerLine, record, comment interspersed with
    // each other.
    hdrRec = samHeader.getNextHeaderRecord();
    assert(hdrRec != NULL);
    assert(strcmp(hdrRec->getTypeString(), "HD") == 0);
    assert(hdrRec->getType() == SamHeaderRecord::HD);
    assert(strcmp(hdrRec->getTagValue("SO"), "queryname") == 0);
    assert(strcmp(hdrRec->getTagValue("VN"), "3.1") == 0);
    assert(samHeader.getNextHeaderLine(hdrline));
    assert(hdrline == "@PG\tID:newID\n");
    assert(samHeader.getNextHeaderLine(hdrline)); 
    assert(hdrline == "@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n");   
    hdrRec = samHeader.getNextHeaderRecord();
    assert(hdrRec != NULL);
    assert(strcmp(hdrRec->getTypeString(), "RG") == 0);
    assert(hdrRec->getType() == SamHeaderRecord::RG);
    assert(strcmp(hdrRec->getTagValue("ID"), "rgID") == 0);
    hdrRec = samHeader.getNextHeaderRecord();
    assert(hdrRec != NULL);
    assert(strcmp(samHeader.getNextComment(), "My Comment") == 0);
    assert(strcmp(hdrRec->getTypeString(), "RG") == 0);
    assert(hdrRec->getType() == SamHeaderRecord::RG);
    assert(strcmp(hdrRec->getTagValue("ID"), "rgID1") == 0);
    assert(strcmp(hdrRec->getTagValue("SM"), "sample1") == 0);
    assert(samHeader.getNextHeaderLine(hdrline));
    assert(hdrline == "@RG\tID:rgID2\tSM:sm2\n");
    hdrRec = samHeader.getNextHeaderRecord();
    assert(hdrRec == NULL);
    assert(samHeader.getNextHeaderLine(hdrline));
    assert(hdrline == "@CO\tMy Comment\n");
    hdrRec = samHeader.getNextHeaderRecord();
    assert(hdrRec == NULL);
    assert(samHeader.getNextHeaderLine(hdrline));
    assert(hdrline == "@CO\tMy Comment2\n");
    assert(!samHeader.getNextHeaderLine(hdrline));
    assert(hdrline == "");
    assert(strcmp(samHeader.getNextComment(), "My Comment2") == 0);
    assert(!samHeader.getNextHeaderLine(hdrline));
    assert(hdrline == "");
    hdrRec = samHeader.getNextHeaderRecord();
    assert(hdrRec == NULL);
    assert(strcmp(samHeader.getNextComment(), "") == 0);
    assert(strcmp(samHeader.getNextComment(), "") == 0);

    samOut.WriteHeader(samHeader);

    // Reset the header.
    samHeader.resetHeader();
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "");
    assert(!samHeader.getNextHeaderLine(hdrline)); 
    assert(hdrline == "");
    assert(strcmp(samHeader.getHDTagValue("SO"), "") == 0);
    assert(strcmp(samHeader.getHDTagValue("VN"), "") == 0);

    // Try adding a key to the HD tag.
    hd = new SamHeaderHD();
    assert(hd->addKey("3.1") == false);
    assert(strcmp(hd->getTagValue("VN"), "") == 0);
    assert(hd->isActiveHeaderRecord() == false);

    assert(hd->setTag("VN", "3.1") == true);
    assert(hd->isActiveHeaderRecord() == true);
    assert(strcmp(hd->getTagValue("VN"), "3.1") == 0);

    // Verify the copied header did not change.
    assert(newHeader.getHeaderString(headerString) == true);
    assert(headerString ==
           "@HD\tSO:queryname\tVN:3.1\n@PG\tID:newID\n@SQ\tSN:newName\tLN:111\tAS:HG18\tSP:species\n@RG\tID:rgID\tSM:sm1\n@RG\tID:rgID1\tSM:sample1\n@RG\tID:rgID2\tSM:sm2\n@CO\tMy Comment\n@CO\tMy Comment2\n");
    // Verify it was added to the parsed header.
    assert(strcmp(newHeader.getSQTagValue("LN", "chr20"), "") == 0);
}


void testWriteCopiedHeader(const char* fileName)
{
    SamFile samIn;
    assert(samIn.OpenForRead(fileName));

    SamFile samOut;
    assert(samOut.OpenForWrite("results/MyTestOut2.bam"));
    SamFile samOut2;
    assert(samOut2.OpenForWrite("results/MyTestOut2.sam"));

    // Read the sam header.
    SamFileHeader samHeader;
    assert(samIn.ReadHeader(samHeader));
    validateHeader(samHeader);

    SamFileHeader newHeader;

    std::string hdrLine;
    assert(samHeader.getNextHeaderLine(hdrLine));
    newHeader.addHeaderLine("@HD\tVN:1.02");
    bool hdrStatus = true;
    while(hdrStatus)
    {
        newHeader.addHeaderLine(hdrLine.c_str());
        hdrStatus = samHeader.getNextHeaderLine(hdrLine);
    }

    // Write the sam header.
    assert(samOut.WriteHeader(newHeader));
    assert(samOut2.WriteHeader(newHeader));

    SamRecord samRecord;

    // Keep reading records until ReadRecord returns false.
    while(samIn.ReadRecord(samHeader, samRecord))
    {
        // Successfully read a record from the file, so write it.
        assert(samOut.WriteRecord(newHeader, samRecord));
        assert(samOut2.WriteRecord(newHeader, samRecord));
    }

    assert(samIn.GetStatus() == SamStatus::NO_MORE_RECS);

    // Close the output files.
    samOut.Close();
    samOut2.Close();

    SamFileReader bamRead("results/MyTestOut2.bam");
    SamFileReader samRead("results/MyTestOut2.sam");
    
    // Read and check the header.
    assert(samRead.ReadHeader(samHeader));
    validateHeaderFields(samHeader);
    std::string headerString = "";
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tVN:1.02\n@SQ\tSN:1\tLN:247249719\n@SQ\tSN:2\tLN:242951149\n@SQ\tSN:3\tLN:199501827\n@SQ\tSN:4\tLN:191273063\n@SQ\tSN:5\tLN:180857866\n@SQ\tSN:6\tLN:170899992\n@SQ\tSN:7\tLN:158821424\n@SQ\tSN:8\tLN:146274826\n@SQ\tSN:9\tLN:140273252\n@SQ\tSN:10\tLN:135374737\n@SQ\tSN:11\tLN:134452384\n@SQ\tSN:12\tLN:132349534\n@SQ\tSN:13\tLN:114142980\n@SQ\tSN:14\tLN:106368585\n@SQ\tSN:15\tLN:100338915\n@SQ\tSN:16\tLN:88827254\n@SQ\tSN:17\tLN:78774742\n@SQ\tSN:18\tLN:76117153\n@SQ\tSN:19\tLN:63811651\n@SQ\tSN:20\tLN:62435964\n@SQ\tSN:21\tLN:46944323\n@SQ\tSN:22\tLN:49691432\n@SQ\tSN:X\tLN:154913754\n@RG\tID:myID\tLB:library\tSM:sample\n@RG\tID:myID2\tSM:sample2\tLB:library2\n@CO\tComment 1\n@CO\tComment 2\n");
    
    assert(bamRead.ReadHeader(samHeader));
    validateHeaderFields(samHeader);
    headerString = "";
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@HD\tVN:1.02\n@SQ\tSN:1\tLN:247249719\n@SQ\tSN:2\tLN:242951149\n@SQ\tSN:3\tLN:199501827\n@SQ\tSN:4\tLN:191273063\n@SQ\tSN:5\tLN:180857866\n@SQ\tSN:6\tLN:170899992\n@SQ\tSN:7\tLN:158821424\n@SQ\tSN:8\tLN:146274826\n@SQ\tSN:9\tLN:140273252\n@SQ\tSN:10\tLN:135374737\n@SQ\tSN:11\tLN:134452384\n@SQ\tSN:12\tLN:132349534\n@SQ\tSN:13\tLN:114142980\n@SQ\tSN:14\tLN:106368585\n@SQ\tSN:15\tLN:100338915\n@SQ\tSN:16\tLN:88827254\n@SQ\tSN:17\tLN:78774742\n@SQ\tSN:18\tLN:76117153\n@SQ\tSN:19\tLN:63811651\n@SQ\tSN:20\tLN:62435964\n@SQ\tSN:21\tLN:46944323\n@SQ\tSN:22\tLN:49691432\n@SQ\tSN:X\tLN:154913754\n@RG\tID:myID\tLB:library\tSM:sample\n@RG\tID:myID2\tSM:sample2\tLB:library2\n@CO\tComment 1\n@CO\tComment 2\n");
    
    assert(samHeader.getNextHeaderLine(hdrLine));
    std::string expectedString = "@HD\tVN:1.02\n";

    assert(expectedString == hdrLine);
    

    // TODO - validate reading these written records back in.

}

