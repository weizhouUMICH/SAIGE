/*
 *  Copyright (C) 2013  Regents of the University of Michigan
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

#include "PhoneHome.h"
#include "knetfile.h"

#include <time.h>
#include <iostream>
#include <string.h>

int PhoneHome::allThinning = 50;
int PhoneHome::ourNumber = -1;

bool PhoneHome::ourEnableCompletionStatus = false;
std::string PhoneHome::ourBaseURL = "http://csgph.sph.umich.edu/ph/";
std::string PhoneHome::ourURL = ourBaseURL;
char PhoneHome::ourPrefixChar = '?';
String PhoneHome::ourReturnString = "";
String PhoneHome::ourToolName = "";

void PhoneHome::enableCompletionStatus(const char* programName)
{
    if(programName != NULL)
    {
        add("pgm", programName);
    }
    ourEnableCompletionStatus = true;
}


void PhoneHome::disableCompletionStatus()
{
    ourEnableCompletionStatus = false;
}


bool PhoneHome::checkVersion(const char* programName, const char* version,
                             const char* params)
{
    enableCompletionStatus();
    add("pgm", programName);
    add("vsn", version);
    add("args", params);

    connect();

    // Look for this program in the returned string.
    int start = ourReturnString.Find(ourToolName+"\t");
    if(start < 0)
    {
        // Parse the toolName, and check for the program name
        // just up to a ':'
        int colStart = ourToolName.FastFindChar(':');
        if(colStart >= 0)
        {
            ourToolName.SetLength(colStart);
            start = ourReturnString.Find(ourToolName+"\t");
        }
    }

    if(start < 0)
    {
        // This program name was not found in the version file,
        // so it is a program for which version is not tracked,
        // just return true.
        return(true);
    }

    // Found this program, so extract the version.
    start += ourToolName.Length();
    while((start < ourReturnString.Length()) && 
           isspace(ourReturnString[start]))
    {
        // Consume whitespace
        ++start;
    }

    // Start now contains the position of the start of the version
    String thisVersion = version;
    String latestVersion;
    int end = start;
    while((end < ourReturnString.Length()) && 
          !isspace(ourReturnString[end]))
    {
        latestVersion += ourReturnString[end];
        ++end;
    }

    //    std::cerr << "latest version = " << latestVersion << "\nthis version = " << thisVersion.c_str() << "\n";

    if(latestVersion.FastCompare(thisVersion) > 0)
    {
        std::cerr << "\n**************************************************************************************\n"
                  << "A new version, " << latestVersion 
                  << ", of " << ourToolName
                  << " is available (currently running " 
                  << thisVersion.c_str() << ")"
                  << "\n**************************************************************************************\n\n";
        return(false);
    }
    return(true);
}

void PhoneHome::completionStatus(const char* status, const char* programName)
{
    if(programName != NULL)
    {
        add("pgm", programName);
        enableCompletionStatus();
    }
    if(ourEnableCompletionStatus)
    {
        add("status", status);
        connect();
    }
}


void PhoneHome::resetURL()
{
    ourURL = ourBaseURL;
    ourPrefixChar = '?';
}


void PhoneHome::add(const char* name, const char* val)
{
    if((name != NULL) && (strlen(name) != 0) &&
       (val != NULL) && (strlen(val) != 0))
    {
        // Check if the value is already set.
        if(ourURL.find(name) != std::string::npos)
        {
            // value already set, so do not set it.
            return;
        }

        // A value was passed in, so add it to the URL.
        ourURL += ourPrefixChar;
        ourURL += name;
        ourURL += '=';
        // If it is a tool name, trim anything before the last '/'
        if(strstr(name, "pgm") != NULL)
        {
            // toolname, so trim the val.
            const char* toolVal = strrchr(val, '/');
            if(toolVal != NULL)
            {
                toolVal++;
            }
            else
            {
                toolVal = val;
            }
            ourURL.append(toolVal);
            ourToolName =  toolVal;
        }
        else
        {
            ourURL += val;
        }
        ourPrefixChar = '&';
    }
}


bool PhoneHome::connect()
{
    if(ourNumber == -1)
    {
        srand (time(NULL));
        ourNumber = rand();
        String numString;
        numString = ourNumber;
        String thinningString;
        thinningString = allThinning;
        add("uniqNum", numString);
        add("thinning", thinningString);
    }
    if((ourNumber % 100) >= allThinning)
    {
        // Skip phoneHome.
        return(true);
    }

    // std::cerr << "url = " << ourURL << std::endl;
    ourReturnString.Clear();
    //  return(true);
#ifndef _NO_PHONEHOME
    knet_silent(1);
    knetFile *file = knet_open(ourURL.c_str(), "r");
    if (file == 0) return(false);

    const int BUF_SIZE = 100;
    char buf[BUF_SIZE];

    ssize_t readLen = BUF_SIZE-1;
    ssize_t numRead = readLen;
    while(numRead == readLen)
    {
        numRead = knet_read(file, buf, readLen);
        buf[numRead] = '\0';
        ourReturnString += buf;
    }

    knet_close(file);
    knet_silent(0);
    // std::cerr << "PhoneHome URL = " << ourReturnString.c_str() << std::endl;
#endif
    return(true);
}
