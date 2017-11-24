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

#ifndef __PHONEHOME_H__
#define __PHONEHOME_H__

#include <string>
#include "StringBasics.h"


// By default, CompletionStatus PhoneHome is disabled.
// To enable it:
//    1) call "enableCompletionStatus"
//    2) call checkVersion
//    3) call completionStatus with the program name passed in
//      
class PhoneHome
{
public:
    // Public method that can be set to control the thinning of version checks.
    static int allThinning;

    // Enable Completion Status PhoneHome, it is disabled by default.
    // It can also be enabled by:
    //    * calling checkVersion
    //    * calling completionStatus with the program name passed in
    // Program name must be specified in order to log completionStatus
    static void enableCompletionStatus(const char* programName = NULL);
    
    // Disable Completion Status PhoneHome. (It is already disabled by default.)
    static void disableCompletionStatus();

    // Check the version, printing a message if a newer version is available.
    // Enables CompletionStatus PhoneHome
    // Returns false if there is a new version available, otherwise true.
    static bool checkVersion(const char* programName,
                             const char* version,
                             const char* params = NULL);

    // If completionStatus is enabled, send the completion status.
    // completionStatus is enabled if:
    //     1) enableCompletionStatus was called
    //     2) checkVersion was called
    //     3) programName is passed in
    // ProgramName is ignored if it has previously been set.
    static void completionStatus(const char* status,
                                 const char* programName = NULL);

    static void setURL(const char* url);
    static void resetURL();

protected:
private:
    static void add(const char* name, const char* val);
    static bool connect();

    static bool ourEnableCompletionStatus;
    static std::string ourBaseURL;
    static std::string ourURL;
    static char ourPrefixChar;
    static int ourNumber;
    static String ourToolName;

    static String ourReturnString;
};

#endif
