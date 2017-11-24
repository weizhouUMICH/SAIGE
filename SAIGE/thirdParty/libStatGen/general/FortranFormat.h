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

#ifndef __FORTRAN_FORMAT__
#define __FORTRAN_FORMAT__

#include "StringBasics.h"
#include "IntArray.h"

class FortranFormat
{
public:
    // This class reads a user specified input file, one line at a time,
    // and returns individual fields according to a user specified format
    // statement
    FortranFormat();

    // Set the fortran format statement
    void SetFormat(const String & formatString);

    // Set the input file
    void SetInputFile(IFILE & file);

    // Read one field from input file
    void GetNextField(String & field);
    int  GetNextInteger();
    char GetNextCharacter();

    // Process a token in format statement and return true
    // if token corresponds to input field. Return false if
    // token led to processing of white-space or input line
    // positioning
    bool ProcessToken(String & field);

    // Flush the pattern -- this finishes processing the current
    // pattern and ensures that all trailing new-lines, etc. are
    // handled correctly
    void Flush();

private:
    // The input line and current position along it
    String inputLine;
    int inputPos;

    // The Fortran format statement and current position along it
    String format;
    int formatPos;

    // The position of the pattern we are repeating, if any
    int repeatCount;

    // Returns an integer from the current format statement, if any
    int GetIntegerFromFormat();

    // These functions check the next character in format string
    bool DigitFollows();
    bool CharacterFollows();

    // This function finish the input field
    void FinishField(bool haveSlash = false);

    // Reject width were appropriate
    void RejectWidth(char type);

    // The input file
    IFILE input;

    // Stacks to keep track of nested parenthesis
    IntArray bracketStack;
    IntArray bracketCount;
    IntArray bracketCounter;

    int lastBracket;
    int lastCount;

    // Buffer for reading fields
    String buffer;

    // Flag that indicates whether we have reached end-of-pattern
    bool   endOfPattern;
};

#endif


