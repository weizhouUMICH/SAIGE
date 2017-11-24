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

#include "FortranFormat.h"
#include "Error.h"

FortranFormat::FortranFormat()
{
    inputPos = -1;
    endOfPattern = false;
}

void FortranFormat::SetInputFile(IFILE & file)
{
    input = file;
    inputPos = -1;
    endOfPattern = false;
}

void FortranFormat::SetFormat(const String & formatString)
{
    format = formatString;

    inputPos = -1;
    endOfPattern = false;

    repeatCount = 0;

    format.Clear();

    // Remove blank spaces from format statement and extract
    // the first bracketed expression
    int level = 0;
    for (int i = 0; i < formatString.Length(); i++)
    {
        if (formatString[i] == ' '  || formatString[i] == '\t' ||
                formatString[i] == '\n' || formatString[i] == '\r')
            continue;

        if (formatString[i] == '(')
            level++;

        if (formatString[i] == ')')
            level--;

        format += formatString[i];

        if (level == 0) break;
    }

    if (format[0] != '(' || format[format.Length() - 1] != ')')
        error("Invalid FORTRAN format statement\n\n"
              "The statement \"%s\" is not bracketed correctly.\n",
              (const char *) formatString);

    lastBracket = 1;
    lastCount = 0;

    formatPos = 1;
    repeatCount = 0;

    bracketStack.Clear();
    bracketCounter.Clear();
    bracketCount.Clear();
}

int FortranFormat::GetNextInteger()
{
    GetNextField(buffer);

    return buffer.AsInteger();
}

char FortranFormat::GetNextCharacter()
{
    GetNextField(buffer);

    return buffer[0];
}

void FortranFormat::GetNextField(String & field)
{
    while (!ProcessToken(field))
        ;
}

bool FortranFormat::ProcessToken(String & field)
{
    // This flag only gets set if we encounter the final bracket or a ':'
    endOfPattern = false;

    // Read input from file, if appropriate
    if (inputPos == -1)
    {
        inputLine.ReadLine(input);
        inputPos = 0;
    }

    // First read repeat count specifier
    if (repeatCount == 0)
        repeatCount = GetIntegerFromFormat();

    // By default, the repeat count should be 1
    if (repeatCount == 0)
        repeatCount = 1;

    int repeatPos = formatPos;

    // Check if this is a new bracketed grouping
    if (format[formatPos] == '(')
    {
        formatPos++;

        bracketStack.Push(formatPos);
        bracketCounter.Push(repeatCount);
        bracketCount.Push(repeatCount);

        repeatCount = 0;

        return false;
    }

    // Check if this an 'X' field
    if (format[formatPos] == 'X')
    {
        formatPos++;

        // No width specifier allowed for these fields
        RejectWidth('X');

        // Skip appropriate number of characters
        inputPos += repeatCount;

        // Reset repeat count
        repeatCount = 0;

        FinishField();

        return false;
    }

    // Check if this is a '/' (vertical tab field)
    if (format[formatPos] == '/')
    {
        formatPos++;

        // No width specifier allowed for these fields
        RejectWidth('/');

        // Skip the appropriate number of lines
        while (repeatCount--)
            inputLine.ReadLine(input);

        inputPos = 0;

        // Separators are optional, so we might already be at the next field
        if (format[formatPos] == ',' || format[formatPos] || ')')
            FinishField();

        return false;
    }

    // Check that we haven't encountered a rare, but unsupported input type
    if (format[formatPos] == 'Q' || format[formatPos] == 'P' || format[formatPos] == 'B')
    {
        formatPos++;

        int problemStart = formatPos;

        while (format[formatPos] != ',' && format[formatPos] != ')' && format[formatPos] != '/')
            formatPos++;

        error("Unsupported pattern in FORMAT statement\n\n"
              "Statement \"%s\" includes unsupporterd pattern '%s'\n",
              (const char *) format,
              (const char *) format.SubStr(problemStart, formatPos - problemStart));
    }

    if (format[formatPos] == ':')
    {
        formatPos++;

        if (format[formatPos] == ',' || format[formatPos] || ')')
            FinishField();

        repeatCount = 0;

        endOfPattern = true;

        return false;
    }

    // All the other types we recognize include a width specifier

    // Identify the location of the type specifier
    int typeStart = formatPos;

    while (CharacterFollows())
        formatPos++;

    int typeLen = formatPos - typeStart;

    // Retrieve the field width
    int width = GetIntegerFromFormat();

    if (width == 0)
        error("Unrecognized FORMAT statement\n\n"
              "Statement \"%s\" is missing a width specifier for a field of type '%s'\n",
              (const char *) format, (const char *) format.SubStr(typeStart, typeLen));

    // Check for horizontal tab character
    if (format[typeStart] == 'T')
    {
        // Move left by a specified number of characters
        if (format[typeStart + 1] == 'L')
            inputPos = width > inputPos ? 0 : inputPos - width;
        // Move right by a specified number of characters
        else if (format[typeStart + 1] == 'R')
            inputPos += width;
        // Or simply set the appropriate horizontal position
        else
            inputPos = width;

        repeatCount--;

        if (repeatCount)
            formatPos = repeatPos;
        else
            FinishField();

        return false;
    }

    // Assume that if we got here, we are looking at a data field!
    field.Copy(inputLine, inputPos, width);
    field.Trim();

    inputPos += width;

    repeatCount--;

    if (repeatCount)
        formatPos = repeatPos;
    else
        FinishField();

    return true;
}

int FortranFormat::GetIntegerFromFormat()
{
    int result = 0;

    while (DigitFollows())
        result = result * 10 + (int)(format[formatPos++] - '0');

    return result;
}

bool FortranFormat::DigitFollows()
{
    return (format[formatPos] >= '0') && (format[formatPos] <= '9');
}

bool FortranFormat::CharacterFollows()
{
    return (format[formatPos] >= 'A') && (format[formatPos] <= 'Z');
}

void FortranFormat::RejectWidth(char ch)
{
    // No width allowed for field types 'X' and '\'
    if (DigitFollows())
        error("Unrecognized FORTRAN format statement\n\n"
              "The statement \"%s\" includes width specifier for field of type '%c'.\n",
              (const char *) format, ch);
}

void FortranFormat::FinishField(bool)
{
    // Find the next field separator
    while (format[formatPos] != ',' && format[formatPos] != ')')
    {
        if (format[formatPos] == '/')
            return;

        formatPos++;
    }

    // Skip commas
    if (format[formatPos] == ',')
    {
        formatPos++;
        return;
    }

    // If we found a bracket, then it is either the end of the statement
    // (if bracketStack is empty) or we finish an internal grouping
    if (bracketStack.Length())
    {
        // Retrieve information about this grouping
        lastBracket = bracketStack.Pop();
        lastCount = bracketCount.Pop();
        int lastCounter = bracketCounter.Pop() - 1;

        // Loop if required
        if (lastCounter)
        {
            bracketStack.Push(lastBracket);
            bracketCount.Push(lastCount);
            bracketCounter.Push(lastCounter);

            formatPos = lastBracket;
        }
        else
            // Otherwise find the next separator
        {
            formatPos++;
            FinishField();
            return;
        }
    }
    else
    {
        // If we finished the input line, then activate reset input counter
        inputPos = -1;
        endOfPattern = true;

        // And re-use input tokens starting at the last bracket
        formatPos = lastBracket;

        if (lastBracket == 1)
            return;

        // With appropriate repeat counts
        bracketStack.Push(lastBracket);
        bracketCounter.Push(lastCount);
        bracketCount.Push(lastCount);
    }
}

void FortranFormat::Flush()
{
    while (!endOfPattern)
        ProcessToken(buffer);

    inputPos = -1;

    lastBracket = 1;
    lastCount = 0;

    formatPos = 1;
    repeatCount = 0;

    bracketStack.Clear();
    bracketCounter.Clear();
    bracketCount.Clear();
}
