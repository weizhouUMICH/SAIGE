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

#include "STLUtilities.h"

namespace STLUtilities
{

//
// Split the string input into words delimited by the character
// delimiter.  For a given number of input delimiters, result.size()
// will not change, regardless of the data in between the delimiters.
//
// Refactor this to pre-allocate the word that we place data into,
// then we have minimal data copy.
//
int Tokenize(std::vector<std::string> &result, const char *input, char delimiter)
{
    if (*input=='\0')
    {
        result.clear();
        result.resize(1);   // one word, and it is empty
        return 0;
    }

    size_t wordCount = 1;

    // since input is non-empty, we know we will have at least
    // one word, so we allocate it here, and begin to fill it in
    if (result.size()<wordCount) result.resize(1);
    else result[0].clear();

    std::string *word = &result[0];

    while (*input)
    {
        if (*input==delimiter)
        {
            // we got a delimeter, and since an empty word following
            // a delimeter still counts as a word, we allocate it here
            wordCount++;
            if (result.size()<wordCount) result.resize(wordCount);
            else
            {
                result[wordCount-1].clear();
            }
            word = &result[wordCount-1];
        }
        else
        {
            // save the char in this word
            word->push_back(*input);
        }
        input++;
    }

    if (wordCount < result.size()) result.resize(wordCount);  // potentially truncate to wordCount elements

    return result.size();
}

} // end of namespace STLUtilities
