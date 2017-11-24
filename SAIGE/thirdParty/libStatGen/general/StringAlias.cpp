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

#include "StringAlias.h"
#include "InputFile.h"

void StringAlias::SetAlias(String & string, String & alias)
{
    int index = lookup.Integer(string);

    if (index < 0)
    {
        aliases.Push(alias);
        lookup.SetInteger(string, aliases.Length() - 1);
    }
    else
        aliases[index] = alias;
}

const String & StringAlias::GetAlias(const String & string) const
{
    int index = lookup.Integer(string);

    if (index < 0)
        return string;
    else
        return aliases[index];
}


int StringAlias::GetAliases(StringArray & list) const
{
    if(lookup.Entries() == 0)
    {
        return 0;
    }

    int edits = 0;
    for(int i = 0; i < list.Length(); i++)
    {
        int index = lookup.Integer(list[i]);
        if(index >= 0)
        {
            list[i] = aliases[index];
            edits++;
        }
    }
    return edits;
}


bool StringAlias::ReadFromFile(const char * filename)
{
    IFILE input = ifopen(filename, "rt");

    if (input == NULL)
        return false;

    ReadFromFile(input);

    ifclose(input);

    return true;
}

bool StringAlias::ReadFromFile(IFILE & input)
{
    StringArray lines, tokens;
    lines.Read(input);

    for (int j = 0; j < lines.Length(); j++)
    {
        tokens.ReplaceTokens(lines[j]);

        if (tokens.Length() != 2) continue;

        SetAlias(tokens[0], tokens[1]);
    }

    return true;
}
