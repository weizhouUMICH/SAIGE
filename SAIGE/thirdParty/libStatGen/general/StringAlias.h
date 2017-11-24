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

#ifndef __STRINGALIAS_H__
#define __STRINGALIAS_H__

#include "StringArray.h"
#include "StringHash.h"

class StringAlias
{
public:
    StringAlias()               {}
    virtual ~StringAlias()      {}

    void SetAlias(String & string, String & alias);

    const String & GetAlias(const String & string) const;
    int            GetAliases(StringArray & list) const;

    bool  ReadFromFile(const char * filename);
    bool  ReadFromFile(IFILE & input);

private:
    StringIntHash  lookup;
    StringArray    aliases;
};

#endif

