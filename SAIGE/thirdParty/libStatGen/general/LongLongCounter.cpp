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

#include "LongLongCounter.h"

LongCounter::LongCounter() : LongHash<int>()
{
    SetAllowDuplicateKeys(false);
}

void LongCounter::IncrementCount(long long key)
{
    unsigned int slot = Find(key);

    if (slot == LH_NOTFOUND)
        Add(key, 1);
    else if (Object(slot) == -1)
        Delete(slot);
    else
        Object(slot)++;
}

void LongCounter::DecrementCount(long long key)
{
    unsigned int slot = Find(key);

    if (slot == LH_NOTFOUND)
        Add(key, -1);
    else if (Object(slot) == 1)
        Delete(slot);
    else
        Object(slot)--;
}

int LongCounter::GetCount(long long key)
{
    unsigned int slot = Find(key);

    if (slot == LH_NOTFOUND)
        return 0;
    else
        return Object(slot)--;
}


