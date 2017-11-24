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

#ifndef __PACKEDVECTOR_H
#define __PACKEDVECTOR_H


// STL:
#include <ostream>
#include <sstream>
#include <string>

#include "Generic.h"


//
// This file implements a packed vector template based on the
// getter/setter code used in MemoryMapArray.h
//


template <
uint32_t accessorFunc(std::vector<uint8_t> &base, uint32_t index),
void setterFunc(std::vector<uint8_t> &base, uint32_t index, uint32_t value),
size_t elementCount2BytesFunc(uint32_t elementCount)
>
class PackedVector
{
protected:
    std::vector<uint8_t> m_data;
    size_t              m_elementCount;
    double              m_growthRateMultiplier;
    double              m_growthRateAdder;
public:
    PackedVector() :
        m_elementCount(0),
        m_growthRateMultiplier(1.20),
        m_growthRateAdder(128) {;}

    // accessing
    inline uint32_t operator[](uint32_t i)
    {
        return accessorFunc(m_data, i);
    }
    inline void set(uint32_t i, uint32_t v)
    {
        setterFunc(m_data, i, v);
    }

    size_t getElementCount() const
    {
        return m_elementCount;
    }

    double getUtilization() {
        return elementCount2BytesFunc(m_elementCount) / (double) m_data.capacity();
    }

    void reserve(uint32_t reserveElements) {
        m_data.reserve(elementCount2BytesFunc(reserveElements));
    }

    size_t size() {return m_elementCount;}

    void resize(uint32_t newSize) {
        m_elementCount = newSize;
        m_data.resize(elementCount2BytesFunc(m_elementCount));
    }

    // it's a bit of a challenge to optimize this...
    void push_back(uint32_t value) {
        m_elementCount++;
        if(elementCount2BytesFunc(m_elementCount) >= m_data.size()) {

            if( (elementCount2BytesFunc(m_elementCount)) > m_data.capacity())
            {
                size_t newCapacity = (size_t) (m_data.capacity() * m_growthRateMultiplier);

                // for small capacities, small fractional multipliers don't work,
                // so we check and do a linear increase in those cases:
                if(newCapacity == m_data.capacity()) {
                    newCapacity = (size_t) (m_data.capacity() + m_growthRateAdder);
                }

                m_data.reserve(newCapacity);
            }

        }
        m_data.resize(elementCount2BytesFunc(m_elementCount));
        set(m_elementCount-1, value);
    }
};

typedef PackedVector<
PackedAccess_1Bit,
PackedAssign_1Bit,
Packed1BitElementCount2Bytes
> PackedVectorBool_t;

typedef PackedVector<
PackedAccess_2Bit,
PackedAssign_2Bit,
Packed2BitElementCount2Bytes
> PackedVector2Bit_t;

typedef PackedVector<
PackedAccess_4Bit,
PackedAssign_4Bit,
Packed4BitElementCount2Bytes
> PackedVector4Bit_t;

#endif
