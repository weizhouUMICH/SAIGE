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

#ifndef __INTARRAY_H__
#define __INTARRAY_H__

#include <stdio.h>

class IntArray
{
private:
    int * items;
    int size, count;

    void Grow(int new_size);
    static int Compare(int * a, int * b);

public:
    static int alloc;

    IntArray(int start_size = 0);
    IntArray(const IntArray & source);
    ~IntArray();

    IntArray & operator = (const IntArray & rhs);

    int & operator [](int index)
    {
        return items[index];
    }
    int   operator [](int index) const
    {
        return items[index];
    }

    // Suggested by Anthony Berno, 12/28/06, to avoid "ambiguities" that
    // Visual Studio encountered when handling implicit conversions ...
    int & operator [](char index)
    {
        return items[int(index)];
    }
    int   operator [](char index) const
    {
        return items[int(index)];
    }
    // ... who knows whether Visual Studio makes C++ annoying to encourage C#?

    int & operator [](double fraction)
    {
        return items[(int)(count * fraction)];
    }
    int   operator [](double fraction) const
    {
        return items[(int)(count * fraction)];
    }

    int  Append(int value);
    int  Append(const IntArray & rhs);

    void Push(int value)
    {
        Append(value);
    }
    int  Pop()
    {
        return items[--count];
    }
    int  Peek() const
    {
        return items[count - 1];
    }
    int &Last() const
    {
        return items[count - 1];
    }

    void PushIfNew(int value);    // used for maintaining list without duplicates

    int  Delete(int index);
    void InsertAt(int index, int value);

    int  Find(int value) const;
    int  FastFind(int value) const
    {
        return BinarySearch(value);
    }
    int  BinarySearch(int value) const;
    void Sort();
    void Sort(IntArray & freeRider); // Sorts two arrays simultaneously

    void Zero();
    void Set(int value);
    void SetSequence(int start = 0, int increment = 1);

    int  Length() const
    {
        return count;
    }
    void Dimension(int new_count)
    {
        Grow(new_count);
        count = new_count;
    }
    void Clear()
    {
        count = 0;
    }

    int  Sum() const
    {
        return Sum(0, count - 1);
    }
    int  Sum(int start) const
    {
        return Sum(start, count - 1);
    }
    int  Sum(int start, int end) const;

    double  dSum() const
    {
        return dSum(0, count - 1);
    }
    double  dSum(int start) const
    {
        return dSum(start, count - 1);
    }
    double  dSum(int start, int end) const;

    int     SumProduct(const IntArray & weight) const;
    double  dSumProduct(const IntArray & weight) const;

    int  Max() const
    {
        return Max(0, count - 1);
    }
    int  Max(int start) const
    {
        return Max(start, count - 1);
    }
    int  Max(int start, int end) const;

    int  Min() const
    {
        return Min(0, count - 1);
    }
    int  Min(int start) const
    {
        return Min(start, count - 1);
    }
    int  Min(int start, int end) const;

    int  Count() const
    {
        return count;
    }
    int  CountIfGreater(int treshold) const;
    int  CountIfGreaterOrEqual(int treshold) const;

    void Swap(int i, int j)
    {
        int tmp = items[i];
        items[i] = items[j];
        items[j] = tmp;
    }

    void Reverse();

    operator int *()
    {
        return items;
    }

    void Add(int term);
    void Subtract(int term)
    {
        Add(-term);
    }
    void Multiply(int factor);
    void Divide(int denominator);

    void Add(const IntArray & rhs);

    IntArray & operator += (int rhs)
    {
        Add(rhs);
        return *this;
    }

    IntArray & operator += (const IntArray & rhs)
    {
        Add(rhs);
        return *this;
    }

    IntArray & operator *= (int rhs)
    {
        Multiply(rhs);
        return *this;
    }

    IntArray & operator -= (int rhs)
    {
        Add(-rhs);
        return *this;
    }

    IntArray & operator /= (int rhs)
    {
        Divide(rhs);
        return *this;
    }

    int  InnerProduct(IntArray & v);

    bool operator == (const IntArray & rhs) const;
    bool operator != (const IntArray & rhs) const;

    bool isAscending();
    bool isDescending();

    void Stack(const IntArray & rhs);

    void Swap(IntArray & rhs);

    void Print()
    {
        Print(stdout);
    }
    void Print(const char * label)
    {
        Print(stdout, label);
    }
    void Print(FILE * output);
    void Print(FILE * output, const char * label);

    int    Product();
    double DoubleProduct();

    int  Hash(int initval = 0);
};

#endif


