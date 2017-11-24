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

#include "IntArray.h"
#include "Error.h"
#include "Hash.h"
#include "Sort.h"

#include <string.h>

int IntArray::alloc = 4;

IntArray::IntArray(int start_size)
{
    count = start_size;
    size = (count + alloc) / alloc * alloc;
    items = new int [size];
}

IntArray::IntArray(const IntArray & source)
{
    count = source.count;
    size = source.size;
    items = new int [size];

    for (int i = 0; i < count; i++)
        items[i] = source.items[i];
}

IntArray::~IntArray()
{
    delete [] items;
}

void IntArray::Grow(int new_size)
{
    if (new_size > size)
    {
        if ((new_size >> 1) >= size)
            size = (new_size + alloc) / alloc * alloc;
        else
        {
            size = alloc;
            while (size <= new_size)
                size *= 2;
        }

        int * new_items = new int [size];
        for (int i = 0; i < count; i++)
            new_items[i] = items[i];
        delete [] items;
        items = new_items;
    }
}

int IntArray::Append(int value)
{
    Grow(count + 1);
    items[count++] = value;
    return count;
}

int IntArray::Append(const IntArray & rhs)
{
    Grow(count + rhs.count);
    for (int i = 0; i < rhs.count; i++)
        items[count + i] = rhs.items[i];
    count += rhs.count;
    return count;
}

void IntArray::Set(int value)
{
    for (int i = 0; i < count; i++)
        items[i] = value;
}

void IntArray::SetSequence(int start, int increment)
{
    for (int i = 0; i < count; i++, start += increment)
        items[i] = start;
}

int IntArray::Delete(int index)
{
    count--;
    if (count - index)
        memmove(items + index, items + index + 1, sizeof(int) *(count - index));
    return count;
}

void IntArray::InsertAt(int index, int value)
{
    Grow(count + 1);
    if (count - index)
        memmove(items + index + 1, items + index, sizeof(int) *(count - index));
    items[index] = value;
    count++;
}

IntArray & IntArray::operator = (const IntArray & rhs)
{
    Grow(rhs.count);
    count = rhs.count;
    for (int i = 0; i < count; i++)
        items[i] = rhs.items[i];
    return *this;
}

int IntArray::Sum(int start, int end) const
{
    int result = 0;

    for (int i = start; i <= end; i++)
        result += items[i];

    return result;
}

double IntArray::dSum(int start, int end) const
{
    double result = 0;

    for (int i = start; i <= end; i++)
        result += items[i];

    return result;
}

int IntArray::Max(int start, int end) const
{
    if (start >= count) return 0;

    int result = items[start];

    for (int i = start + 1; i <= end; i++)
        if (result < items[i])
            result = items[i];

    return result;
}

int IntArray::Min(int start, int end) const
{
    if (start >= count) return 0;

    int result = items[start];

    for (int i = start + 1; i <= end; i++)
        if (result > items[i])
            result = items[i];

    return result;
}

int IntArray::Find(int value) const
{
    for (int i = 0; i < count; i++)
        if (value == items[i])
            return i;
    return -1;
}

int IntArray::BinarySearch(int value) const
{
    int start = 0;
    int stop = count - 1;

    while (start <= stop)
    {
        int mid = (start + stop) / 2;

        if (items[mid] == value)
            return mid;

        if (items[mid] > value)
            stop = mid - 1;
        else
            start = mid + 1;
    }

    return -1;
}

void IntArray::Zero()
{
    for (int i = 0; i < count; i++)
        items[i] = 0;
}

int IntArray::Compare(int * a, int * b)
{
    return *a - *b;
}

void IntArray::Sort()
{
    QuickSort(items, count, sizeof(int), COMPAREFUNC Compare);
}

void IntArray::Sort(IntArray & freeRider)
{
    QuickSort2(items, freeRider.items, count, sizeof(int), COMPAREFUNC Compare);
}


void IntArray::Reverse()
{
    for (int i = 0, j = count - 1; i < j; i++, j--)
        Swap(i, j);
}

int IntArray::CountIfGreater(int threshold) const
{
    int result = 0;

    for (int i = 0; i < count; i++)
        if (items[i] > threshold)
            result++;

    return result;
}

int IntArray::CountIfGreaterOrEqual(int treshold) const
{
    int result = 0;

    for (int i = 0; i < count; i++)
        if (items[i] >= treshold)
            result++;

    return result;
}

void IntArray::Add(int term)
{
    for (int i = 0; i < count; i++)
        items[i] += term;
}

void IntArray::Multiply(int factor)
{
    for (int i = 0; i < count; i++)
        items[i] *= factor;
}

void IntArray::Divide(int denominator)
{
    for (int i = 0; i < count; i++)
        items[i] /= denominator;
}

void IntArray::Stack(const IntArray & a)
{
    int end = count;

    Dimension(count + a.count);

    for (int i = 0; i < a.count; i++)
        items[i + end] = a[i];
}

bool IntArray::operator == (const IntArray & rhs) const
{
    if (count != rhs.count)
        return false;

    for (int i = 0; i < rhs.count; i++)
        if (items[i] != rhs.items[i])
            return false;

    return true;
}

bool IntArray::operator != (const IntArray & rhs) const
{
    return !(*this == rhs);
}

// Check if all values are in ascending or descending order
//

bool IntArray::isAscending()
{
    for (int i = 1; i < count; i++)
        if (items[i] < items[i - 1])
            return false;
    return true;
}

bool IntArray::isDescending()
{
    for (int i = 1; i < count; i++)
        if (items[i] > items[i - 1])
            return false;
    return true;
}

void IntArray::Add(const IntArray & v)
{
    if (Length() != v.Length())
        error("IntArray::Add - vectors have different lengths\n"
              "IntArrays     - Left[%d] += Right[%d] ",
              Length(), v.Length());

    for (int i = 0; i < Length(); i++)
        items[i] += v[i];
}

int IntArray::InnerProduct(IntArray & v)
{
    if (Length() != v.Length())
        error("IntArray::InnerProduct - vectors have different dimensions\n"
              "IntArrays              - Left[%d] * Right[%d] ",
              Length(), v.Length());

    int sum = 0;
    for (int i = 0; i < Length(); i++)
        sum += items[i] * v[i];

    return sum;
}

void IntArray::Swap(IntArray & rhs)
{
    int * temp = rhs.items;
    rhs.items = items;
    items = temp;

    int swap = rhs.count;
    rhs.count = count;
    count = swap;

    swap = rhs.size;
    rhs.size = size;
    size = swap;
}

void IntArray::Print(FILE * output)
{
    Print(output, "Array of Integers");
}

void IntArray::Print(FILE * output, const char * label)
{
    fprintf(output, "%s [%d elements]: ", label, count);

    for (int i = 0; i < count; i++)
        fprintf(output, "%d ", items[i]);

    fprintf(output, "\n");
}

void IntArray::PushIfNew(int value)
{
    for (int i = 0; i < count; i++)
        if (items[i] == value)
            return;

    Push(value);
}

int IntArray::Product()
{
    int product = 1;

    for (int i = 0; i < count; i++)
        product *= items[i];

    return product;
}

double IntArray::DoubleProduct()
{
    double product = 1.0;

    for (int i = 0; i < count; i++)
        product *= items[i];

    return product;
}

int IntArray::Hash(int initval)
{
    return hash((unsigned char *) items, sizeof(int) * count, initval);
}

int IntArray::SumProduct(const IntArray & weight) const
{
    if (count != weight.count)
        error("IntArray::SumProduct called with different sized arrays\n");

    int sum = 0;
    for (int i = 0; i < count; i++)
        sum += items[i] * weight[i];

    return sum;
}

double IntArray::dSumProduct(const IntArray & weight) const
{
    if (count != weight.count)
        error("IntArray::dSumProduct called with different sized arrays\n");

    double sum = 0.0;
    for (int i = 0; i < count; i++)
        sum += items[i] * weight[i];

    return sum;
}

