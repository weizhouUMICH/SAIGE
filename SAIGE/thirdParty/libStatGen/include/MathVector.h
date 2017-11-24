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

#ifndef __MATHVECTOR_H__
#define __MATHVECTOR_H__

#include "StringBasics.h"

#include <stdio.h>
#include <assert.h>

class Matrix;

class Vector
{
public:
    int         dim, size;
    double *    data;
    String      label;

    Vector()
    {
        Init();
    }
    Vector(Vector & v)
    {
        Init();
        Copy(v);
    }
    Vector(int d)
    {
        Init();
        Dimension(d);
    }
    Vector(const char * text)
    {
        Init();
        label = text;
    }
    Vector(const char * text, int d)
    {
        Init();
        label = text;
        Dimension(d);
    }
    Vector(const char * text, Vector & v)
    {
        Init();
        label = text;
        Copy(v);
    }

    ~Vector();

    void   Dimension(int d);
    void   Dimension(int d, double value);

    void   GrowTo(int d)
    {
        Dimension(d > dim ? d : dim);
    }
    void   GrowTo(int d, double value)
    {
        Dimension(d > dim ? d : dim, value);
    }

    int    Length() const
    {
        return dim;
    }

    void   SetLabel(const char * text)
    {
        label = text;
    }

    void   Zero();
    void   Set(double k);
    void   Set(Vector & v)
    {
        Copy(v);
    };
    void   SetMultiple(double k, Vector & v);

    void   Negate();
    void   Add(double n);
    void   Multiply(double k);

    double InnerProduct(Vector & v);
    void   Copy(const Vector & v);
    void   Add(Vector & v);
    void   AddMultiple(double k, Vector & v);
    void   Subtract(Vector & v);

    void Product(Matrix & m, Vector & v);

    double & operator [](int n)
    {
        assert(n < dim);
        return data[n];
    }
    double operator [](int n) const
    {
        assert(n < dim);
        return data[n];
    }

    double operator [](double fraction)
    {
        return data[(int)(dim * fraction)];
    }
    double & operator [](double fraction) const
    {
        return data[(int)(dim * fraction)];
    }

    Vector & operator = (const Vector & v);
    bool operator == (const Vector & v) const;
    bool operator != (const Vector & v) const
    {
        return !(*this == v);
    }

    void Swap(int i, int j)
    {
        double swap = data[i];
        data[i] = data[j];
        data[j] = swap;
    }
    void Swap(Vector & rhs);

    Vector & operator *= (double rhs)
    {
        Multiply(rhs);
        return *this;
    }
    Vector & operator += (double rhs)
    {
        Add(rhs);
        return *this;
    }
    Vector & operator -= (double rhs)
    {
        return *this += -rhs;
    }
    Vector & operator /= (double rhs)
    {
        return *this *= 1/rhs;
    }

    Vector & operator += (Vector & rhs)
    {
        Add(rhs);
        return * this;
    }
    Vector & operator -= (Vector & rhs)
    {
        Subtract(rhs);
        return * this;
    }

    void DeleteDimension(int n);
    void Delete(int n)
    {
        DeleteDimension(n);
    }
    void Insert(int n, double value);

    // Calculates average and variance
    void   AveVar(double & ave, double & var) const;
    double Average() const;
    double Var() const;
    double StandardDeviation() const;

    double Average(double returnIfNull);
    double Var(double returnIfNull);
    double StandardDeviation(double returnIfNull);

    // Common descriptive functions
    double Sum() const;
    double SumSquares() const;
    double Product() const;

    // Find extreme values
    double Min() const;
    double Max() const;

    // Return the number of elements in a subset
    int  CountIfGreater(double treshold) const;
    int  CountIfGreaterOrEqual(double treshold) const;

    // Append another vector to the end
    void Stack(const Vector & v);

    void Print(int maxDim = -1)
    {
        Print(stdout, maxDim);
    }
    void Print(FILE * output, int maxDim = -1);

    // Routines for creating and searching through sorted vectors
    void Sort();
    void Reverse();
    void Sort(Vector & freeRider);
    int  BinarySearch(double value);
    int  FastFind(double value)
    {
        return BinarySearch(value);
    }

    // Remove consecutive duplicate elements from vector
    void RemoveDuplicates();

    // Query first and last elements
    //

    double & First()
    {
        return data[0];
    }
    double & Last()
    {
        return data[dim - 1];
    }

    // Routines for using a vector as a stack of doubles
    //

    void   Clear()
    {
        dim = 0;
    }
    void   Push(double value);
    double Pop()
    {
        return data[--dim];
    }
    double Peek() const
    {
        return data[dim-1];
    }

    // This routine adds items to a sorted list
    //

    void   InsertInSortedList(int item);

    static int alloc;

    bool   isAscending();
    bool   isDescending();

    // Routines for dealing with vectors that include missing data
    //

    int SafeCount() const;
    double SafeMin() const;
    double SafeMax() const;

private:
    static int CompareDouble(const double * a, const double * b);
    void Init();
};



class VectorFunc
// Wrapper for multi-dimensional functions
// so that they can be used as parameters
// and keep private data
{
private:
    double(*f)(Vector &);

public:
    // Constructors
    VectorFunc();
    VectorFunc(double(*func)(Vector &));

    // Virtual destructor ensures that dynamic objects are
    // handled correctly
    virtual ~VectorFunc() { }

    virtual double Evaluate(Vector & v);

    // Calculate derivatives along each direction. Delta is a guess value
    // for the initial stepsize in numerical derivation
    virtual void   Derivative(Vector & point, Vector & d, double delta = 1.0);

    // Minimum function value found while evaluating derivative
    // and its location...
    double dfmin;
    Vector dpmin;
};

#endif



