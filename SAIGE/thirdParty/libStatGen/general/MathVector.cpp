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

#include "MathVector.h"
#include "MathMatrix.h"
#include "MathConstant.h"
#include "Sort.h"
#include "Error.h"

#ifdef  _MSC_VER
#define _USE_MATH_DEFINES
#endif

#include <string.h>
#include <math.h>

int Vector::alloc = 32;

void Vector::Init()
{
    dim = size = 0;
    label = "Unknown";
    data = NULL;
}

Vector::~Vector()
{
    // printf(" Deleting vector %s ...\n", (const char *) label);
    if (data != NULL) delete [] data;
}

void Vector::Dimension(int d)
{
    if (d > size)
    {
        if (size < 1024)
        {
            size = (d + alloc) / alloc * alloc;
            double * newData = new double [size];
            if (data != NULL)
            {
                for (int i = 0; i < dim; i++)
                    newData[i] = data[i];
                delete [] data;
            }
            data = newData;
        }
        else
        {
            while (size <= d)
                size *= 2;

            double * newData = new double [size];
            if (data != NULL)
            {
                for (int i = 0; i < dim; i++)
                    newData[i] = data[i];
                delete [] data;
            }
            data = newData;
        }
    }
    dim = d;
}

void Vector::Dimension(int d, double value)
{
    int original = dim;

    Dimension(d);

    for (int i = original; i < dim; i++)
        data[i] = value;
}

void Vector::Negate()
{
    for (int i = 0; i < dim; i++)
        data[i] = -data[i];
}

void Vector::Add(double n)
{
    for (int i = 0; i< dim; i++)
        data[i] += n;
}

void Vector::Multiply(double k)
{
    for (int i = 0; i < dim; i++)
        data[i] *= k;
}

void Vector::Copy(const Vector & v)
{
    Dimension(v.dim);

    if (v.data != NULL)
        for (int i=0; i < dim; i++)
            data[i] = v.data[i];
}

Vector & Vector::operator = (const Vector & rhs)
{
    Copy(rhs);
    return *this;
}

void Vector::Add(Vector & v)
{
    if (dim != v.dim)
        error("Vector::Add - vectors have different dimensions\n"
              "Vectors     - %s [%d] + %s [%d] ",
              (const char *) label, dim, (const char  *) v.label, v.dim);

    for (int i = 0; i < dim; i++)
        data[i] += v.data[i];
}

void Vector::AddMultiple(double k, Vector & v)
{
    if (dim != v.dim)
        error("Vector::AddMultiple - vectors are incompatible\n"
              "Vectors             - %s [%d] + %s [%d] ",
              (const char  *) label, dim, (const char  *) v.label, v.dim);

    for (int i = 0; i < dim; i++)
        data[i] += k * v.data[i];
}


void Vector::Subtract(Vector & v)
{
    if (dim != v.dim)
        error("Vector::Subtract - vectors have different dimensions\n"
              "Vectors          - %s [%d] + %s [%d] ",
              (const char  *) label, dim, (const char  *) v.label, v.dim);

    for (int i = 0; i < dim; i++)
        data[i] -= v.data[i];
}


void Vector::Zero()
{
    for (int i = 0; i < dim; i++)
        data[i] = 0.0;
}

void Vector::Set(double k)
{
    for (int i = 0; i < dim; i++)
        data[i] = k;
}

void Vector::SetMultiple(double k, Vector & v)
{
    Dimension(v.dim);

    for (int i = 0; i < dim; i++)
        data[i] = k * v[i];
}

double Vector::InnerProduct(Vector & v)
{
    if (dim != v.dim)
        error("Vector::InnerProduct - vectors have different dimensions\n"
              "Vectors              - %s[%d] * %s[%d] ",
              (const char  *) label, dim, (const char  *) v.label, v.dim);

    double sum = 0.0;
    for (int i = 0; i < dim; i++)
        sum += data[i] * v.data[i];

    return sum;
}

void Vector::Insert(int n, double value)
{
    Dimension(dim + 1);

    for (int i = dim - 1; i > n; i--)
        data[i] = data[i - 1];
    data[n] = value;
}

void Vector::DeleteDimension(int n)
{
    for (int i = n; i < dim - 1; i++)
        data[i] = data[i + 1];
    dim --;
}

void Vector::Product(Matrix & m, Vector & v)
{
    if (m.cols != v.dim)
        error("Vector::Product - Cannot Multiply Matrix by Vector\n"
              "Vectors         - %s [%d, %d] * %s [%d]\n",
              (const char  *) m.label, m.rows, m.cols,
              (const char  *) v.label, v.dim);

    Dimension(m.rows);
    Zero();

    for (int i = 0; i < m.rows; i++)
        for (int j = 0; j < m.cols; j++)
            data[i] += m[i][j] * v[j];
}

double Vector::Average() const
{
    if (dim == 0)
        error("Average undefined for null vector %s",
              (const char  *) label);

    return Sum() / dim;
}

double Vector::Product() const
{
    double product = 1.0;

    for (int j = 0; j < dim; j++)
        product *= data[j];

    return product;
}

double Vector::Sum() const
{
    double sum = 0.0;

    for (int j=0; j<dim; j++)
        sum += data[j];

    return sum;
}

double Vector::SumSquares() const
{
    double sum = 0.0;

    for (int j=0; j<dim; j++)
        sum += data[j] * data[j];

    return sum;
}

void Vector::AveVar(double & ave, double & var) const
{
    // uses a two pass method to correct for
    // round-off errors

    if (dim == 0)
        error("Average and Variance undefined for null vector %s",
              (const char  *) label);

    double s, ep;

    ave = var = ep = 0.0;

    for (int j=0; j<dim; j++)
        ave += data[j];

    ave /= dim;

    for (int j=0; j<dim; j++)
    {
        s = data[j] - ave;
        ep += s;
        var += s*s;
    }

    if (dim > 1)
        var = (var - ep*ep/dim)/(dim-1);
}

double Vector::Var() const
{
    double mean, var;
    AveVar(mean, var);
    return var;
}

double Vector::StandardDeviation() const
{
    double var = Var();

    if (var < 0.0) var = 0.0;

    return sqrt(var);
}

void Vector::Print(FILE * f, int d)
{
    if (d == -1 || d > dim) d = dim;

    fprintf(f, "%.15s : ", (const char  *) label);
    for (int i = 0; i < d; i++)
        fprintf(f, "%7.3f ", data[i]);
    fprintf(f, "\n");
}

int Vector::CompareDouble(const double * a, const double * b)
{
    if (*a < *b) return -1;
    if (*a > *b) return 1;
    return 0;
}

void Vector::Sort()
{
    QuickSort(data, dim, sizeof(double), COMPAREFUNC CompareDouble);
}

void Vector::Sort(Vector & freeRider)
{
    QuickSort2(data, freeRider.data, dim, sizeof(double),
               COMPAREFUNC CompareDouble);
}

int Vector::BinarySearch(double element)
{
    void * pointer = ::BinarySearch
                     (&element, data, dim, sizeof(double), COMPAREFUNC CompareDouble);

    if (pointer == NULL)
        return -1;

    return ((double *) pointer) - data;
}

void Vector::RemoveDuplicates()
{
    int out = 0;

    for (int in = 1; in < Length(); in++)
        if (data[in] != data[out])
            data[++out] = data[in];

    Dimension(out + 1);
}

bool Vector::operator == (const Vector & rhs) const
{
    if (rhs.dim != dim) return false;

    for (int i = 0; i < dim; i++)
        if (data[i] != rhs[i])
            return false;
    return true;
}

// These functions are useful for simulation
//

int Vector::CountIfGreater(double threshold) const
{
    int count = 0;

    for (int i = 0; i < dim; i++)
        if (data[i] > threshold)
            count++;

    return count;
}

int Vector::CountIfGreaterOrEqual(double treshold) const
{
    int count = 0;

    for (int i = 0; i < dim; i++)
        if (data[i] >= treshold)
            count++;

    return count;
}

// Min and max functions
//

double Vector::Min() const
{
    if (dim == 0)
        return 0.0;

    double min = data[0];

    for (int i = 1; i < dim; i++)
        if (data[i] < min)
            min = data[i];

    return min;
}

double Vector::Max() const
{
    if (dim == 0)
        return 0.0;

    double max = data[0];

    for (int i = 1; i < dim; i++)
        if (data[i] > max)
            max = data[i];

    return max;
}

// Push and Pop functions for using vector as a stack
//

void Vector::Push(double value)
{
    Dimension(dim + 1);
    data[dim - 1] = value;
}

void Vector::Stack(const Vector & v)
{
    int end = dim;

    Dimension(dim + v.dim);

    for (int i = 0; i < v.dim; i++)
        data[i + end] = v[i];
}

// Check if all values are in ascending or descending order
//

bool Vector::isAscending()
{
    for (int i = 1; i < dim; i++)
        if (data[i] < data[i - 1])
            return false;
    return true;
}

bool Vector::isDescending()
{
    for (int i = 1; i < dim; i++)
        if (data[i] > data[i - 1])
            return false;
    return true;
}

// VectorFunc class
//

VectorFunc::VectorFunc()
{
    f = NULL;
}

VectorFunc::VectorFunc(double(*func)(Vector &))
{
    f = func;
}

double VectorFunc::Evaluate(Vector & v)
{
    return f(v);
}

#ifndef   M_SQRT2
#define   M_SQRT2        1.41421356
#endif

#define   MAXROUNDS      10
#define   SQRT_HALF      (1.0/M_SQRT2)
#define   TWO            (M_SQRT2 * M_SQRT2)

void VectorFunc::Derivative(Vector & x, Vector & d, double h_start)
{
    double a[MAXROUNDS][MAXROUNDS];

    // Calculate derivatives along each direction ...
    for (int k = 0; k < x.dim; k++)
    {
        double left, right;
        double save_x = x[k];
        double h = h_start;

        // Evaluate function to the left of x along direction k
        x[k] = save_x - h;
        left  = Evaluate(x);

        // Initialize or update dfmin if appropriate...
        if (k == 0 || left < dfmin)
            dfmin = left, dpmin = x;

        // Evaluate function to the right of x along direction k
        x[k] = save_x + h;
        right = Evaluate(x);

        // Update dfmin
        if (right < dfmin)
            dfmin = left, dpmin = x;

        // Initial crude estimate
        a[0][0] = (right - left) / (2.0 * h);

        // Initial guess of error is large
        double err = 1e30;

        // At each round, update Neville tableau with smaller stepsize and higher
        // order extrapolation ...
        for (int i = 1; i < MAXROUNDS; i++)
        {
            // Decrease h
            h *= SQRT_HALF;

            // Re-evaluate function and update dfmin as required
            x[k]  = save_x - h;
            left  = Evaluate(x);
            if (left < dfmin) dfmin = left, dpmin = x;
            x[k]  = save_x + h;
            right = Evaluate(x);
            if (right < dfmin) dfmin = right, dpmin = x;

            // Improved estimate of derivative
            a[0][i] = (right - left) / (2.0 * h);

            // Calculate extrapolations of various orders ...
            double factor = TWO;

            for (int j = 1; j <= i; j++)
            {
                a[j][i] = (a[j-1][i] * factor - a[j-1][i-1])/(factor - 1.0);

                factor *= TWO;

                double error = max(fabs(a[j][i] - a[j-1][i]), fabs(a[j][i] - a[j-1][i-1]));

                // Did we improve solution?
                if (error < err)
                {
                    err = error;
                    d[k] = a[j][i];
                }
            }

            // Stop if solution is deteriorating ...
            if (fabs(a[i][i] - a[i-1][i-1]) >= 2.0 * err)
            {
                x[k] = save_x;
                break;
            }
        }

        x[k] = save_x;
    }
}

int Vector::SafeCount() const
{
    int nonMissing = dim;

    for (int i = 0; i < dim; i++)
        if (data[i] == _NAN_)
            nonMissing--;

    return nonMissing;
}

double Vector::SafeMin() const
{
    double min = _NAN_;
    int i;

    for (i = 0; i < dim; i++)
        if (data[i] != _NAN_)
        {
            min = data[i];
            break;
        }

    for (; i < dim; i++)
        if (data[i] != _NAN_ && data[i] < min)
            min = data[i];

    return min;
}

double Vector::SafeMax() const
{
    double max = _NAN_;
    int i;

    for (i = 0; i < dim; i++)
        if (data[i] != _NAN_)
        {
            max = data[i];
            break;
        }

    for (; i < dim; i++)
        if (data[i] != _NAN_ && data[i] > max)
            max = data[i];

    return max;
}

void Vector::Reverse()
{
    for (int i = 0, j = dim - 1; i < j; i++, j--)
        Swap(i, j);
}

void Vector::InsertInSortedList(int value)
{
    // Skip through large elements
    int pos = dim - 1;

    while (pos >= 0 && data[pos] > value)
        pos--;

    // If the value is already in the list, we are done
    if (pos >= 0 && data[pos] == value)
        return;

    // Otherwise we need to grow array
    Dimension(dim + 1);

    // And then shift larger elements to the right
    pos++;
    for (int i = dim - 1; i > pos; i--)
        data[i] = data[i - 1];

    data[pos] = value;
}

void Vector::Swap(Vector & rhs)
{
    double * temp = rhs.data;
    rhs.data = data;
    data = temp;

    int swap = rhs.dim;
    rhs.dim = dim;
    dim = swap;

    swap = rhs.size;
    rhs.size = size;
    size = swap;
}

double Vector::Average(double returnIfNull)
{
    if (Length() == 0)
        return returnIfNull;

    return Average();
}

double Vector::Var(double returnIfNull)
{
    if (Length() == 0)
        return returnIfNull;

    return Var();
}

double Vector::StandardDeviation(double returnIfNull)
{
    if (Length() == 0)
        return returnIfNull;

    return StandardDeviation();
}
