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

#include "MathMatrix.h"
#include "MathVector.h"
#include "MathConstant.h"
#include "Sort.h"
#include "Error.h"

#include <string.h>
#include <math.h>
#include <stdio.h>

int Matrix::alloc = 2;

Matrix::~Matrix()
{
    // printf("Deleting Matrix %s...\n", (const char *) label);

    for (int i=0; i<size; i++)
        delete data[i];

    if (size)
        delete [] data;

    if (extraSize)
        delete [] extras;
}

void Matrix::Init()
{
    rows = cols = extraSize = size = 0;
    data = NULL;
    extras = NULL;
    label = "[Matrix]";
}

void Matrix::SetLabel(const char * name)
{
    label = '[';
    label += name;
    label += ']';
}

void Matrix::Dimension(int m, int n)
{
    if (n == cols && m == rows)
        return;

    if (n > extraSize)
    {
        int newSize = (n + alloc) / alloc  * alloc;
        ColumnExtras * newExtras = new ColumnExtras [newSize];

        if (extras != NULL)
            for (int i = 0; i < extraSize; i++)
                newExtras[i] = extras[i];

        if (extraSize)
            delete [] extras;

        extraSize = newSize;
        extras = newExtras;
    }

    if (m > size)
    {
        int       newSize = (m + alloc) / alloc * alloc;
        Vector ** newData = new Vector * [newSize];

        if (data != NULL)
            for (int i = 0; i < size; i++)
                newData[i] = data[i];

        for (int i = size; i < newSize; i++)
            newData[i] = new Vector(n);

        if (size)
            delete [] data;

        size = newSize;
        data = newData;
    }

    if (cols != n)
        for (int i = 0; i < rows; i++)
            data[i]->Dimension(n);

    if (rows != m)
        for (int i = rows; i < m; i++)
            data[i]->Dimension(n);

    rows = m;
    cols = n;
}

void Matrix::Dimension(int m, int n, double value)
{
    int originalRows = rows;
    int originalColumns =  cols;

    Dimension(m, n);

    if (rows > originalRows)
        for (int i = originalRows; i < rows; i++)
            data[i]->Set(value);

    if (cols > originalColumns)
        for (int i = 0; i < originalRows; i++)
            for (int j = originalColumns; j < cols; j++)
                data[i]->data[j] = value;
}

void Matrix::Zero()
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            (*(data[i]))[j] = 0.0;
}

void Matrix::Identity()
{
    if (rows != cols)
        error("Matrix.Identity - Identity matrices must be square");

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            if (i == j)
                (*(data[i]))[j] = 1.0;
            else
                (*(data[i]))[j] = 0.0;
}

void Matrix::Set(double k)
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            (*(data[i]))[j] = k;
}

void Matrix::Negate()
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            (*(data[i]))[j] = -(*(data[i]))[j];
}

void Matrix::Copy(const Matrix & m)
{
    Dimension(m.rows, m.cols);

    if (m.data != NULL)
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                (*this)[i][j] = m[i][j];
}

void Matrix::Transpose(const Matrix & m)
{
    Dimension(m.cols, m.rows);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            (*(data[i]))[j] = m[j][i];
}

void Matrix::Add(double k)
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            (*(data[i]))[j] += k;
}

void Matrix::Multiply(double k)
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            (*(data[i]))[j] *= k;
}

void Matrix::Add(const Matrix & m)
{
    if ((rows != m.rows) && (cols != m.cols))
        error("Matrix.Add - Attempted to add incompatible matrices\n"
              "Matrices   - %s [%d, %d] + %s [%d, %d]\n",
              (const char  *) label, rows, cols,
              (const char  *) m.label, m.rows, m.cols);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            (*(data[i]))[j] += m[i][j];
}

void Matrix::AddMultiple(double k, const Matrix & m)
{
    if ((rows != m.rows) && (cols != m.cols))
        error("Matrix.AddMultiple - Attempted to add incompatible matrices\n"
              "Matrices           - %s [%d, %d] + k * %s [%d, %d]\n",
              (const char  *) label, rows, cols,
              (const char  *) m.label, m.rows, m.cols);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            (*(data[i]))[j] += k * m[i][j];
}


void Matrix::Product(const Matrix & l, const Matrix & r)
{
    if (l.cols != r.rows)
        error("Matrix.Multiply - Attempted to multiply incompatible matrices\n"
              "Matrices        - %s [%d, %d] + %s [%d, %d]\n",
              (const char  *) l.label, l.rows, l.cols,
              (const char  *) r.label, r.rows, r.cols);

    Dimension(l.rows, r.cols);
    Zero();

    for (int k = 0; k < l.cols; k++)
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                (*(data[i]))[j] += l[i][k] * r[k][j];
}

void Matrix::AddRows(double k, int r1, int r2)
{
    Vector v(*(data[r1]));

    v.Multiply(k);
    data[r2]->Add(v);
}

void Matrix::MultiplyRow(int r1, double k)
{
    data[r1]->Multiply(k);
}

void Matrix::AddRows(int r1, int r2)
{
    data[r2]->Add(*(data[r1]));
}

void Matrix::Reduce(double tol)
{
    double pivot;
    int pivotr = 0;      // Initializing pivotr is not necessary, but avoids warnings
    int r = 0;           // the row we are currently reducing

    for (int j = 0; j < cols; j++)
    {
        if (r > rows)
            return;

        pivot = 0.0;
        for (int i = r; i < rows; i++)
            if (fabs((*this)[i][j]) > pivot)
            {
                pivot = fabs((*this)[i][j]);
                pivotr = i;
            }

        if (pivot <= tol)
        {
            for (int i = r; i < rows; i++)
                (*this)[i][j] = 0.0;
            continue;
        }

        SwapRows(pivotr, r);

        double scale = (*this)[r][j];

        (*this)[r][j] = 1.0;
        for (int k = j+1; k < cols; k++)
            (*this)[r][k] /= scale;

        for (int i = r + 1; r < rows; i++)
        {
            scale = (*this)[i][j];
            (*this)[i][j] = 0.0;
            for (int k = j+1; k < cols; k++)
                (*this)[i][k] -= (*this)[r][k] * scale;
        }

        r++;
    }
}

void Matrix::DeleteRow(int r)
{
    Vector * temp = data[r];

    for (int i = r + 1; i < rows; i++)
        data[i-1] = data[i];

    data[rows - 1] = temp;
    rows--;
}

void Matrix::DeleteColumn(int c)
{
    for (int i = 0; i < rows; i++)
        data[i] -> DeleteDimension(c);

    for (int i = c + 1; i < cols; i++)
        extras[i-1] = extras[i];

    cols--;
}

void Matrix::SwapColumns(int c1, int c2)
{
    double temp;

    for (int i = 0; i < rows; i++)
    {
        temp = (*data[i])[c1];
        (*data[i])[c1] = (*data[i])[c2];
        (*data[i])[c2] = temp;
    }

    extras[c1].Swap(extras[c2]);
}

void Matrix::Read(FILE * f)
{
    int  r, c;
    char buffer[100];
    int numItems = 0;

    numItems = fscanf(f, " %s =", buffer);
    if(numItems != 1) { }
    buffer[strlen(buffer) - 1] = 0;
    SetLabel(buffer);

    numItems = fscanf(f, " [ %d x %d ]", &r, &c);
    if(numItems != 2) { }
    Dimension(r, c);

    for (int c = 0; c < cols; c++)
    {
        numItems = fscanf(f, " %s", buffer);
        if(numItems != 1) { }
        SetColumnLabel(c, buffer);
    }

    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++)
        {
            numItems = fscanf(f, " %lf", &((*this)[r][c]));
            if(numItems != 1) { }
        }
}


void Matrix::Print(FILE * f, int r, int c)
{
    if (r == -1 || r > rows) r = rows;
    if (c == -1 || c > cols) c = cols;

    char dimensions[30];

    sprintf(dimensions, "[%d x %d]", r, c);

    int columnZero = label.Length() > 15 ? label.Length() : 15;

    fprintf(f, "\n%*s =\n%*s ", columnZero, (const char  *) label,
            columnZero, dimensions);

    int * precision = new int [c + 1];
    int * width = new int [c + 1];

    for (int j = 0; j < c; j++)
    {
        precision[j] = extras[j].GetPrecision();
        width[j] = extras[j].GetWidth();
        fprintf(f, "%*s ", width[j], (const char  *) extras[j].label);
    }

    for (int i = 0; i < r; i++)
    {
        fprintf(f, "\n%*s ", columnZero, (const char  *) data[i]->label);
        for (int j = 0; j < c; j++)
            fprintf(f, "%*.*f ", width[j], precision[j], (*this)[i][j]);
    }

    fprintf(f, "\n");

    delete [] precision;
    delete [] width;
}

void Matrix::CopyLabels(Matrix & M)
{
    for (int i = 0; i < rows; i++)
        if (i < M.rows)
            data[i]->SetLabel(M[i].label);

    for (int i = 0; i < cols; i++)
        if (i < M.cols)
            SetColumnLabel(i, M.GetColumnLabel(i));
}

// ColumnExtras class
//

void ColumnExtras::Init()
{
    label = "column";
    dirty = true;
    precision = 3;
    width = 7;
}

ColumnExtras::~ColumnExtras()
{ }

void ColumnExtras::SetLabel(const char * name)
{
    label = name;
}

int ColumnExtras::GetWidth()
{
    if (dirty)
    {
        if (precision + 2 > width)
            width = precision + 2;
        if (label.Length() > width)
            width = label.Length();
        dirty = false;
    }
    return width;
}

void ColumnExtras::Copy(ColumnExtras & c)
{
    width = c.width;
    precision = c.precision;
    dirty = c.dirty;
    label = c.label;
}

#define SWAP(a,b)      {int swap=(a); (a)=(b); (b)=swap;}
#define SWAPBOOL(a,b)  {bool swap=(a); (a)=(b); (b)=swap;}

void ColumnExtras::Swap(ColumnExtras & c)
{
    SWAP(c.width, width);
    SWAP(c.precision, precision);
    SWAPBOOL(c.dirty, dirty);
    c.label.Swap(label);
}

int Matrix::CompareRows(Vector ** row1, Vector ** row2)
{
    if ((**row1)[0] < (**row2)[0]) return -1;
    if ((**row1)[0] > (**row2)[0]) return 1;
    return 0;
}

void Matrix::Sort()
{
    QuickSort(data, rows, sizeof(Vector *), COMPAREFUNC CompareRows);
}

bool Matrix::operator == (const Matrix & rhs) const
{
    if (rhs.rows != rows || rhs.cols != cols) return false;

    for (int i = 0; i < rows; i++)
        if ((*this)[i] != rhs[i])
            return false;
    return true;
}

void Matrix::StackBottom(const Matrix & m)
{
    if (m.cols != cols)
        error("Attempted to stack matrices with different number of columns");

    int end = rows;

    Dimension(rows + m.rows, cols);

    for (int i = 0; i < m.rows; i++)
        *(data[i + end]) = m[i];
}

void Matrix::StackLeft(const Matrix & m)
{
    if (m.rows != rows)
        error("Attempted to stack matrics with different numbers of rows");

    for (int i = 0; i < rows; i++)
        data[i]->Stack(m[i]);

    Dimension(rows, cols + m.cols);
}

void Matrix::Swap(Matrix & m)
{
    label.Swap(m.label);

    ColumnExtras * tmpExtras = extras;
    extras = m.extras;
    m.extras = tmpExtras;

    int swap;
    swap = rows;
    rows = m.rows;
    m.rows = swap;
    swap = cols;
    cols = m.cols;
    m.cols = swap;
    swap = size;
    size = m.size;
    m.size = swap;
    swap = extraSize;
    extraSize = m.extraSize;
    m.extraSize = swap;

    Vector ** tmpData = data;
    data = m.data;
    m.data = tmpData;
}

double Matrix::Min() const
{
    if (rows == 0 || cols == 0)
        return 0.0;

    double minimum = data[0]->Min();

    for (int i = 1; i < rows; i++)
        minimum = min(data[i]->Min(), minimum);

    return minimum;
}

double Matrix::Max() const
{
    if (rows == 0 || cols == 0)
        return 0.0;

    double maximum = data[0]->Max();

    for (int i = 1; i < rows; i++)
        maximum = max(data[i]->Max(), maximum);

    return maximum;
}

double Matrix::Mean() const
{
    if (rows == 0 || cols == 0)
        return 0.0;

    double sum = data[0]->Sum();

    for (int i = 1; i < rows; i++)
        sum += data[i]->Sum();

    return sum / (rows * cols);
}

double Matrix::SafeMin() const
{
    double lo = (rows > 0 && cols > 0) ? _NAN_ : 0.0;

    int i, j;

    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
            if (data[i]->data[j] != _NAN_)
            {
                lo = data[i]->data[j];
                break;
            }
        if (j != cols) break;
    }

    for (; i < rows; i++, j = 0)
        for (; j < cols; j++)
            if (data[i]->data[j] < lo && data[i]->data[j] != _NAN_)
                lo = data[i]->data[j];

    return lo;
}

double Matrix::SafeMax() const
{
    double hi = (rows > 0 && cols > 0) ? _NAN_ : 0.0;

    int i, j;

    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
            if (data[i]->data[j] != _NAN_)
            {
                hi = data[i]->data[j];
                break;
            }
        if (j != cols) break;
    }

    for (; i < rows; i++, j = 0)
        for (; j < cols; j++)
            if (data[i]->data[j] > hi && data[i]->data[j] != _NAN_)
                hi = data[i]->data[j];

    return hi;
}

double Matrix::SafeMean() const
{
    double sum = 0.0;
    int count = 0;

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            if ((*this)[i][j] != _NAN_)
            {
                sum += (*this)[i][j];
                count ++;
            }

    return (count) ? sum / count : 0.0;
}

int Matrix::SafeCount() const
{
    int total = 0;

    for (int i = 0; i < rows; i++)
        total += data[i]->SafeCount();

    return total;
}

void Matrix::PrintUpper(FILE * f, int r, int c, bool print_diag)
{
    int columnZero;
    int * precision = NULL, * width = NULL; // Initialization avoids compiler warnings

    SetupPrint(f, r, c, columnZero, precision, width);

    int upper = print_diag ? 0 : 1;
    for (int i = 0; i < r ; i++)
    {
        fprintf(f, "\n%*s ", columnZero, (const char  *) data[i]->label);

        for (int j = 0; j < upper; j++)
            fprintf(f, "%*.*s ", width[j], precision[j], " ");
        for (int j = upper; j < c; j++)
            fprintf(f, "%*.*f ", width[j], precision[j], (*this)[i][j]);

        upper++;
    }

    fprintf(f, "\n");

    delete [] precision;
    delete [] width;
}

void Matrix::PrintLower(FILE * f, int r, int c, bool print_diag)
{
    if (r == -1 || r > rows) r = rows;
    if (c == -1 || c > cols) c = cols;

    String dimensions;
    dimensions.printf("[%d x %d]", r, c);

    int columnZero = label.Length() > 15 ? label.Length() : 15;

    fprintf(f, "\n%*s =\n%*s ", columnZero, (const char  *) label,
            columnZero, (const char *) dimensions);

    int * precision = new int [c + 1];
    int * width = new int [c + 1];

    for (int j = 0; j < c; j++)
    {
        precision[j] = extras[j].GetPrecision();
        width[j] = extras[j].GetWidth();
        fprintf(f, "%*s ", width[j], (const char  *) extras[j].label);
    }

    int upper = print_diag ? 1 : 0;

    for (int i = 0; i < r ; i++)
    {
        fprintf(f, "\n%*s ", columnZero, (const char  *) data[i]->label);
        for (int j = 0; j < upper; j++)
            fprintf(f, "%*.*f ", width[j], precision[j],(*this)[i][j]);
        for (int j = upper; j < c; j++)
            fprintf(f, "%*.*s ", width[j], precision[j]," ");

        upper++;
    }

    fprintf(f, "\n");

    delete [] precision;
    delete [] width;
}


void Matrix::SetupPrint(FILE *f, int r, int c, int & column_zero, int * precision, int * width)
{
    if (r == -1 || r > rows) r = rows;
    if (c == -1 || c > cols) c = cols;

    String dimensions;
    dimensions.printf("[%d x %d]", r, c);

    column_zero = label.Length() > 15 ? label.Length() : 15;

    fprintf(f, "\n%*s =\n%*s ", column_zero, (const char  *) label,
            column_zero, (const char *) dimensions);

    precision = new int [c + 1];
    width = new int [c + 1];

    for (int j = 0; j < c; j++)
    {
        precision[j] = extras[j].GetPrecision();
        width[j] = extras[j].GetWidth();
        fprintf(f, "%*s ", width[j], (const char  *) extras[j].label);
    }
}
