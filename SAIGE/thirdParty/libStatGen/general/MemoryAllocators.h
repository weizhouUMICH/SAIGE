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

#ifndef __MEMORY_ALLOCATORS_H__
#define __MEMORY_ALLOCATORS_H__

#include <stdlib.h>

template <class T> T** AllocateMatrix(int rows, int cols);
template <class T> T** AllocateMatrix(int rows, int cols, T value);
template <class T> void FreeMatrix(T ** & matrix, int rows);

char **  AllocateCharMatrix(int rows, int cols);
void     FreeCharMatrix(char ** & matrix, int rows);

float ** AllocateFloatMatrix(int rows, int cols);
void     FreeFloatMatrix(float ** & matrix, int rows);

double ** AllocateDoubleMatrix(int rows, int cols);
void     FreeDoubleMatrix(double ** & matrix, int rows);

int  **  AllocateIntMatrix(int rows, int cols);
void     FreeIntMatrix(int ** & matrix, int rows);

short ** AllocateShortMatrix(int rows, int cols);
void     FreeShortMatrix(short ** & matrix, int rows);

char *** AllocateCharCube(int n, int rows, int cols);
void     FreeCharCube(char *** & matrix, int n, int rows);


// Template definitions follow ...
//

template <class T> T** AllocateMatrix(int rows, int cols)
{
    T ** matrix = new T * [rows];

    // Stop early if we are out of memory
    if (matrix == NULL)
        return NULL;

    for (int i = 0; i < rows; i++)
    {
        matrix[i] = new T [cols];

        // Safely unravel allocation if we run out of memory
        if (matrix[i] == NULL)
        {
            while (i--)
                delete [] matrix[i];

            delete [] matrix;

            return NULL;
        }
    }

    return matrix;
};

template <class T> T** AllocateMatrix(int rows, int cols, T value)
{
    T ** matrix = AllocateMatrix<T>(rows, cols);

    if (matrix != NULL)
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                matrix[i][j] = value;

    return matrix;
};

template <class T> void FreeMatrix(T ** & matrix, int rows)
{
    if (matrix == NULL)
        return;

    for (int i = 0; i < rows; i++)
        delete [] matrix[i];

    delete [] matrix;

    matrix = NULL;
};

#endif

