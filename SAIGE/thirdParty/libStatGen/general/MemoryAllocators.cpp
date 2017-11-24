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

#include "MemoryAllocators.h"

#include <stdlib.h>

char *** AllocateCharCube(int n, int rows, int cols)
{
    char *** cube = new char ** [n];

    // Stop early if we are out of memory
    if (cube == NULL)
        return NULL;

    for (int i = 0; i < n; i++)
    {
        cube[i] = AllocateCharMatrix(rows, cols);

        // Safely unravel allocation if we run out of memory
        if (cube[i] == NULL)
        {
            while (i--)
                FreeCharMatrix(cube[i], rows);

            delete [] cube;

            return NULL;
        }
    }

    return cube;
}

int ** AllocateIntMatrix(int rows, int cols)
{
    int ** matrix = new int * [rows];

    // Stop early if we are out of memory
    if (matrix == NULL)
        return NULL;

    for (int i = 0; i < rows; i++)
    {
        matrix[i] = new int [cols];

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
}

char ** AllocateCharMatrix(int rows, int cols)
{
    char ** matrix = new char * [rows];

    // Stop early if we are out of memory
    if (matrix == NULL)
        return NULL;

    for (int i = 0; i < rows; i++)
    {
        matrix[i] = new char [cols];

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
}

float ** AllocateFloatMatrix(int rows, int cols)
{
    float ** matrix = new float * [rows];

    // Stop early if we are out of memory
    if (matrix == NULL)
        return NULL;

    for (int i = 0; i < rows; i++)
    {
        matrix[i] = new float [cols];

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
}

void FreeCharCube(char *** & cube, int n, int rows)
{
    if (cube == NULL)
        return;

    for (int i = 0; i < n; i++)
        FreeCharMatrix(cube[i], rows);

    delete [] cube;

    cube = NULL;
}

void FreeCharMatrix(char ** & matrix, int rows)
{
    if (matrix == NULL)
        return;

    for (int i = 0; i < rows; i++)
        delete [] matrix[i];

    delete [] matrix;

    matrix = NULL;
}

void FreeFloatMatrix(float ** & matrix, int rows)
{
    if (matrix == NULL)
        return;

    for (int i = 0; i < rows; i++)
        delete [] matrix[i];

    delete [] matrix;

    matrix = NULL;
}

void FreeIntMatrix(int ** & matrix, int rows)
{
    if (matrix == NULL)
        return;

    for (int i = 0; i < rows; i++)
        delete [] matrix[i];

    delete [] matrix;

    matrix = NULL;
}

short ** AllocateShortMatrix(int rows, int cols)
{
    short ** matrix = new short * [rows];

    // Stop early if we are out of memory
    if (matrix == NULL)
        return NULL;

    for (int i = 0; i < rows; i++)
    {
        matrix[i] = new short [cols];

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
}

void FreeShortMatrix(short ** & matrix, int rows)
{
    if (matrix == NULL)
        return;

    for (int i = 0; i < rows; i++)
        delete [] matrix[i];

    delete [] matrix;

    matrix = NULL;
}

double ** AllocateDoubleMatrix(int rows, int cols)
{
    double ** matrix = new double * [rows];

    // Stop early if we are out of memory
    if (matrix == NULL)
        return NULL;

    for (int i = 0; i < rows; i++)
    {
        matrix[i] = new double [cols];

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
}

void FreeDoubleMatrix(double ** & matrix, int rows)
{
    for (int i = 0; i < rows; i++)
        delete [] matrix[i];

    delete [] matrix;

    matrix = NULL;
}


