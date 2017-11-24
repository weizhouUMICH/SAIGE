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

#ifndef _SIMPLESTATS_H_
#define _SIMPLESTATS_H_

#include <math.h>   // for sqrt
#include <iostream>

//
// see http://www.johndcook.com/standard_deviation.html
// or Donald Knuth's Art of Computer Programming, Vol 2, page 232, 3rd edition
//
class RunningStat
{
public:
    RunningStat() : m_n(0), m_oldM(0), m_newM(0), m_oldS(0), m_newS(0) {}

    void Clear()
    {
        m_n = 0;
    }

    void Push(double x)
    {
        m_n++;

        // See Knuth TAOCP vol 2, 3rd edition, page 232
        if (m_n == 1)
        {
            m_oldM = x;
            m_oldS = 0.0;
            m_newM = x;
            m_newS = 0.0;
        }
        else
        {
            m_newM = m_oldM + (x - m_oldM)/m_n;
            m_newS = m_oldS + (x - m_oldM)*(x - m_newM);

            // set up for next iteration
            m_oldM = m_newM;
            m_oldS = m_newS;
        }
    }

    int NumDataValues() const
    {
        return m_n;
    }

    double Mean() const
    {
        return (m_n > 0) ? m_newM : 0.0;
    }

    double Variance() const
    {
        return ((m_n > 1) ? m_newS/(m_n - 1) : 0.0);
    }

    double StandardDeviation() const
    {
        return sqrt(Variance());
    }

private:
    uint64_t m_n;
    double m_oldM, m_newM, m_oldS, m_newS;
};

//
// helpers for Tabulate template
//
inline bool operator == (RunningStat &r, int i)
{
    return r.NumDataValues() == i;
}

inline std::ostream &operator << (std::ostream &stream, RunningStat &s)
{
    stream  << "N: " << s.NumDataValues()
    << " Mean: " << s.Mean()
    << " Standard Deviation: " << s.StandardDeviation();

    return stream;
}


#endif
