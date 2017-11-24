/*
 *  Copyright (C) 2011  Regents of the University of Michigan,
 *                      Hyun Min Kang, Matthew Flickenger, Matthew Snyder,
 *                      and Goncalo Abecasis
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

#include <vector>
#include <stdexcept>

#ifndef __REUSABLE_VECTOR_H__
#define __REUSABLE_VECTOR_H__

/// Create a vector of DATA_TYPE that reuses created objects to save on
/// memory reallocations.  DATA_TYPE must have a function called clear()
/// that is used to reset it for reuse.
template <class DATA_TYPE>
class ReusableVector
{
public:
    ReusableVector(): myCont(), myNextEmpty(0) {}
    virtual ~ReusableVector();

    /// Clear the vector contents.
    void reset();
    /// Clear the vector contents.
    void clear() {reset();}

    /// Get a reference to a new entry to be populated so the user can
    /// directly populate it rather than having to copy into it.
    DATA_TYPE& getNextEmpty();

    /// Get a reference to the data at the specified index.
    /// Throws an exception if the index is out of range.
    DATA_TYPE& get(unsigned int index) const;

    /// Return the number of populated entries in the vector.
    // The next empty position is the same as the size.
    int size() const {return(myNextEmpty);}

    void rmLast();

protected:
    std::vector<DATA_TYPE*> myCont;
    unsigned int myNextEmpty;
private:
    ReusableVector& operator=(const ReusableVector& rv);
    ReusableVector(const ReusableVector& rv);
};


/////////////////////////////////////////////////////////////
// ReusableVector
template <class DATA_TYPE>
ReusableVector<DATA_TYPE>::~ReusableVector()
{
    for(unsigned int i = 0; i < myCont.size(); i++)
    {
        // Delete all the entries.
        delete myCont[i];
        myCont[i] = NULL;
    }
    myCont.clear();
    myNextEmpty = 0;
}


template <class DATA_TYPE>
void ReusableVector<DATA_TYPE>::reset()
{
    // Set the next empty element to be the first one on the list.
    // That means there are none used.
    myNextEmpty = 0;
}
        

template <class DATA_TYPE>
DATA_TYPE& ReusableVector<DATA_TYPE>::getNextEmpty()
{
    if(myNextEmpty == myCont.size())
    {
        // We are at the end of the available entries, so add a new one.
        myCont.resize(myCont.size() + 1);

        // Create a new entry.
        myCont[myNextEmpty] = new DATA_TYPE;
    }
    else
    {
        // myNextEmpty is an element, and not the end.
        // So, clear out the data.
        myCont[myNextEmpty]->clear();
    }

    DATA_TYPE* returnVal = myCont[myNextEmpty];

    // Increment next empty to the next element.
    ++myNextEmpty;
    // return the element to be used.
    return(*returnVal);
}


template <class DATA_TYPE>
DATA_TYPE& ReusableVector<DATA_TYPE>::get(unsigned int index) const
{
    if((index < myNextEmpty) && (index >= 0))
    {
        // index is a valid position, so return that data.
        if(myCont[index] == NULL)
        {
            throw(std::runtime_error("ReusableVector::get BUG, found a null pointer."));
        }
        return(*myCont[index]);
    }

    // Not set in the vector, so throw an exception.
    throw(std::runtime_error("ReusableVector::get called with out of range index."));
    // return(myCont[0]);
}

template <class DATA_TYPE>
void ReusableVector<DATA_TYPE>::rmLast()
{
    if(myNextEmpty > 0)
    {
        --myNextEmpty;
    }
}

#endif
