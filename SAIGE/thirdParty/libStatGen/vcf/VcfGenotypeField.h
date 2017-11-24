/*
 *  Copyright (C) 2012  Regents of the University of Michigan
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


#ifndef __VCF_GENOTYPE_FIELD_H__
#define __VCF_GENOTYPE_FIELD_H__

#include "VcfRecordField.h"
#include "ReusableVector.h"

class VcfGenotypeField
{
public:
    static const int GENOTYPE_INDEX_NA = -1;

    /// Default Constructor, initializes the variables.
    VcfGenotypeField()
        : myGenotypeSubFields()
    {}

    /// Destructor
    virtual ~VcfGenotypeField()
    {}

    bool write(IFILE filePtr);

    inline void reset() { myGenotypeSubFields.reset(); internal_reset(); }
    void clear() {reset(); }

    /// Get the number of genotype format fields there are.
    /// \return the number of genotype format fields.
    inline int getNumFields() { return(myGenotypeSubFields.size()); }

protected:
    enum SUBFIELD_READ_STATUS {
        MORE_SUBFIELDS = 0,
        END_OF_FIELD = 1,
        END_OF_RECORD = 2,
    };

    // Specify null if the field is not to be stored.
    static SUBFIELD_READ_STATUS readGenotypeSubField(IFILE filePtr, 
                                                     std::string* stringDest);
    // The stopCharacters must be in the same order as, "\n\t:".
    static SUBFIELD_READ_STATUS getReadStatus(int stopChar);

    virtual void internal_reset() = 0;

    ReusableVector<std::string> myGenotypeSubFields;
private:
    VcfGenotypeField(const VcfGenotypeField& vcfGenotypeField);
    VcfGenotypeField& operator=(const VcfGenotypeField& vcfGenotypeField);
};

#endif
