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

#ifndef __MEMORYMAPARRAY_H
#define __MEMORYMAPARRAY_H

#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _WIN32
#include <unistd.h> // for gethostname()
#endif

#include <string>
#include <sys/types.h>
#include <time.h>

// STL:
#include <ostream>
#include <sstream>

#include "Generic.h"
#include "MemoryMap.h"


//
// This file defines a template for generating memory map backed arrays
// of different types of values.
//
// The template assumes that the mapped files are broken into two parts,
// first, a header (MemoryMapArrayHeader), then followed by the data
// in the array.
//
// typedefs are used to declare various types of arrays beforehand,
// since there will be only a few.
//
// They are:
//      mmapArrayUint32_t;
//      mmapArrayBool_t;
//      mmapArray4Bit_t;
//
// XXX consider adding env("USER"), argv[0], date/time creation, etc.
//
class MemoryMapArrayHeader
{
public:
    void constructorClear()
    {
        memset(this, 0, sizeof(*this));
    }
    uint32_t     typeCookie;
    uint32_t     typeVersion;
    uint32_t     contentCookie;
    uint32_t     contentVersion;
    size_t       headerSize;

    // file generation info
    time_t      creationDate;
    char        creationUser[32];
    char        creationHost[32];
    char        application[32];
    // now describe the data:
    size_t      elementCount;
    void debugPrint(FILE *);
    size_t getHeaderSize(int i)
    {
        return sizeof(*this);
    }

    void setApplication(const char *s)
    {
        strncpy(application, s, sizeof(application)-1);
        application[sizeof(application)-1] = '\0';
    }
    void setCreationUser(const char *s)
    {
        strncpy(creationUser, s, sizeof(creationUser)-1);
        creationUser[sizeof(creationUser)-1] = '\0';
    }
    void setCreationHost(const char *s)
    {
        strncpy(creationHost, s, sizeof(creationHost)-1);
        creationHost[sizeof(creationHost)-1] = '\0';
    }
};

//
// stream output for header information
//
std::ostream &operator << (std::ostream &stream, MemoryMapArrayHeader &h);

//
// This class object represents the application specific information that doesn't
// fit in the general header above.  Since it is only allocated via an mmap operation,
// as part of the mapped file, the destructor must never be called.  The virtual
// destructor is declared to eliminate gcc warnings.
//
// For many arrays, this will be empty.
//
struct MemoryMapGenericHeader
{
protected:
    size_t  headerSize;     // set in ::create and ::open only
public:
    size_t getHeaderSize()
    {
        return headerSize;
    }
    // other stuff follows...
};

template <
class elementT,
typename indexT,
unsigned int cookieVal,
unsigned int versionVal,
elementT accessorFunc(char *base, indexT),
void setterFunc(char *base, indexT, elementT),
size_t elementCount2BytesFunc(indexT),
class arrayHeaderClass>
class MemoryMapArray : public MemoryMap
{
protected:
    arrayHeaderClass    *header;
    char                *data;
    std::string         errorStr;
public:
    void constructorClear()
    {
        header = NULL;
        data = NULL;
//      errorStr = "";
    }
    MemoryMapArray()
    {
        constructorClear();
    }
    ~MemoryMapArray()
    {
        if (data) close();
    }

    const std::string &getErrorString()
    {
        return errorStr;
    }

    arrayHeaderClass &getHeader()
    {
        return *header;
    }

    void setContentCookie(uint32_t c)
    {
        header->contentCookie = c;
    }
    void setContentVersion(uint32_t v)
    {
        header->contentVersion = v;
    }

    // accessing
    inline elementT operator[](indexT i)
    {
        return accessorFunc(data, i);
    }
    inline void set(indexT i, elementT v)
    {
        setterFunc(data, i, v);
    }



    /// Create a vector with elementCount memebers.
    //
    /// Does administrative setup of the header and populating this
    /// class members.  User will need to finish populating the
    /// contents of the metaData and data sections.
    ///
    /// If file==NULL, the underlying allocation is done via malloc(),
    /// so that the results of write access to this vecor are not
    /// saved in a file.
    ///
    /// If file!=NULL, a file will be created on disk, and all
    /// write accesses done via the method ::set will be persistent
    /// in that file.
    ///
    int create(const char *file, indexT elementCount, int optionalHeaderCount = 0)
    {
        size_t len = elementCount2BytesFunc(elementCount) +
                     header->getHeaderSize(optionalHeaderCount);
        int rc;
        rc = MemoryMap::create(file, len);
        if (rc)
        {
            std::ostringstream buf;
            buf << file << ": failed to create file";
            errorStr = buf.str();
            close();
            return rc;
        }
        header = (arrayHeaderClass *) MemoryMap::data;
        header->constructorClear();
        header->typeCookie = cookieVal;
        header->typeVersion = versionVal;
        header->headerSize = header->getHeaderSize(optionalHeaderCount);
        header->elementCount = elementCount;
        data = (char *)((char *) MemoryMap::data + header->headerSize);

        const char *env;
        char hostname[256];
        env = getenv("USER");
        if (env) header->setCreationUser(env);
        header->creationDate = time(NULL);
#if defined(_WIN32)
        hostname[0] = '\0';
#else
        gethostname(hostname, sizeof(hostname));
#endif
        header->setCreationHost(hostname);
        return 0;
    }

    /// allow anonymous (malloc) create.
    ///
    /// we do this when we don't expect to save the results.
    ///
    /// The single use case so far is in GenomeSequence::populateDBSNP.
    ///
    int create(indexT elementCount, int optionalHeaderCount = 0)
    {
        return create(NULL, elementCount, optionalHeaderCount);
    }

    //
    // Open the given filename.  flags may be set to
    // O_RDONLY or O_RDWR, and allows the file to be
    // condtionally written to.
    //
    // Several sanity checks are done:
    //   compare the expected cookie value to the actual one
    //   compare the expected version value to the actual one
    //
    // if either condition is not met, the member errorStr is
    // set to explain why, and true is returned.
    //
    // If there were no errors, false is returned.
    //
    bool open(const char *file, int flags = O_RDONLY)
    {
        int rc = MemoryMap::open(file, flags);
        if (rc)
        {
            std::ostringstream buf;
            buf << file << ": open() failed (error=" << strerror(errno) << ").";
            errorStr = buf.str();
            return true;
        }
        header = (arrayHeaderClass *) MemoryMap::data;
        data = (char *)((char *) MemoryMap::data + header->headerSize);
        if (header->typeCookie!=cookieVal)
        {
            std::ostringstream buf;
            buf << file << ": wrong type of file (expected type "
            << cookieVal << " but got " << header->typeCookie << ")";
            errorStr = buf.str();
            // XXX insert better error handling
            close();
            return true;
        }
        if (header->typeVersion!=versionVal)
        {
            std::ostringstream buf;
            buf << file << ": wrong version of file (expected version "
            << versionVal << " but got " << header->typeVersion << ")";
            errorStr = buf.str();
            // XXX insert better error handling
            close();
            return true;
        }
        return false;
    }

    bool close()
    {
        constructorClear();
        return MemoryMap::close();
    }
    void debugPrint(FILE *f)
    {
        if (header) header->debugPrint(f);
    }

    size_t getElementCount() const
    {
        return header->elementCount;
    }

};

struct emptyGenericHeader : public MemoryMapGenericHeader
{
public:
    size_t getHeaderSize()
    {
        return sizeof(*this);
    }
};

//
// define the uint32 array type:
//
inline uint32_t mmapUint32Access(char *base, uint32_t index)
{
    return ((uint32_t *)base)[index];
}
inline void mmapUint32Set(char *base, uint32_t index, uint32_t v)
{
    ((uint32_t *)base)[index] = v;
}
inline size_t mmapUint32elementCount2Bytes(uint32_t i)
{
    return sizeof(uint32_t) * i;
}

typedef MemoryMapArray<
uint32_t,
uint32_t,
0x16b3816c,
20090109,
mmapUint32Access,
mmapUint32Set,
mmapUint32elementCount2Bytes,
MemoryMapArrayHeader
> mmapArrayUint32_t;

//
// define the boolean memory mapped array type.
// NB: it is limited to 2**32 elements
//

typedef MemoryMapArray<
uint32_t,
uint32_t,
0xac6c1dc7,
20090109,
PackedAccess_1Bit,
PackedAssign_1Bit,
Packed1BitElementCount2Bytes,
MemoryMapArrayHeader
> mmapArrayBool_t;

//
// define the two bit memory mapped array type:
//

typedef MemoryMapArray<
uint32_t,
uint32_t,
0x25b3ea5f,
20090109,
PackedAccess_2Bit,
PackedAssign_2Bit,
Packed2BitElementCount2Bytes,
MemoryMapArrayHeader
> mmapArray2Bit_t;

typedef MemoryMapArray<
uint32_t,
uint32_t,
0x418e1874,
20090109,
PackedAccess_4Bit,
PackedAssign_4Bit,
Packed4BitElementCount2Bytes,
MemoryMapArrayHeader
> mmapArray4Bit_t;

#if 0
// XXX this is example code I want to use to define arrays of genome wide match values
class   baseRecord
{
    unsigned int base:4;
    unsigned int qScore:7;
    unsigned int conflicts:5;   // how many cases of poorer matches that disagree
};

//
// define the baseRecord array type:
//
inline baseRecord &mmapBaseRecordAccess(void *base, uint32_t index)
{
    return *((baseRecord *)((char *)base + index*sizeof(baseRecord)));
}
inline void mmapBaseRecordSet(void *base, uint32_t index, baseRecord &v)
{
    mmapBaseRecordAccess(base, index) = v;
}
inline size_t mmapBaseRecordElementCount2Bytes(uint32_t i)
{
    return sizeof(baseRecord) * i;
}

typedef MemoryMapArray<
baseRecord &,
uint32_t,
0x12341234,
0xdeadbeef,
&mmapBaseRecordAccess,
mmapBaseRecordSet,
mmapBaseRecordElementCount2Bytes,
MemoryMapArrayHeader
> mmapArrayBaseRecord_t;
#endif

#endif
