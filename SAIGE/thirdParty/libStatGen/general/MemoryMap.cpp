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

#include <sys/stat.h>
#include <sys/types.h>
#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <iostream>
#include <stdio.h>
#include <stdexcept>
#include <stdlib.h>
#include <string>

#include "MemoryMap.h"

#ifndef _WIN32
#include <sys/mman.h>
#include <unistd.h>
#endif

#ifndef MAP_POPULATE
#define MAP_POPULATE 0x0000
#endif
#ifndef MAP_NONBLOCK
#define MAP_NONBLOCK 0x0000
#endif

MemoryMap::MemoryMap()
{
    constructor_clear();
#if defined(_WIN32)
    SYSTEM_INFO sysinfo = {0};
    ::GetSystemInfo(&sysinfo);
    DWORD cbView = sysinfo.dwAllocationGranularity;
#else
    page_size = sysconf(_SC_PAGE_SIZE);
#endif
}

MemoryMap::~MemoryMap()
{
    destructor_clear();
};

void MemoryMap::debug_print()
{
#if defined(_WIN32)
    std::cout << "fd = " << file_handle << std::endl;
#else
    std::cout << "fd = " << fd << std::endl;
#endif
    std::cout << "data = 0x" << std::hex << data << std::endl;
    std::cout << "offset = 0x" << std::hex << offset << std::endl;
    std::cout << "mapped_length = 0x" << std::hex << mapped_length << std::endl;
    std::cout << "total_length = 0x" << std::hex << total_length << std::endl;
    std::cout << "page_size = 0x" << std::hex << page_size << std::endl;
};

void MemoryMap::constructor_clear()
{
#if defined(_WIN32)
    file_handle = NULL;
    map_handle = NULL;
#else
    fd = -1;
#endif
    data = (void *) NULL;
    offset = 0;
    mapped_length = 0;
    total_length = 0;
    useMemoryMapFlag = true;
};

void MemoryMap::destructor_clear()
{
#if defined(_WIN32)
    if (data!=NULL)
    {
        // free windows mapped object
        ::UnmapViewOfFile((LPVOID) data);
    }
    if (map_handle != NULL)
       ::CloseHandle(map_handle);
    if (file_handle != NULL)
       ::CloseHandle(file_handle);
#else
    if (data!=NULL)
    {
        // free unix mapped object
        munmap(data, mapped_length);
    }
    // free unix resources
    if (fd!=-1)
    {
        ::close(fd);
    }
#endif

    constructor_clear();
}


bool MemoryMap::allocate()
{
    data = (void *) malloc(mapped_length);

    if (data == NULL)
    {
#ifdef __WIN32__
        ::CloseHandle(file_handle);
#else
        ::close(fd);
#endif
        perror("MemoryMap::open");
        constructor_clear();
        return true;
    }

#ifdef __WIN32__
    DWORD resultSize = 0;
    ReadFile(file_handle, data, mapped_length, &resultSize, NULL);
#else
    size_t resultSize = read(fd, data, mapped_length);
#endif

    if ( resultSize != mapped_length)
    {
#ifdef __WIN32__
        ::CloseHandle(file_handle);
#else
        ::close(fd);
#endif
        perror("MemoryMap::open");
        constructor_clear();
        return true;
    }
    return false;
}


bool MemoryMap::open(const char * file, int flags)
{
   const char * message = "MemoryMap::open - problem opening file %s";
#if defined(_WIN32)
    file_handle = CreateFile(file,
                             (flags==O_RDONLY) ? GENERIC_READ : (GENERIC_READ | GENERIC_WRITE),
                             FILE_SHARE_READ | FILE_SHARE_WRITE, // subsequent opens may either read or write
                             NULL,
                             OPEN_EXISTING,
                             FILE_ATTRIBUTE_NORMAL,
                             NULL);

    if(file_handle == INVALID_HANDLE_VALUE)
    {
        fprintf(stderr, message, file);
        constructor_clear();
        return true;
    }

    LARGE_INTEGER file_size = {0};
    ::GetFileSizeEx(file_handle, &file_size);
    mapped_length = total_length = file_size.QuadPart;

#else
    struct stat buf;
    fd = ::open(file, flags);
    if ((fd==-1) || (fstat(fd, &buf) != 0))
    {
        fprintf(stderr, message, file);
        constructor_clear();
        return true;
    }
    mapped_length = total_length = buf.st_size;
#endif

    if(!useMemoryMapFlag)
    {
        return allocate();
    }

#if defined(_WIN32)
    assert(offset == 0);

    map_handle = CreateFileMapping(file_handle, NULL,
                                   (flags==O_RDONLY) ? PAGE_READONLY : PAGE_READWRITE,
                                   file_size.HighPart, // upper 32 bits of map size
                                   file_size.LowPart,  // lower 32 bits of map size
                                   NULL);

    if(map_handle == NULL)
    {
        ::CloseHandle(file_handle);
        fprintf(stderr, message, file);
        constructor_clear();
        return true;
    }

    data = MapViewOfFile(map_handle,
                         (flags == O_RDONLY) ? FILE_MAP_READ : FILE_MAP_ALL_ACCESS,
                         0, 0, mapped_length);

    if (data == NULL)
    {
        CloseHandle(map_handle);
        CloseHandle(file_handle);

        fprintf(stderr, message, file);
        constructor_clear();
        return true;
    }
#else
    data = ::mmap(NULL, mapped_length,
                  (flags == O_RDONLY) ? PROT_READ : PROT_READ | PROT_WRITE,
                  MAP_SHARED, fd, offset);

    if (data == MAP_FAILED)
    {
        ::close(fd);
        fprintf(stderr, message, file);
        constructor_clear();
        return true;
    }
#endif
    return false;
}


bool MemoryMap::create(const char *file, size_t size)
{
    if (file==NULL)
    {
        data = calloc(size, 1);
        return(data==NULL);
    }

    const char * message = "MemoryMap::create - problem creating file %s";

#ifdef __WIN32__
    file_handle = CreateFile(file,
                             GENERIC_READ | GENERIC_WRITE,
                             FILE_SHARE_READ | FILE_SHARE_WRITE,
                             NULL,
                             CREATE_ALWAYS,
                             FILE_ATTRIBUTE_NORMAL,
                             NULL);

    if (file_handle == INVALID_HANDLE_VALUE)
    {
        fprintf(stderr, message, file);
        constructor_clear();
        return true;
    }

    SetFilePointer(file_handle, size - 1, NULL, FILE_BEGIN);
    char dummy = 0;
    DWORD check = 0;
    WriteFile(file_handle, &dummy, 1, &check, NULL);

    if (check != 0)
    {
        CloseHandle(file_handle);
        DeleteFile(file);
        fprintf(stderr, message, file);
        constructor_clear();
        return true;
    }
    CloseHandle(file_handle);
    open(file, O_RDWR);
#else
    fd = ::open(file, O_RDWR|O_CREAT|O_TRUNC, 0666);
    if(fd == -1)
    {
        fprintf(stderr, message, file);
        constructor_clear();
        return true;
    }

    lseek(fd, (off_t) size - 1, SEEK_SET);
    char dummy = 0;
    if(write(fd, &dummy, 1)!=1)
    {
        fprintf(stderr, message, file);
        constructor_clear();
        return true;
    }

    data = ::mmap(NULL, size, PROT_READ|PROT_WRITE,
                  MAP_SHARED, fd, offset);

    if (data == MAP_FAILED)
    {
        ::close(fd);
        unlink(file);
        fprintf(stderr, message, file);
        constructor_clear();
        return true;
    }
   mapped_length = total_length = size;
#endif
    return false;
}


bool MemoryMap::create(size_t size)
{
    return create(NULL, size);
}

bool MemoryMap::close()
{
    destructor_clear();
    return false;
}

void MemoryMap::test()
{
    int result;

    result = this->open("test/test_memmap_data.txt");
    assert(result == 0);
    assert(data!=NULL);
    assert(mapped_length == 183);   // length of the above file
    close();

    // now try non memory mapped (direct slow file I/O)
    useMemoryMap(false);
    result = this->open("test/test_memmap_data.txt");
    assert(result == 0);
    assert(data!=NULL);
    assert(mapped_length == 183);   // length of the above file
    close();
}

int MemoryMap::prefetch()
{
    int sum = 0;
    size_t i;

    for (i=0; i<mapped_length; i += page_size) sum += *(i + (char *) data);

    return sum;
}

#if defined(TEST)
//
// compile test using:
// g++ -DTEST -o testMemoryMap MemoryMap.cpp Error.o -lz
//

int main(int argc, const char *argv)
{
    MemoryMap map;

    map.test();

//  map.debug_print();

    exit(0);
}
#endif

