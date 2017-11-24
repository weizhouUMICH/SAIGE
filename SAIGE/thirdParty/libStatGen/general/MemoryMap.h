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

#ifndef __MEMORYMAP_H
#define __MEMORYMAP_H
#include <sys/types.h>
#include <fcntl.h>

#if defined(_WIN32)
#include <windows.h>
#endif

///
/// There are a pair of related data structures in the operating system,
/// and also a few simple algorithms that explain why your processes are
/// waiting forever.
///
/// The symptom you have is that they are getting little or no CPU time,
/// as shown in the command 'top'.  The machine will appear to have
/// available CPU time (look at the Cpu(s): parameter - if less than 100%,
/// you have available CPU).  The real key, however, is to look at the
/// 'top' column with the label 'S' - that is the status of the process,
/// and crucial to understanding what is going on.
///
/// In your instance, the 'S' column for your karma jobs is 'D', which
/// means it is waiting for data.  This is because the process is doing
/// something that is waiting for the filesystem to return data to it.
/// Usually, this is because of a C call like read() or write(), but it
/// also happens in large processes where memory was copied to disk and
/// re-used for other purposes (this is called paging).
///
/// So, a bit of background on the operating system... there is a CPU
/// secheduler that takes a list of waiting processes, and picks one to
/// run - if the job is waiting for the disk, there is no point in picking
/// it to run, since it is blocked, waiting for the disk to return data.
/// The scheduler marks the process with 'D' and moves on to the next
/// process to schedule.
///
/// In terms of data structures that we care about for this example, there
/// are two that we care about.  First is a linear list of disk buffers
/// that are stored in RAM and controlled by the operating system.  This
/// is usually called the disk buffer pool.  Usually, when a program asks
/// for data from the disk, this list can be scanned quickly to see if the
/// data is already in RAM - if so, no disk operation needs to take place.
///
/// Now in the case of the normal Unix read() and write() calls, when the
/// operating system is done finding the page, it copies the data into a
/// buffer to be used by the process that requested it (in the case of a
/// read() - a write() is the opposite).  This copy operation is slow and
/// inefficient, but gets the job done.
///
/// So overall, you gain some efficiency in a large memory system by
/// having this disk buffer pool data structure, since you aren't
/// re-reading the disk over and over to get the same data that you
/// already have in RAM.  However, it is less efficient than it might be
/// because of the extra buffer copying.
///
/// Now we come to memory mapped files, and karma.  The underlying system
/// call of interest to us is mmap(), and is in MemoryMap.cpp.  What it
/// does and how it works are important to understanding the benefits of
/// it, and frankly, most people don't care about it because it is
/// seemingly complex.
///
/// Two things are important to know: firstly, there is a data structure
/// in the CPU called the page table, which is mostly contained in the CPU
/// hardware itself.  All memory accesses for normal user processes like
/// karma go through this hardware page table.  Secondly, it is very fast
/// for the operating system to put together a page table that 'connects'
/// a bunch of memory locations in your user programs address space to the
/// disk buffer pool pages.
///
/// The combination of those two facts mean that you can implement a 'zero
/// copy' approach to reading data, which means that the data that is in
/// the disk buffer pool is directly readable by the program without the
/// operating system ever having to actually copy the data, like it does
/// for read() or write().
///
/// So the benefit of mmap() is that when the underlying disk pages are
/// already in the disk buffer pool, a hardware data structure gets built,
/// then the program returns, and the data is available at full processor
/// speed with no intervening copy of the data, or waiting for disk or
/// anything else.  It is as near to instantaneous as you can possibly
/// get.  This works whether it is 100 bytes or 100 gigabytes.
///
/// So, the last part of the puzzle is why your program winds up in 'D'
/// (data wait), and what to do about it.
///
/// The disk buffer pool is a linear list of blocks ordered by the time
/// and date of access.  A process runs every once in awhile to take the
/// oldest of those pages, and free them, during which it also has to
/// update the hardware page tables of any processes referencing them.
///
/// So on wonderland, most file access (wget, copy, md5sum, anything else)
/// is constantly putting new fresh pages at the front of the list, and
/// karma index files, having been opened awhile ago, are prime candidates
/// for being paged out.  The reason they get paged out as far as I know
/// is that in any given second of execution, nowhere near the entire
/// index is getting accessed... so at some point, at least one page gets
/// sent back to disk (well, flushed from RAM).  Once that happens, a
/// cascading effect happens, where the longer it waits, the older the
/// other pages get, then the more that get reclaimed, and the slower it
/// gets, until karma is at a standstill, waiting for pages to be brought
/// back into RAM.
///
/// Now in an ideal world, karma would rapidly recover, and it can...
/// sometimes.  The problem is that your karma job is accessing data all
/// over that index, and it is essentially looking like a pure random I/O
/// to the underlying filesystem.  There is about a 10 to 1 performance
/// difference between accessing the disk sequentially as compared to
/// randomly.
///
/// So to make karma work better, the first thing I do when starting karma
/// is force it to read all of the disk pages in order.  This causes the
/// entire index to be forced into memory in order, so it is forcing
/// sequential reads, which is the best case possible.  There are
/// problems, for example if three karma jobs start at once, the disk I/O
/// is no longer as purely sequential as we would like.  Also, if the
/// filesystem is busy taking care of other programs, even if karma thinks
/// it is forcing sequential I/O, the net result looks more random.  This
/// happens when the system is starting to break down (thrashing) and it
/// will certainly stall, or look very very slow, or crash.
///
/// The upshot of all of this is that when a single reference is shared,
/// it is more likely that all the pages will be in the disk buffer pool
/// to begin with, and thereby reduce startup time to nearly zero.  It is
/// also the ideal situation in terms of sharing the same reference among
/// say 24 copies of karma on wonderland - the only cost is the hardware
/// page table that gets set up to point to all of the disk buffers.
///
/// As I mentioned a paragraph back, the pages can still get swapped out,
/// even with dozens of karma jobs running.  A workaround I created is a
/// program in utilities called mapfile - it simply repeatedly accesses
/// the data in sequential order to help ensure that all of the pages are
/// at the head of the disk buffer pool, and therefore less likely to get
/// swapped out.
///
/// The benefit of such a program (mapfile) is greater on wonderland,
/// where a lot of processes are competing for memory and disk buffers.
///
///
class MemoryMap
{
#if defined(_WIN32)
    HANDLE      file_handle;
    HANDLE      map_handle;
    DWORD       page_size;
#else
    int fd;
    size_t page_size;
#endif
    off_t offset;
    size_t mapped_length;
    size_t total_length;
    bool    useMemoryMapFlag;
public:

    void *data;

    MemoryMap();

    virtual ~MemoryMap();

    void debug_print();

    void constructor_clear();

    void destructor_clear();

    virtual bool allocate();

    /// open a previously created mapped vector
    ///
    /// useMemoryMapFlag will determine whether it
    /// uses mmap() or malloc()/read() to populate
    /// the memory
    virtual bool open(const char * file, int flags = O_RDONLY);

    /// create the memory mapped file on disk
    ///
    /// a file will be created on disk with the header
    /// filled in.  The caller must now populate elements
    /// using (*this).set(index, value).
    //
    virtual bool create(const char * file, size_t size);

    /// store in allocated memory (malloc), not mmap:
    ///
    /// This is for code that needs to more flexibly
    /// the case when an mmap() file _might_ be available,
    /// but if it is not, we want to load it as a convenience
    /// to the user. GenomeSequence::populateDBSNP does exactly this.
    //
    virtual bool create(size_t size);

    bool close();
    void test();
    size_t length()
    {
        return mapped_length;
    }

    char operator[](unsigned int index)
    {
        return ((char *)data)[index];
    };
    int prefetch();     // force pages into RAM

    //
    // set or unset use of mmap() call in ::open().
    // This flag must be set before ::open() is called,
    // if it is called afterwards, it has no effect.
    //
    void useMemoryMap(bool flag=true)
    {
        useMemoryMapFlag = flag;
    }
};

#endif
