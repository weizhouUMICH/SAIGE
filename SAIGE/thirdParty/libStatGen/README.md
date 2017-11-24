Dependencies
------------

On debian type systems (including Ubuntu), add the following packages if they are not already installed (or have your admin add them if you do not have permission):

    sudo apt-get install g++ libssl-dev zlib1g-dev

Building
--------

To compile, from the top level directory, type: `make`.

To compile with debug symbols, type: `make debug`.

To test (after compiling), from the top level directory, type: `make test`.

Under the main statgen repository, there are: 

- `bam` - library code for operating on bam files.
- `copyrights` - copyrights for the library and any code included with it.
- `fastq` - library code for operating on fastq files.
- `general` - library code for general operations
- `glf` - library code for operating on glf files.
- `include` - after compiling, the library headers are linked here
- `Makefiles` - directory containing Makefiles that are used in the library and can be used for developing programs using the library
- `samtools` - library code used from samtools

After Compiling: `libStatGen.a`, `libStatGen_debug.a`, and `libStatGen_profile.a` are created at the top level.

Makefiles
---------
`Makefiles/Makefile.include` should contain the definitions that you need for creating software using this library.

`Makefiles/Makefile.lib` and `Makefiles/Makefile.src` can be used as templates for creating Makefiles for new software.  If possible, just include them within your Makefile.  Just set the proper variables for your program in your Makefile first.  (both Makefiles automatically include `Makefile.include`)

A similar setup should be used for test code, by including `Makefiles/Makefile.test` and defining your test specific variables first.

Other Notes
-----------
* Typically the `.o` files are compiled into their own directory called `obj`.


Compile Issues/Troubleshooting
------------------------------
See <http://genome.sph.umich.edu/wiki/LibStatGen_Troubleshooting> for the latest troubleshooting information.

If you are compiling on OSX and you encounter errors like the following:

> GenomeSequence.cpp: In instantiation of 'std::basic_ostream<_CharT, _Traits>& std::operator<<(std::basic_ostream<_CharT, _Traits>&, const std::basic_string<_CharT, _Traits, _Alloc>&) [with _CharT = char, _Traits = std::char_traits<char>, _Alloc = std::allocator<char>]':

> GenomeSequence.cpp:161: instantiated from here

> GenomeSequence.cpp:161: error: explicit instantiation of 'std::basic_ostream<_CharT, _Traits>& std::operator<<(std::basic_ostream<_CharT, _Traits>&, const std::basic_string<_CharT, _Traits, _Alloc>&) [with _CharT = char, _Traits = std::char_traits<char>, _Alloc = std::allocator<char>]' but no definition available

Try compiling with: `make USER_COMPILE_VARS=-mmacosx-version-min=10.8`

