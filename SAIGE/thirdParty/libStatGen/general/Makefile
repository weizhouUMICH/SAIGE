UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
  # For Linux, add a couple extra warnings
  #
  # trying to strike the right balance on error/warning
  # handling is difficult, but here we make most everthing
  # an error, but allow unused results and variables for the
  # time being.
  #
  USER_WARNINGS ?= -Werror $(shell if [ X$(CCVERSION) \> X4.2.0 ] ; then echo " -Wno-strict-overflow" ; fi)
#-Wno-strict-overflow
# -Wno-unused-variable $(shell if [ X$(CCVERSION) \> X4.2.0 ] ; then echo " -Wno-unused-result" ; fi)
endif
ifeq ($(UNAME), Darwin)
  # Mac OS, nothing to add.
  # -Wall contain -Wunused-variable and -Wunused-result from 'man gcc'
endif

#USER_COMPILE_VARS = -D_NO_PHONEHOME

TOOLBASE=\
	BaseAsciiMap \
	BaseQualityHelper \
	BaseUtilities \
	BasicHash \
	BgzfFileType \
	BgzfFileTypeRecovery \
	CharBuffer \
	Chromosome \
	Cigar \
	CigarRoller \
	Error \
	ErrorHandler \
	FileType \
	FortranFormat \
	GenomeSequence \
	GenotypeLists \
	glfHandler \
	GzipFileType \
	GzipHeader \
	Hash \
	IndexBase \
	Input \
	InputFile \
	IntArray \
	IntHash \
	LongLongCounter \
	MapFunction \
	MathMatrix \
	MathVector \
	MemoryAllocators \
	MemoryInfo \
	MemoryMapArray \
	MemoryMap \
	MiniDeflate \
	NonOverlapRegions \
	Parameters \
	PedigreeAlleleFreq \
	Pedigree \
	PedigreeDescription \
	PedigreeFamily \
	PedigreeGlobals \
	PedigreePerson \
	PhoneHome \
	QuickIndex \
	Random \
	ReferenceSequence \
	SmithWaterman \
	Sort \
	STLUtilities \
	StatGenStatus \
	StringAlias \
	StringArray \
	StringBasics \
	StringHash \
	StringMap \
	Tabix \
	UncompressedFileType

SRCONLY=\
	PedigreeLoader.cpp \
	PedigreeTrim.cpp \
	PedigreeTwin.cpp 

HDRONLY= \
	Constant.h \
	CSG_MD5.h \
	Generic.h \
	GenomeSequenceHelpers.h \
	GreedyTupleAligner.h \
	InplaceMerge.h \
	LongHash.h \
	LongInt.h \
	MathConstant.h \
	PackedVector.h \
	PedigreeAlleles.h \
	Performance.h \
	ReusableVector.h \
	SimpleStats.h \
	TrimSequence.h \
	UnitTest.h

include ../Makefiles/Makefile.lib


CPP_TESTS := $(shell grep -l 'if defined(TEST)' *cpp)
selftest: 
	@for i in $(CPP_TESTS); do \
		if [ "XXX$$i" = XXX ] ;\
		then \
		    continue; \
		fi;\
		(echo "building in self test $$i"; g++ -DTEST -D__STDC_LIMIT_MACROS -o test_$$i $$i  -L. -lcsg -lm -lz -lssl ; ./test_$$i) ; \
		if [ $$? -ne 0 ] ; \
		then \
		    echo "make stopped because of errors." ; \
		    break ; \
		fi \
	done
