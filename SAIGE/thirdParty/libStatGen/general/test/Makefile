
include ../../Makefiles/Makefile.include

.DEFAULT_GOAL := all

SUBDIRS=inputFileTest cigar string memoryMapArrayTest packedVectorTest referenceSequenceTest nonOverlapRegions baseUtilitiesTest trimSequence reusableVector

OPTFLAG?=-O0

.PHONY: $(SUBDIRS) all test clean

all: TARGET = all
test: TARGET = test
clean: TARGET = clean

all test: $(SUBDIRS)

clean: $(SUBDIRS)
	# clean current directory
	rm -f *.o *~
	rm -f phiX-bs.umfa phiX-cs.umfa gtest

$(SUBDIRS):
	@$(MAKE) OPTFLAG="$(OPTFLAG)" -C $@ $(TARGET)

######################################################################

GTEST_COMPILE_FLAG=-DGTEST_HAS_TR1_TUPLE=0 -D GTEST_HAS_RTTI=0

LIBRARY=$(REQ_LIBS)

.c.o :
	$(CXX) $(CFLAGS) $(GTEST_COMPILE_FLAG) -o $@ -c $*.c -DVERSION=\""$(VERSION)\"" -DDATE="\"${DATE}\"" -DNODE="\"${NODE}\"" -DUSER="\"${USER}\"" -I$(GTEST_INC)

.cpp.o : 
	$(CXX) $(CFLAGS) $(GTEST_COMPILE_FLAG) -o $@ -c $*.cpp -DVERSION="\"$(VERSION)\"" -DDATE="\"${DATE}\"" -DNODE="\"${NODE}\"" -DUSER="\"${USER}\"" -I$(GTEST_INC)

# Google test part
GTEST_ROOT=../../../gtest-1.5.0
GTEST_INC = $(GTEST_ROOT)/include
GTEST_SRC = $(GTEST_ROOT)/src/gtest-all.cc
GTEST_OBJ = $(GTEST_ROOT)/gtest-all.o

$(GTEST_OBJ): $(GTEST_SRC)
	g++ -c -o $@ $< -I$(GTEST_INC) -I$(GTEST_ROOT)

#TOOLTEST=$(TOOLBASE:=_test)
TOOLTEST=gtest TrimSequence_test GreedyTupleAligner_test STLUtilities_test GenomeSequence_test Chromosome_test
TOOLTESTHDR=$(TOOLTEST:=.h)
TOOLTESTSRC=$(TOOLTEST:=.cpp)
TOOLTESTOBJ=$(TOOLTESTSRC:.cpp=.o)
TOOLOBJ_NOMAIN=$(TOOLOBJ:Main.o=)

GenomeSequence_test.o: GenomeSequence_test.cpp ../GenomeSequence.h ../GenomeSequence.cpp
	$(CXX) $(CFLAGS) $(GTEST_COMPILE_FLAG) -o $@ -c $< -DVERSION="\"$(VERSION)\"" -DDATE="\"${DATE}\"" -DNODE="\"${NODE}\"" -DUSER="\"${USER}\"" -I$(GTEST_INC)

STLUtilities_test.o: STLUtilities_test.cpp ../STLUtilities.h ../STLUtilities.cpp
	$(CXX) $(CFLAGS) $(GTEST_COMPILE_FLAG) -std=c++0x -o $@ -c $< -DVERSION="\"$(VERSION)\"" -DDATE="\"${DATE}\"" -DNODE="\"${NODE}\"" -DUSER="\"${USER}\"" -I$(GTEST_INC)

TrimSequence_test.o: TrimSequence_test.cpp ../TrimSequence.h
	$(CXX) $(CFLAGS) $(GTEST_COMPILE_FLAG) -o $@ -c $< -DVERSION="\"$(VERSION)\"" -DDATE="\"${DATE}\"" -DNODE="\"${NODE}\"" -DUSER="\"${USER}\"" -I$(GTEST_INC)

GreedyTupleAligner_test.o: GreedyTupleAligner_test.cpp ../GreedyTupleAligner.h
	$(CXX) $(CFLAGS) $(GTEST_COMPILE_FLAG) -o $@ -c $< -DVERSION="\"$(VERSION)\"" -DDATE="\"${DATE}\"" -DNODE="\"${NODE}\"" -DUSER="\"${USER}\"" -I$(GTEST_INC)

Chromosome_test.o: Chromosome_test.cpp ../Chromosome.h
	$(CXX) $(CFLAGS) $(GTEST_COMPILE_FLAG) -o $@ -c $< -DVERSION="\"$(VERSION)\"" -DDATE="\"${DATE}\"" -DNODE="\"${NODE}\"" -DUSER="\"${USER}\"" -I$(GTEST_INC)

gtest_build: $(TOOLTESTOBJ) $(TOOLOBJ) $(GTEST_OBJ)
	g++  $(TOOLTESTOBJ) $(TOOLOBJ_NOMAIN) $(GTEST_OBJ) $(LIBRARY) -lm -lz -lssl -lpthread -lcrypto -o ./gtest ../../libStatGen_debug.a
gtest: gtest_build
	./gtest --gtest_color=yes

include $(MAKEFILES_PATH)/Makefile.footer
