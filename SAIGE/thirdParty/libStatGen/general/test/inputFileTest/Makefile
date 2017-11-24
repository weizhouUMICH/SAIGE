EXE = inputFileTest
TOOLBASE = InputFileTest

ifeq ($(ZLIB_AVAIL), 0)
TEST_COMMAND=	mkdir -p results && \
	./inputFileTest > results/results.log && \
	diff data/InputFileTest.txt results/bgzfFile.glf && \
	diff data/InputFileTest.txt results/defaultFile.glf && \
	diff data/InputFileTest.txt results/defaultFile.gz && \
	diff data/InputFileTest.txt results/gzipFile.glf && \
	diff data/InputFileTest.txt results/InputFileTest.glf && \
	diff data/InputFileTest.txt results/uncompressedFile.glf && \
	diff data/textFile.txt results/textFile.gz && \
	diff data/textFile.txt results/textFile1.gz && \
	diff results/results.log expected/resultsNoZlib.log
else
TEST_COMMAND=	mkdir -p results && \
	./inputFileTest > results/results.log && \
	diff data/InputFileTest.bam results/bgzfFile.glf && \
	diff data/InputFileTest.txt results/defaultFile.glf && \
	diff data/InputFileTest.gz results/defaultFile.gz && \
	diff data/InputFileTest.gz results/gzipFile.glf && \
	diff data/InputFileTest.txt results/InputFileTest.glf && \
	diff data/InputFileTest.txt results/uncompressedFile.glf && \
	diff data/textFile.gz results/textFile.gz && \
	diff data/textFile.gz results/textFile1.gz && \
	diff results/results.log expected/results.log
endif


include ../../../Makefiles/Makefile.test
