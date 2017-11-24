EXE = samTest
TOOLBASE = WriteFiles ValidationTest ReadFiles BamIndexTest ModifyVar Modify SamFileTest TestValidate TestEquals TestFilter ShiftIndels TestPileup TestPosList TestCigarHelper TestSamRecordPool TestSamCoordOutput TestSamRecordHelper
SRCONLY = Main.cpp
ifeq ($(ZLIB_AVAIL), 0)
TEST_COMMAND = ./test.sh noZlib
else
TEST_COMMAND = ./test.sh
endif

include ../../Makefiles/Makefile.test
