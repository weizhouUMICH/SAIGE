EXE = printRefPositions

SRCONLY = Main.cpp
TOOLBASE = PrintRefPositions

TEST_COMMAND = mkdir -p results; ./printRefPositions > results/test.txt; diff results/test.txt expected/test.txt

include ../../../Makefiles/Makefile.test
