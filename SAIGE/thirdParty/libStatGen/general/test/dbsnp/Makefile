PATH_TO_BASE=../../..
EXE = dbsnpTest
SRCONLY = Main.cpp

TEST_COMMAND=@mkdir -p results; ./dbsnpTest 2> results/results.txt && diff expected/results.txt results/results.txt

include ../../../Makefiles/Makefile.test
