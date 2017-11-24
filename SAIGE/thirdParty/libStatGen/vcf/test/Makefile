EXE = vcfTest
TOOLBASE = VcfFileTest VcfHeaderTest
SRCONLY = Main.cpp
TEST_COMMAND = ./vcfTest && diff results/vcfHeader.vcf expected/vcfHeader.vcf && diff results/vcfHeaderAddedFirst.vcf expected/vcfHeader.vcf && diff results/vcfHeaderAddedLast.vcf expected/vcfHeader.vcf && diff results/vcfHeaderAddedMiddle.vcf expected/vcfHeader.vcf && diff results/vcfFile.vcf testFiles/vcfFile.vcf && diff results/vcfFileNoInfo.vcf expected/vcfFileNoInfo.vcf && diff results/vcfFileNoInfoBGZF.vcf expected/vcfFileNoInfoBGZF.vcf && diff results/vcfFileNoInfoKeepGT.vcf expected/vcfFileNoInfoKeepGT.vcf && diff results/vcfFileNoInfoKeepGQHQ.vcf expected/vcfFileNoInfoKeepGQHQ.vcf

include ../../Makefiles/Makefile.test
