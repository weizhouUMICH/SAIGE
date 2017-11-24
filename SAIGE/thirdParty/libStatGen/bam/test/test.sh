expected="expected/"
myDiffs=""
if [ $# -eq 1 ]
then
    expected="$expected$1"
./samTest 2> results/samTest.log && \
diff $expected/testEqWithBases.sam results/outSamEqBases.sam && \
diff $expected/testEqWithEq.sam results/outSamEqEquals.sam && \
diff $expected/testEqWithOrig.sam results/outSamEqOrig.sam && \
diff $expected/testEqWithBases.bam results/outSamEqBases.bam && \
diff $expected/testEqWithEq.bam results/outSamEqEquals.bam && \
diff $expected/testEqWithOrig.bam results/outSamEqOrig.bam && \
diff $expected/updateTagFromSam.sam results/updateTagFromSam.sam && \
diff $expected/updateTag.bam results/updateTagFromSam.bam && \
diff $expected/samTest.log results/samTest.log && \
diff $expected/addedTagToSam.bam results/addedTagToSam.bam && \
diff $expected/addedTagToSam.sam results/addedTagToSam.sam && \
diff $expected/testShift.sam results/testShift.sam && \
diff $expected/testShift.bam results/testShiftFromSam.bam && \
diff $expected/TestSamCoordOutput.sam results/TestSamCoordOutput.sam
else
./samTest 2> results/samTest.log && \
diff $expected/testEqWithBases.sam results/outSamEqBases.sam && \
diff $expected/testEqWithEq.sam results/outSamEqEquals.sam && \
diff $expected/testEqWithOrig.sam results/outSamEqOrig.sam && \
diff $expected/testEqWithBases.bam results/outSamEqBases.bam && \
diff $expected/testEqWithEq.bam results/outSamEqEquals.bam && \
diff $expected/testEqWithOrig.bam results/outSamEqOrig.bam && \
diff $expected/testEqWithBases.sam results/outBamEqBases.sam && \
diff $expected/testEqWithEq.sam results/outBamEqEquals.sam && \
diff $expected/testEqWithOrig.sam results/outBamEqOrig.sam && \
diff $expected/testEqWithBases.bam results/outBamEqBases.bam && \
diff $expected/testEqWithEq.bam results/outBamEqEquals.bam && \
diff $expected/testEqWithOrig.bam results/outBamEqOrig.bam && \
diff $expected/updateTagFromBam.sam results/updateTagFromBam.sam && \
diff $expected/updateTagFromSam.sam results/updateTagFromSam.sam && \
diff $expected/updateTag.bam results/updateTagFromBam.bam && \
diff $expected/updateTag.bam results/updateTagFromSam.bam && \
diff $expected/samTest.log results/samTest.log && \
diff $expected/addedTagToSam.bam results/addedTagToSam.bam && \
diff $expected/addedTagToSam.sam results/addedTagToSam.sam && \
diff $expected/addedTagToBam.sam results/addedTagToBam.sam && \
diff $expected/addedTagToBam.bam results/addedTagToBam.bam && \
diff $expected/testShift.sam results/testShift.sam && \
diff $expected/testShift.bam results/testShift.bam && \
diff $expected/testShift.sam results/testShiftFromBam.sam && \
diff $expected/testShift.bam results/testShiftFromSam.bam && \
diff $expected/TestSamCoordOutput.sam results/TestSamCoordOutput.sam
fi