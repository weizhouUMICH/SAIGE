ERROR=false

./fastqTest > results/fastqTest.txt
diff  results/fastqTest.txt expectedResults/ExpectedResultsFastqTestResults.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

if($ERROR == true)
then
  exit 1
fi

