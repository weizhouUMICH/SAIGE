#!/bin/bash
dir=$(dirname ${0})
basedir=${dir}/../..
echo ${dir}
pybot -o  output.html -l log.html -r NONE -d ${basedir}/build/test/functional/test-reports -x xunit.xml ${dir}/tests.txt
