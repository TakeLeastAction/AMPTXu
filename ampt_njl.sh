#!/bin/sh
echo "welcome jobs.process.csh"
Index=$1
nrandom=`date +%N` 
cp -r ampt/code  ./ampt_test/code$Index
cd ./ampt_test/code$Index
sh compile.sh
echo $nrandom > nseed_runtime
./ampt < nseed_runtime
sleep 20
