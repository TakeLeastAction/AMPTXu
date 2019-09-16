#!/bin/bash
rm -rf data0
mkdir data0
cp -r ampt/input data0

cd submit_ampt
condor_submit ampt.con 


