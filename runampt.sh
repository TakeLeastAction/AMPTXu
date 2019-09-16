#!/bin/sh
rm -rf ampt_test
mkdir ampt_test
cp -r ampt/input ampt_test

condor_submit job* 


