#!/bin/bash

source /home/comp/root53424/bin/thisroot.sh
export LD_LIBRARY_PATH=/usr/lib:/usr/lib64:/lib:/lib64:/usr/local/lib:/usr/local/lib64:$LD_LIBRARY_PATH

Index=$(($1+1))

ad=00
ad1=0

if [ ${Index} -lt 10 ]
then
    echo ${Index}
    Index=${ad}${Index} 
    echo ${Index}
elif [ ${Index} -lt 100 ]
then
    Index=${ad1}${Index}
fi



#Index=250
#Index1=$((Index/5))
#Index2=$((Index%5))
#cho $Index1 $Index2
#path=/home/zhaoxinli/AMPT/test/amptree1
cd ..
s=pwd
path=data0
fileout=ampt-${Index}
fileout1=data0/$fileout
#fileout2=/home/lwchen/KJSUN/AMPTNew/Amptnew-Pu-bl-2-FengLi/data/$fileout
#cd $path

#if [ -f $fileout ]
#then
#    rm -rf $fileout
#fi
#mkdir $fileout
#cd ..
#ifort -o coal CoalNucl.f90
cp -r coal $fileout1
cd $fileout1

rm yield*.dat
./coal
tput bel
#rm $fileout1/*
sleep 1
