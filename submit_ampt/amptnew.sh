#!/bin/bash


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
path=data0
fileout=ampt-${Index}
fileout1=data0/$fileout
#fileout2=/home/lwchen/KJSUN/AMPTNew/Amptnew-Pu-bl-2-FengLi/data/$fileout
cd $path

if [ -f $fileout ]
then
    rm -rf $fileout
fi
mkdir $fileout
cd ..
cp -r ampt/code/* $fileout1
cd $fileout1


system=`uname -s`
#nrandom=`date '+%d%H%M%S'`
nrandom="`date +%N`"
case $system in
Linux|LINUX|UNIX|Unix|Darwin|DARWIN)

echo $nrandom > nseed_runtime
        ;;
    OSF*)
        echo $nrandom > nseed_runtime
        ;;
    *)
        echo "You are running operating system other than Linux/UNIX/OSF"
        echo "Modify this file first"
        exit 1;
        ;;
esac
#make -f Makefile
cp input.ampt ana/
echo "#  AMPT started at " `date` > start.time
./ampt < nseed_runtime > nohup.out
uname -n >> nohup.out
cat start.time >> nohup.out
#rm -f start.time
#cp -r ana/* $fileout2
#rm -r $fileout1
echo "#  AMPT Program finished at " `date` >> nohup.out
tput bel
#rm $fileout1/*
sleep 1
