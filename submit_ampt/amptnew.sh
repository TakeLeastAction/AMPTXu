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
cp -r ampt/* $fileout1
cd $fileout1

EFRM=2.5
FRAME=CMS
PROJ=A
TARG=A
IAP=197
IZP=79
IAT=197
IZT=79
NEVENT=100
BMIN=0
BMAX=3
ISOFT=4
NTMAX=150
DT=0.2
PARJ1=0.55
PARJ2=0.15
POPCORN=1
PARJ5=1.0
SHF=1
QUF=0
QUPA=2.0
P0=2.0
PASCMA=2.265d0
IZPC=0
ALPHA=0.33d0
DPCOAL=1d6
DRCOAL=1d6
IHJSED=11
RANSEED=13150909
RANSEEP=8
FLAGK=0
FLAGPHI=1
FLAGPI0=0
OSCAR=0
FLAGDEU=0
INTEDEU=1
CROSEC=1
PT=-7
MAXMISS=1000
FLAGRADIA=3
FLAGKT=1
FLAGQUA=0
PX=7.
PY=0.
X=0.
Y=0.
NSEMBD=1
PSEMBD=5.
TMAXEMBD=0.
FLAGMOD=0
FACTOR=1.d0
FLAGORIEN=0
#define parameters

#infile=input.ampt
#if [ -f ${infile} ]
#then
#    rm -f ${infile}
#fi

#echo -e $EFRM'\n'$FRAME'\n'$PROJ'\n'$TARG'\n'$IAP              >>${infile}
#echo -e $IZP'\n'$IAT'\n'$IZT'\n'$NEVENT'\n'$BMIN               >>${infile}
#echo -e $BMAX'\n'$ISOFT'\n'$NTMAX'\n'$DT'\n'$PARJ1             >>${infile}
#echo -e $PARJ2'\n'$POPCORN'\n'$PARJ5'\n'$SHF'\n'$QUF           >>${infile}
#echo -e $QUPA'\n'$P0'\n'$PASCMA'\n'$IZPC'\n'$ALPHA             >>${infile}
#echo -e $DPCOAL'\n'$DRCOAL'\n'$IHJSED'\n'$RANSEED'\n'$RANSEEP  >>${infile}
#echo -e $FLAGK'\n'$FLAGPHI'\n'$FLAGPI0'\n'$OSCAR'\n'$FLAGDEU   >>${infile}
#echo -e $INTEDEU'\n'$CROSEC'\n'$PT'\n'$MAXMISS'\n'$FLAGRADIA   >>${infile}
#echo -e $FLAGKT'\n'$FLAGQUA'\n'$PX,$PY'\n'$X,$Y                >>${infile}
#echo -e $NSEMBD,$PSEMBD,$TMAXEMBD'\n'$FLAGMOD'\n'$FACTOR'\n'$FLAGORIEN'\n'   >>${infile}
#parameters input infile


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
