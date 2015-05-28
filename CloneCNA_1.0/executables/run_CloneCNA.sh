#!/bin/bash
# script for execution of CloneCNA algorithm
# 28/03/2015 by Zhenhua, yzh163@mail.ustc.edu.cn

echo 

function usage() {
	echo -e "$0 <MCRDIR> <TumorCountFile> <NormalCountFile> <TumorDepthFile> <gcFile> <outputDir> <plotDir>\n"
}

if [ $# -ne 7 ]; then
	echo -e "7 parameters are required!\n"
	echo "Usage:" 
	usage
	exit 0
fi

MCRDIR=$1
countFile=$2
refFile=$3
depthFile=$4
gcFile=$5
outputDir=$6
plotDir=$7

MCRDIR=${MCRDIR%/}
countFile=${countFile%/}
refFile=${refFile%/}
depthFile=${depthFile%/}
gcFile=${gcFile%/}
outputDir=${outputDir%/}
plotDir=${plotDir%/}

[ ! -e $MCRDIR ] && echo "Error: MCR directory $MCRDIR does not exist!" && exit 0
[ ! -e $countFile ] && echo "Error: Tumor read counts file $countFile does not exist!" && exit 0
[ ! -e $refFile ] && echo "Error: Normal read counts file $refFile does not exist!" && exit 0
[ ! -e $depthFile ] && echo "Error: Tumor allelic depths file $depthFile does not exist!" && exit 0
[ ! -e $gcFile ] && echo "Error: GC file $gcFile does not exist!" && exit 0

[ ! -e $outputDir ] && echo "Create the directory $outputDir!" && mkdir -p $outputDir
[ ! -e $plotDir ] && echo "Create the directory $plotDir!" && mkdir -p $plotDir

echo "Setting up environment variables"
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRDIR}/runtime/glnxa64:${MCRDIR}/bin/glnxa64:${MCRDIR}/sys/os/glnxa64:${MCRDIR}/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:${MCRDIR}/sys/java/jre/glnxa64/jre/lib/amd64/server:${MCRDIR}/sys/java/jre/glnxa64/jre/lib/amd64
#LD_LIBRARY_PATH=.:${MCRDIR}/runtime/glnxa64
#LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRDIR}/bin/glnxa64
#LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRDIR}/sys/os/glnxa64
XAPPLRESDIR=${MCRDIR}/X11/app-defaults
export XAPPLRESDIR
export LD_LIBRARY_PATH

eval ./CloneCNA $countFile $refFile $depthFile $gcFile $outputDir $plotDir

exit 0
