#!/bin/bash
## utility script to compile and run the macro for EFTFitter

###
if [ -z ${1} ]; then
    echo 'No macro detected! Aborting...'
    echo 'Usage: ./exec/execMacro.sh EFTFITTER_EXEC_MACRO'
    echo 'To log and watch simultaneously: exec/execMacro.sh EFTFITTER_EXEC_MACRO |& tee LOGFILE'
    exit
fi

if [ ! -f ${1} ]; then
    echo 'No macro detected. Aborting...'
    echo 'Usage: ./exec/execMacro.sh EFTFITTER_EXEC_MACRO'
    echo 'To log and watch simultaneously: exec/execMacro.sh EFTFITTER_EXEC_MACRO |& tee LOGFILE'
    exit
fi

## quick check that ROOT is present
if [ -z ${ROOTSYS} ]; then
    echo "ROOT is not present, please ensure that it is! Aborting..."
    exit
fi

## and that it is the right version
rooVer=$(root-config --version)
if [ $(echo "${rooVer%%/*} < 6.00" | bc -l) -eq 1 ]; then
    echo "Detected ROOT version < 6.00, while a newer version is required. Aborting..."
    exit
fi

## ok, do stuff
eftSrc="/nfs/dust/cms/user/afiqaize/cms/rand/eftRivet_290118/EFTFitter/src/"
txtMacro=$(readlink -f ${1})
exeMacro=${txtMacro%.*}

## set the compiler options to be used
compileOpt=(`root-config --cflags --evelibs`)
compileOpt+=(-std=c++14 -O3 -Wall -Wextra -Wpedantic -Werror)

## eh these are fine
compileOpt+=(-Wno-float-equal -Wno-sign-compare)

g++ "${compileOpt[@]}" -o ${exeMacro} ${txtMacro} ${eftSrc}/EFTFitter.cc
if ((${?} == 0)); then
    ${exeMacro}
    sleep 1
    rm ${exeMacro}
fi
