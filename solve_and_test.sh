#!/bin/bash

DIR_SOLVER=.
DIR_CHECKER=./test
SOLVER=${DIR_SOLVER}/for-ch_solver
CHECKER=${DIR_CHECKER}/pedibus_checker.py

rm -f *.sol && \
    ${SOLVER} $@ && \
    FILE_SOL=$(ls | grep *.sol) && \
    FILE_DAT=${DIR_CHECKER}/$(echo ${FILE_SOL} | sed 's/\.sol/\.dat/g') && \
    echo "Solution file: ${FILE_SOL}" && \
    echo "Data file: ${FILE_DAT}" && \
    python2 ${CHECKER} ${FILE_DAT} ${DIR_SOLVER}/${FILE_SOL}
