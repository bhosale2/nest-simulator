#!/bin/bash

SETTINGS=" -study NEST"
OPTIONS=${SETTINGS}
EXECNAME="run"

make clean
make -j8
rm -fr ../${EXECNAME}
mkdir ../${EXECNAME}
cp ../makefiles/nest ../${EXECNAME}/executable  
cd ../${EXECNAME}
cp -r ../source/ ./source 
./executable ${OPTIONS}
