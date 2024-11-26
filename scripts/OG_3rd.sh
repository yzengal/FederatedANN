#!/bin/bash

ORIGINAL_DIR=$(pwd)
cd ../src/build
silo_num=2
port=50049
ip="localhost"

./3rd_party --ip=$ip --port=$port --silo-num=$silo_num

if [ $? -ne 0 ]; then  
    echo "[OPT-GREEDY] Untrust third party ${silo_id} FAIL"  
    exit 1  
fi 

cd "$ORIGINAL_DIR"