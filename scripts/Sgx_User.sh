#!/bin/bash

ORIGINAL_DIR=$(pwd)
cd ../src/build
port="50059"
ip_path="../../scripts/ip_local.txt"
output_path="../../scripts/siftsmall_answer_linear.txt"
truth_path="/home/dataset/siftsmall/siftsmall_groundtruth.ivecs"

./sgx_user --ip=localhost --port=$port --silo-ip=$ip_path --output-path=$output_path --truth-path=$truth_path

if [ $? -ne 0 ]; then  
    echo "[SGX-GREEDY] Query user FAIL"  
    exit 1  
fi 

cd "$ORIGINAL_DIR"
