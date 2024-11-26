#!/bin/bash

ORIGINAL_DIR=$(pwd)
cd ../src/build

query_path="/home/dataset/siftsmall/siftsmall_query.fvecs"
ip_path="../../scripts/ip_local.txt"
output_path="../../scripts/siftsmall_answer_linear.txt"
truth_path="/home/dataset/siftsmall/siftsmall_groundtruth.ivecs"
option="greedyO"
query_k=100

./opt_server --query-path=$query_path --silo-ip=$ip_path --output-path=$output_path --truth-path=$truth_path --option=$option --query-k=$query_k

if [ $? -ne 0 ]; then  
    echo "[OPT-GREEDY] Query user FAIL"  
    exit 1  
fi 

cd "$ORIGINAL_DIR"
