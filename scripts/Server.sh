#!/bin/bash

ORIGINAL_DIR=$(pwd)
cd ../src/build

query_path="/home/dataset/wikipedia-1M/query.fvecs"
ip_path="../../scripts/ip.txt"
output_path="../../scripts/equality_wiki_answer_hnsw_128_5.txt"
truth_path="/home/dataset/wikipedia-1M/wikipedia_base.bin_gt_128_5.ivecs"
option="equality"
query_k=128

./server --query-path=$query_path --silo-ip=$ip_path --output-path=$output_path --truth-path=$truth_path --option=$option --query-k=$query_k

if [ $? -ne 0 ]; then  
    echo "server FAIL"  
    exit 1  
fi 

cd "$ORIGINAL_DIR"
