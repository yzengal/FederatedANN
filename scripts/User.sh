#!/bin/bash

ORIGINAL_DIR=$(pwd)
cd ../src/build

query_path="/home/dataset/audio/query.fvecs"
ip_path="../../scripts/ip_local.txt"
output_path="../../scripts/audio_answer_hnsw.txt"
truth_path="/home/dataset/audio/audio_gt_128_5.ivecs"
option="public"

./server --query-path=$query_path --silo-ip=$ip_path --output-path=$output_path --truth-path=$truth_path --option=$option

if [ $? -ne 0 ]; then  
    echo "server FAIL"  
    exit 1  
fi 

cd "$ORIGINAL_DIR"
