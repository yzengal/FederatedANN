#!/bin/bash

ORIGINAL_DIR=$(pwd)
cd ../src/build
silo_id=2
silo_sub_id=4
port=$((50050 + ${silo_id} + ${silo_sub_id} * 5))
data_path="/home/dataset/wiki50silo/wiki50silo_${silo_id}_${silo_sub_id}.fivecs"
index_type="HNSW"
#data_path="/home/yzengal/dataset/sift/sift_base_${silo_id}.fivecs"

./silo --id=$silo_id --ip=172.21.0.16 --port=$port --data-path=$data_path --index-type=$index_type

if [ $? -ne 0 ]; then
    echo "Silo ${silo_id} FAIL"
    exit 1
fi

cd "$ORIGINAL_DIR"