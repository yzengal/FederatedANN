#!/bin/bash

ORIGINAL_DIR=$(pwd)
cd ../src/build
silo_id=0
port=$((50050 + ${silo_id}))
data_path="/home/dataset/audio/audio_base_5.fivecs"
#data_path="/home/yzengal/dataset/sift/sift_base_${silo_id}.fivecs"

./silo --id=$silo_id --ip=localhost --port=$port --data-path=$data_path

if [ $? -ne 0 ]; then  
    echo "Silo ${silo_id} FAIL"  
    exit 1  
fi 

cd "$ORIGINAL_DIR"