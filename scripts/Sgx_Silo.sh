#!/bin/bash

silo_id=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --id=*)
            silo_id="${1#*=}"
            ;;
        --id)
            shift
            if [[ $# -gt 0 && "$1" != --* ]]; then
                silo_id="$1"
            else
                echo "Error: --id requires an argument"
                exit 1
            fi
            ;;
        *)
            silo_id="$1"
            ;;
    esac
    shift
done

ORIGINAL_DIR=$(pwd)
cd ../src/build
silo_id=$((silo_id))
port=$((50050 + ${silo_id}))
data_path="/home/dataset/siftsmall/2/siftsmall_base_${silo_id}.fivecs"
index_type="Linear"

./sgx_silo --id=$silo_id --ip=localhost --port=$port --data-path=$data_path --index-type=$index_type

if [ $? -ne 0 ]; then  
    echo "[SGX-GREEDY] Data Silo ${silo_id} FAIL"  
    exit 1  
fi 

cd "$ORIGINAL_DIR"