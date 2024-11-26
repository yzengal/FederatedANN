#!/bin/bash

ORIGINAL_DIR=$(pwd)
cd ../src/build

if [ -f /opt/intel/sgxsdk/environment ]; then
    source /opt/intel/sgxsdk/environment
else
    echo "Error: /opt/intel/sgxsdk/environment does not exist."
    exit 1
fi

user_ipaddr="localhost:50059"
query_path="/home/dataset/siftsmall/siftsmall_query.fvecs"
ip_path="../../scripts/ip_local.txt"
option="greedyO"
query_k=100

./sgx_broker --query-path=$query_path --user-ip=$user_ipaddr --silo-ip=$ip_path --option=$option --query-k=$query_k

if [ $? -ne 0 ]; then  
    echo "[SGX-GREEDY] Data brok FAIL"  
    exit 1  
fi 

cd "$ORIGINAL_DIR"
