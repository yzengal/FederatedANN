# Federated ANN: ANN Search Over Federated Vector Databases
This repository develops an accurate, efficient and secure solution to approximate nearest neighbor search over federated vector databases. 

## Environment

OS: Ubuntu 18.04 LTS   
GCC/G++: >= 8.4.0   
CMake: >= 3.19.1   
[gRPC](https://grpc.io): >= 1.66.0 
[Faiss library](https://github.com/facebookresearch/faiss): >= 1.9.0  
[Boost C++ library](https://www.boost.org/): >= 1.85.0  
[Microsoft SEAL](https://github.com/microsoft/SEAL): >= 4.1.0  
 
## How to install the third-party

### I. Install gRPC

gRPC is a modern, open source, high-performance remote procedure call (RPC) framework that can run anywhere. gRPC enables client and server applications to communicate transparently, and simplifies the building of connected systems.

Before compiling this project, you need to install the [gRPC](https://github.com/grpc/grpc) first by following the [guideline](https://grpc.io/docs/languages/cpp/quickstart/) as follows.

#### 1. Setup

Choose a directory to hold locally installed packages. This page assumes that the environment variable ``MY_INSTALL_DIR`` holds this directory path. For example:
```
export MY_INSTALL_DIR=$HOME/.local
mkdir -p $MY_INSTALL_DIR
```

#### 2. Install cmake
You need version 3.13 or later of cmake. Install it by following these instructions:
```
sudo apt install -y cmake
```

#### 3. Install other required tools
Install the basic tools required to build gRPC:
```
sudo apt install -y build-essential autoconf libtool pkg-config
```

#### 4. Clone the **grpc** repo
Clone the **grpc** repo and its submodules:
```
git clone --recurse-submodules -b v1.66.0 --depth 1 --shallow-submodules https://github.com/grpc/grpc
```

#### 5. Build and install **gRPC** and **Protocol Buffers**
The following commands build and locally install **gRPC** and **Protocol Buffers**:
```
cd grpc
mkdir -p cmake/build
pushd cmake/build
cmake -DgRPC_INSTALL=ON \
      -DgRPC_BUILD_TESTS=OFF \
      -DCMAKE_INSTALL_PREFIX=$MY_INSTALL_DIR \
      ../..
make -j 4
make install
popd
```

### II. Install Boost C++ Libraries

Boost provides peer-reviewed and widely useful C++ libraries that work well with the Standard Library. 
Before compiling this project, you also need to install the [Boost](https://www.boost.org/) first by following the [guideline](https://www.boost.org/doc/libs/1_85_0/more/getting_started/unix-variants.html) as follows.

#### 1. Download the **Boost**

In Ubuntu, you can use the following commands to download Boost 1.85.0
```
wget https://archives.boost.io/release/1.85.0/source/boost_1_85_0.tar.gz
tar -xzvf /boost_1_85_0.tar.gz
```

#### 2. Setup

Choose a directory to hold locally installed packages. This page assumes that the environment variable ``MY_INSTALL_DIR`` holds this directory path. For example:
```
export MY_INSTALL_DIR=$HOME/.local
mkdir -p $MY_INSTALL_DIR
```

#### 3. Build and install **Boost**
The following commands build and locally install **Boost**:
```
cd boost_1_85_0
./bootstrap.sh --prefix=$MY_INSTALL_DIR
./b2 install
```

### III. Install Faiss libraries

Facebook AI Similarity Search ([Faiss](https://github.com/facebookresearch/faiss)) is a library for efficient similarity search and clustering of dense vectors. It contains algorithms that search in sets of vectors of any size, up to ones that possibly do not fit in RAM. It also contains supporting code for evaluation and parameter tuning. Faiss is written in C++ with complete wrappers for Python/numpy. Some of the most useful algorithms are implemented on the GPU. It is developed primarily at Meta's [Fundamental AI Research](https://ai.facebook.com/) group.

#### 1. Install other required tools
Install the basic tools required to build faiss:
```
sudo apt-get install libblas-dev liblapack-dev
```

#### 2. Clone the **faiss** repo
Clone the **faiss** repo and its submodules:
```
git clone https://github.com/facebookresearch/faiss
```

#### 2. Setup

Choose a directory to hold locally installed packages. This page assumes that the environment variable ``MY_INSTALL_DIR`` holds this directory path. For example:
```
export MY_INSTALL_DIR=$HOME/.local
mkdir -p $MY_INSTALL_DIR
```

#### 3. Build and install **faiss**
The following commands build and locally install **faiss**:
```
cd faiss
cmake -B build . -DCMAKE_INSTALL_PREFIX=$MY_INSTALL_DIR -DFAISS_ENABLE_GPU=OFF -DFAISS_ENABLE_PYTHON=OFF -DFAISS_ENABLE_RAFT=OFF -DBUILD_TESTING=OFF
make -C build -j4 faiss
make -C build install
```
Note: if your **cmake** version is lower, you can directly change the cmake version in ``CMakeLists.txt``.

#### 4. Linking with **faiss** through **CMake**
It is very easy to link your own applications and libraries with **faiss** if you use **CMake**. Simply add the following to your ``CMakeLists.txt``:
```
find_package(faiss REQUIRED)
target_link_libraries(<your target> faiss)
```
If **faiss** was installed globally, the above ``find_package`` command will likely find the library automatically. 

### IV. Install SEAL (Optional)

Microsoft SEAL is an easy-to-use open-source (MIT licensed) homomorphic encryption library developed by the Cryptography and Privacy Research Group at Microsoft. Microsoft SEAL is written in modern standard C++ and is easy to compile and run in many different environments.

#### 1. Setup
Choose a directory to hold locally installed packages. This page assumes that the environment variable ``MY_INSTALL_DIR`` holds this directory path. For example:
```
export MY_INSTALL_DIR=$HOME/.local
mkdir -p $MY_INSTALL_DIR
```

#### 2. Clone the **SEAL** repo
Clone the **SEAL** repo and switch to version v4.1.1:
```
git clone https://github.com/microsoft/SEAL.git
cd SEAL
git checkout v4.1.1
```

#### 3. Build and install **SEAL**
The following commands build and locally install **SEAL**:
```
cd SEAL
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$MY_INSTALL_DIR
cmake --build build
cmake --install build
```

#### 4. Linking with **SEAL** through **CMake**
It is very easy to link your own applications and libraries with **SEAL** if you use **CMake**. Simply add the following to your ``CMakeLists.txt``:
```
find_package(SEAL 4.1 REQUIRED)
target_link_libraries(<your target> SEAL::seal)
```
If **SEAL** was installed globally, the above ``find_package`` command will likely find the library automatically. 

## Compile and run our algorithms

### Compile the algorithms

Execute the following commands to compile **data silo** and **query user**:
```
cd scripts
chmod +x *.sh
./compile.sh
```

### Enable the data silo of the algorithms

Execute the following command in one terminal to enable a data silo:
```
./Silo.sh
```
The detailed commands in the ``Silo.sh`` is as follows:
```
#!/bin/bash

ORIGINAL_DIR=$(pwd)
cd ../src/build

silo_id=1
port=$((50050 + ${silo_id}))
data_path="/home/yzengal/dataset/siftsmall/siftsmall_base_${silo_id}.fivecs"

./silo --id=$silo_id --ip=localhost --port=$port --data-path=$data_path

if [ $? -ne 0 ]; then  
    echo "Silo ${silo_id} FAIL"  
    exit 1  
fi 

cd "$ORIGINAL_DIR"
```
Here, we only need to revise the ``silo_id``, ``port`` of the IP address, and ``data_path``.

### Enable the query user of the algorithms

Similar to the script for data silos, the query user can be enabled with the following command:
```
./Server.sh
```
The detailed commands in the ``Server.sh`` is as follows:
```
#!/bin/bash

ORIGINAL_DIR=$(pwd)
cd ../src/build

query_path="/home/yzengal/dataset/siftsmall/siftsmall_query.fvecs"
ip_path="../../scripts/ip.txt"
output_path="../../scripts/siftsmall_answer.txt"
truth_path="/home/yzengal/dataset/siftsmall/siftsmall_groundtruth.ivecs"
option="greedy"

./server --query-path=$query_path --silo-ip=$ip_path --output-path=$output_path --truth-path=$truth_path --option=$option

if [ $? -ne 0 ]; then  
    echo "server FAIL"  
    exit 1  
fi 

cd "$ORIGINAL_DIR"
```
Here, we only need to revise the ``query_path`` (for query file name), ``ip_path`` (for all silo's IP addresses), ``output_path`` (for log files), ``truth_path`` (for ground truth file name), and ``option`` (for algorithm setting).
In our solution, ``option = equality`` indicates the **equality greedy algorithm (EGA)**,
``option = public`` indicates the **plaintext** baseline,
and ``option = greedX`` indicates the **optimized greedy algorithm (OGA)**.

Similarly, for the **Patching OGA algorithm (PGA)**, you can run the scripts and executable programs that are related to ``Opt_Server`` and ``Opt_Silo``.
