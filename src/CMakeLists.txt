cmake_minimum_required(VERSION 3.12)

project(FederatedVectorQuery VERSION 1.0 DESCRIPTION "A bechmark project for query processing over federated vector databases" LANGUAGES CXX)  
set(CMAKE_VERBOSE_MAKEFILE OFF)

set(LOCAL_DEBUG OFF CACHE BOOL "Enable debugging log information") 
If (LOCAL_DEBUG)
    message(STATUS "Enable debugging log information")  
else()
    message(STATUS "Disable debugging log information")  
endif()

if(USE_SGX STREQUAL "ON")
    message(STATUS "Compiling SGX MPC:")

    set(Enclave_Test_Key "${CMAKE_BINARY_DIR}/Enclave_private_test.pem")  
    if(NOT EXISTS ${Enclave_Test_Key})  
        message(STATUS "There is no enclave test key<${Enclave_Test_Key}>.")  
        message(STATUS "The project will generate a key<${Enclave_Test_Key}> for test.")  
        execute_process(  
            COMMAND openssl genrsa -out ${Enclave_Test_Key} -3 3072  
            RESULT_VARIABLE openssl_result  
            OUTPUT_QUIET  
            ERROR_QUIET  
        )  
        if(NOT openssl_result EQUAL 0)  
            message(FATAL_ERROR "Failed to generate RSA key.")  
        endif()  
    endif()

    list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/sgx/cmake")
    find_package(SGX REQUIRED)

    message("CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

    set(EDL_SEARCH_PATHS ${CMAKE_SOURCE_DIR}/sgx/Enclave)
    set(E_SRCS ${CMAKE_SOURCE_DIR}/sgx/Enclave/Enclave.cpp)

    if(SGX_HW AND SGX_MODE STREQUAL "Release")
        set(LDS ${CMAKE_SOURCE_DIR}/sgx/Enclave/Enclave.lds)
    else()
        set(LDS ${CMAKE_SOURCE_DIR}/sgx/Enclave/Enclave_debug.lds)
    endif()

    add_enclave_library(enclave 
        SRCS ${E_SRCS} 
        EDL ${EDL_SEARCH_PATHS}/Enclave.edl 
        EDL_SEARCH_PATHS ${EDL_SEARCH_PATHS} 
        LDSCRIPT ${LDS})

    message(STATUS "Enclave_Test_Key=${Enclave_Test_Key}")
    enclave_sign(enclave 
        KEY ${Enclave_Test_Key} 
        CONFIG ${EDL_SEARCH_PATHS}/Enclave.config.xml)

    set(SRCS ${CMAKE_SOURCE_DIR}/sgx/mpc/SgxMpc.cpp)
    add_untrusted_library(SgxMpc 
        STATIC 
        SRCS ${SRCS} 
        EDL ${EDL_SEARCH_PATHS}/Enclave.edl
        EDL_SEARCH_PATHS ${EDL_SEARCH_PATHS})
endif()

set(CMAKE_CXX_STANDARD 17)  
set(CMAKE_CXX_STANDARD_REQUIRED True)  
   
find_package(Boost 1.85.0 COMPONENTS program_options system REQUIRED)
If (NOT Boost_FOUND)
    message(FATAL_ERROR "Could not find Boost version 1.85.0 or higher with the require components")
endif()

find_package(OpenMP REQUIRED)
If (NOT OpenMP_FOUND)
    message(FATAL_ERROR "Could not find OpenMP")  
endif()

find_package(BLAS REQUIRED)
If (NOT BLAS_FOUND)
    message(FATAL_ERROR "Could not find BLAS")  
endif()

find_package(LAPACK REQUIRED)
If (NOT LAPACK_FOUND)
    message(FATAL_ERROR "Could not find LAPACK")  
endif()

find_package(faiss REQUIRED)
If (NOT faiss_FOUND)
    message(FATAL_ERROR "Could not find faiss")  
endif() 

include(./common.cmake)
include_directories("${CMAKE_CURRENT_BINARY_DIR}")

set(project_protos
    "./proto/FedVectorDB.proto"
    "./proto/FedVectorDBO.proto"
    "./proto/UntrustThirdParty.proto"
    "./proto/SgxSilo.proto"
    "./proto/SgxUser.proto"
) 

foreach(proto_file IN LISTS project_protos)
    message(STATUS "Processing proto file: ${proto_file}")

    get_filename_component(proto_name ${proto_file} NAME_WE) 
    get_filename_component(database_proto ${proto_file} ABSOLUTE)
    get_filename_component(database_proto_path "${database_proto}" PATH)

    #compile proto file
    set(database_proto_srcs "${CMAKE_CURRENT_BINARY_DIR}/${proto_name}.pb.cc")
    set(database_proto_hdrs "${CMAKE_CURRENT_BINARY_DIR}/${proto_name}.pb.h")
    set(database_grpc_srcs "${CMAKE_CURRENT_BINARY_DIR}/${proto_name}.grpc.pb.cc")
    set(database_grpc_hdrs "${CMAKE_CURRENT_BINARY_DIR}/${proto_name}.grpc.pb.h")

    add_custom_command(
        OUTPUT "${database_proto_srcs}" "${database_proto_hdrs}" "${database_grpc_srcs}" "${database_grpc_hdrs}"
        COMMAND ${_PROTOBUF_PROTOC}
        ARGS --grpc_out "${CMAKE_CURRENT_BINARY_DIR}"
            --cpp_out "${CMAKE_CURRENT_BINARY_DIR}"
            -I "${database_proto_path}"
            --plugin=protoc-gen-grpc="${_GRPC_CPP_PLUGIN_EXECUTABLE}"
            "${database_proto}"
        DEPENDS "${database_proto}"
        COMMENT "Generating ${proto_name}.pb.h, ${proto_name}.pb.cc, ${proto_name}.grpc.pb.h, and ${proto_name}.grpc.pb.cc"
    )

    add_library("${proto_name}_proto"
        ${database_grpc_srcs}
        ${database_grpc_hdrs}
        ${database_proto_srcs}
        ${database_proto_hdrs})

    target_link_libraries("${proto_name}_proto"
        ${_REFLECTION}
        ${_GRPC_GRPCPP}
        ${_PROTOBUF_LIBPROTOBUF})
endforeach()

add_library(File_IO
    "./utils/File_IO.h"
    "./utils/File_IO.cpp")

set(_target "silo")
message(STATUS "Compiling executable ${_target}")
add_executable(${_target} "./Silo.cpp")
if(LOCAL_DEBUG)
    target_compile_definitions(${_target} PRIVATE LOCAL_DEBUG)
endif()
target_link_libraries(${_target} PRIVATE
    FedVectorDB_proto
    File_IO
    "pthread"
    ${_REFLECTION}
    ${_GRPC_GRPCPP}
    ${_PROTOBUF_LIBPROTOBUF})
target_link_libraries(${_target} PRIVATE faiss OpenMP::OpenMP_CXX blas lapack Boost::program_options Boost::system)

set(_target "server")
message(STATUS "Compiling executable ${_target}")
add_executable(${_target} "./Server.cpp")
if(LOCAL_DEBUG)
    target_compile_definitions(${_target} PRIVATE LOCAL_DEBUG)
endif()
target_link_libraries(${_target} PRIVATE
    FedVectorDB_proto
    File_IO
    "pthread"
    ${_REFLECTION}
    ${_GRPC_GRPCPP}
    ${_PROTOBUF_LIBPROTOBUF})
target_link_libraries(${_target} PRIVATE Boost::program_options Boost::system)


set(_target "opt_silo")
message(STATUS "Compiling executable ${_target}")
add_executable(${_target} "./Opt_Silo.cpp")
if(LOCAL_DEBUG)
    target_compile_definitions(${_target} PRIVATE LOCAL_DEBUG)
endif()
target_link_libraries(${_target} PRIVATE
    FedVectorDBO_proto
    UntrustThirdParty_proto
    File_IO
    "pthread"
    ${_REFLECTION}
    ${_GRPC_GRPCPP}
    ${_PROTOBUF_LIBPROTOBUF})
target_link_libraries(${_target} PRIVATE faiss OpenMP::OpenMP_CXX blas lapack Boost::program_options Boost::system)

set(_target "opt_server")
message(STATUS "Compiling executable ${_target}")
add_executable(${_target} "./Opt_Server.cpp")
if(LOCAL_DEBUG)
    target_compile_definitions(${_target} PRIVATE LOCAL_DEBUG)
endif()
target_link_libraries(${_target} PRIVATE
    FedVectorDBO_proto
    UntrustThirdParty_proto
    File_IO
    "pthread"
    ${_REFLECTION}
    ${_GRPC_GRPCPP}
    ${_PROTOBUF_LIBPROTOBUF})
target_link_libraries(${_target} PRIVATE Boost::program_options Boost::system)

set(_target "3rd_party")
message(STATUS "Compiling executable ${_target}")
add_executable(${_target} "./Opt_Untrust3Party.cpp")
if(LOCAL_DEBUG)
    target_compile_definitions(${_target} PRIVATE LOCAL_DEBUG)
endif()
target_link_libraries(${_target} PRIVATE
    UntrustThirdParty_proto
    File_IO
    "pthread"
    ${_REFLECTION}
    ${_GRPC_GRPCPP}
    ${_PROTOBUF_LIBPROTOBUF})
target_link_libraries(${_target} PRIVATE Boost::program_options Boost::system)

add_library(AES
    "./crypto/AES.h"
    "./crypto/AES.cpp")

set(_target "sgx_silo")
message(STATUS "Compiling executable ${_target}")
add_executable(${_target} "./Sgx_Silo.cpp")
if(LOCAL_DEBUG)
    target_compile_definitions(${_target} PRIVATE LOCAL_DEBUG)
endif()
target_link_libraries(${_target} PRIVATE
    SgxSilo_proto
    File_IO
    AES
    "pthread"
    ${_REFLECTION}
    ${_GRPC_GRPCPP}
    ${_PROTOBUF_LIBPROTOBUF})
target_link_libraries(${_target} PRIVATE faiss OpenMP::OpenMP_CXX blas lapack Boost::program_options Boost::system)

set(_target "sgx_user")
message(STATUS "Compiling executable ${_target}")
add_executable(${_target} "./Sgx_User.cpp")
if(LOCAL_DEBUG)
    target_compile_definitions(${_target} PRIVATE LOCAL_DEBUG)
endif()
target_link_libraries(${_target} PRIVATE
    SgxSilo_proto
    SgxUser_proto
    File_IO
    "pthread"
    ${_REFLECTION}
    ${_GRPC_GRPCPP}
    ${_PROTOBUF_LIBPROTOBUF})
target_link_libraries(${_target} PRIVATE Boost::program_options Boost::system)

set(_target "sgx_broker")
message(STATUS "Compiling executable ${_target}")
if(USE_SGX STREQUAL "ON")
    set(SRCS "${CMAKE_SOURCE_DIR}/Sgx_Broker.cpp")
    message(STATUS "SRCS=${SRCS}")
    message(STATUS "EDL=${CMAKE_SOURCE_DIR}/sgx/Enclave/Enclave.edl")
    add_untrusted_executable(${_target} 
        SRCS ${SRCS}
        EDL ${EDL_SEARCH_PATHS}/Enclave.edl
        EDL_SEARCH_PATHS ${EDL_SEARCH_PATHS})
    target_compile_definitions(${_target} PRIVATE USE_INTEL_SGX)
    if(LOCAL_DEBUG STREQUAL "ON")
        target_compile_definitions(${_target} PRIVATE LOCAL_DEBUG)
    endif()
    add_dependencies(${_target} enclave-sign)
    target_link_libraries(${_target}
        SgxSilo_proto
        SgxUser_proto
        File_IO
        AES
        SgxMpc
        "pthread"
        ${_REFLECTION}
        ${_GRPC_GRPCPP}
        ${_PROTOBUF_LIBPROTOBUF})
    target_link_libraries(${_target} Boost::program_options Boost::system)
else()
    add_executable(${_target} "./Sgx_Broker.cpp")
    target_compile_definitions(${_target} PRIVATE USE_TRUST_BROKER)
    if(LOCAL_DEBUG)
        target_compile_definitions(${_target} PRIVATE LOCAL_DEBUG)
    endif()
    target_link_libraries(${_target} PRIVATE
        SgxSilo_proto
        SgxUser_proto
        File_IO
        AES
        "pthread"
        ${_REFLECTION}
        ${_GRPC_GRPCPP}
        ${_PROTOBUF_LIBPROTOBUF})
    target_link_libraries(${_target} PRIVATE Boost::program_options Boost::system)
endif()