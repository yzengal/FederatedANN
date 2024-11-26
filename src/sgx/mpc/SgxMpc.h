#ifndef _SGX_MPC_H_
#define _SGX_MPC_H_

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <vector>

#include "sgx_error.h"       /* sgx_status_t */
#include "sgx_eid.h"     /* sgx_enclave_id_t */

#ifndef TRUE
# define TRUE 1
#endif

#ifndef FALSE
# define FALSE 0
#endif

#if   defined(__GNUC__)
# define ENCLAVE_FILENAME "enclave.signed.so"
#endif

extern sgx_enclave_id_t global_eid;    /* global enclave id */

//sgx_status_t SGX_CDECL ocall_print_string(const char* str);

void SgxInitEnclave();

void SgxFreeEnclave();

std::vector<int> SgxPSI(const uint8_t *aes_key, const uint8_t* aes_iv,
                        const unsigned char* a_encrypt_key, size_t a_size, size_t a_num, 
                        const unsigned char* b_encrypt_key, size_t b_size, size_t b_num);

void SgxOblivCreateQueue(int queue_size);

void SgxOblivFreeQueue(void);

void SgxOblivInitQueue(int queue_size);

void SgxEnqueue(const int silo_id, const int silo_beg_idx, const int num,
                const uint8_t *aes_key, const uint8_t* aes_iv,
                const unsigned char* encrypt_data, size_t data_size);

int SgxGetSiloIndex(int silo_id);

int SgxGetDataIndex(int data_id);

void SgxOblivQueueHead(int* idx, int* silo_id);

void SgxOblivDequeue(void);

void SgxOblivDequeue(int num);

#endif // _SGX_MPC_H_
