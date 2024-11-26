/*
 * Copyright (C) 2011-2021 Intel Corporation. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in
 *     the documentation and/or other materials provided with the
 *     distribution.
 *   * Neither the name of Intel Corporation nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
#include <stdarg.h>
#include <stdio.h>      /* vsnprintf */
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "Enclave.h"
#include "Enclave_t.h"
#include "Enclave_t.c"
#include "mbedtls/aes.h"
#include "mbedtls/cipher.h"
#include "mbedtls/entropy.h"
#include "mbedtls/ctr_drbg.h"

#define FAIL_AES	0x2

// #define LOCAL_SGX_DEBUG

extern sgx_status_t SGX_CDECL ocall_printf(const char* msg);

/* 
 * local_printf: 
 *   Invokes OCALL to display the enclave buffer to the terminal.
 *   'local_printf' function is required for sgx protobuf logging module.
 */
 
static int local_printf(const char *fmt, ...)
{
    char buf[BUFSIZ] = {'\0'};
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, BUFSIZ, fmt, ap);
    va_end(ap);
    ocall_print_string(buf);
    return 0;
}

static float UnsignedVectorToFloat(const unsigned char* arr, size_t offset) {
    float result;
    memcpy(&result, arr+offset, sizeof(float));
    return result;
}

static int32_t UnsignedVectorToInt32(const unsigned char* arr, size_t offset) {
    int32_t result;
    memcpy(&result, arr+offset, sizeof(int32_t));
    return result;
}

static int mbedtls_decrypt_data(const uint8_t *aes_key, const uint8_t* aes_iv,
                                const unsigned char* encrypt_key, size_t encrypt_key_size, 
                                size_t key_num, unsigned char*& decrypt_key) {
    const int key_size = 16;
    uint8_t key[key_size], iv[key_size];
    int ret_status = 0;

    memcpy(key, aes_key, sizeof(uint8_t)*key_size);
    memcpy(iv, aes_iv, sizeof(uint8_t)*key_size);

    #ifdef LOCAL_SGX_DEBUG
    local_printf("[FINISH] memcpy key and iv\n");
    local_printf("Key: ");
    for (int i = 0; i < key_size; ++i) {
        local_printf("%02x", key[i]);
    }
    local_printf("\n");
    local_printf(" Iv: ");
    for (int i = 0; i < key_size; ++i) {
        local_printf("%02x", iv[i]);
    }
    local_printf("\n");
    #endif

    mbedtls_aes_context aes_context;
    mbedtls_aes_init(&aes_context);

    #ifdef LOCAL_SGX_DEBUG
    local_printf("[FINISH] mbedtls_aes_init\n");
    #endif

    ret_status = mbedtls_aes_setkey_dec(&aes_context, key, 128);
    if(ret_status != 0){
        #ifdef LOCAL_SGX_DEBUG
        local_printf("failed to set decryption key\n");
        #endif
        mbedtls_aes_free(&aes_context);
        return -1;
    }

    decrypt_key = (unsigned char *)calloc(encrypt_key_size, sizeof(unsigned char));
    if (decrypt_key == NULL) {
        #ifdef LOCAL_SGX_DEBUG
        local_printf("failed to allocate memory for decrypt_key\n");
        #endif
        mbedtls_aes_free(&aes_context);
        return -1;
    }

    ret_status = mbedtls_aes_crypt_cbc(&aes_context, MBEDTLS_AES_DECRYPT, encrypt_key_size, iv, encrypt_key, decrypt_key);
    if (ret_status != 0){
        #ifdef LOCAL_SGX_DEBUG
        local_printf("failed to decrypt data\n");
        #endif
        free(decrypt_key);
        decrypt_key = NULL;
        mbedtls_aes_free(&aes_context);
        return -1;
    }

    mbedtls_aes_free(&aes_context);
    return 0;
}

static int mbedtls_decode(const unsigned char* decrypt_key) {
    int value = 0;  
    for (size_t i = 0; i < sizeof(int); ++i) {  
        int tmp = decrypt_key[i];
        value |= (tmp << (i * 8));
    }  
    return value;
}

int ecall_PSI(const uint8_t *aes_key, const uint8_t* aes_iv,
                const unsigned char* a_encrypt_key, size_t a_size, size_t a_num, 
                const unsigned char* b_encrypt_key, size_t b_size, size_t b_num, 
                int* psi_set, size_t result_buffer_size, size_t* psi_num) {
    
    int ret_status;

    #ifdef LOCAL_SGX_DEBUG
    local_printf("[START] ecall_PSI\n");
    #endif

    unsigned char* a_decrypt_key = NULL;
    ret_status = mbedtls_decrypt_data(aes_key, aes_iv, a_encrypt_key, a_size, a_num, a_decrypt_key);
    if (0 != ret_status) return ret_status;

    #ifdef LOCAL_SGX_DEBUG
    local_printf("[FINISH] decrypt a\n");
    #endif

    unsigned char* b_decrypt_key = NULL;
    ret_status = mbedtls_decrypt_data(aes_key, aes_iv, b_encrypt_key, b_size, b_num, b_decrypt_key);
    if (0 != ret_status) return ret_status;

    #ifdef LOCAL_SGX_DEBUG
    local_printf("[FINISH] decrypt b\n");
    #endif
    
    size_t a_idx = 0, b_idx = 0, res_idx = 0;
    if (a_num>0 && b_num>0) {
        int a_key = mbedtls_decode(a_decrypt_key+a_idx*sizeof(int));
        int b_key = mbedtls_decode(b_decrypt_key+b_idx*sizeof(int));
        while (a_idx<a_num && b_idx<b_num) {
            if (a_key < b_key) {
                if (++a_idx < a_num)
                    a_key = mbedtls_decode(a_decrypt_key+a_idx*sizeof(int));
            } else if (a_key > b_key) {
                if (++b_idx < b_num)
                    b_key = mbedtls_decode(b_decrypt_key+b_idx*sizeof(int));
            } else {
                psi_set[res_idx++] = a_key;
                if (++a_idx < a_num)
                    a_key = mbedtls_decode(a_decrypt_key+a_idx*sizeof(int));
                if (++b_idx < b_num)
                    b_key = mbedtls_decode(b_decrypt_key+b_idx*sizeof(int));
            }
        }
    }
    *psi_num = res_idx;

    free(a_decrypt_key);
    free(b_decrypt_key);

    return ret_status;
}

#define INF_PAIR_FIRST       1e20

struct QueueItem {
	float dist;
	int vid;
    int silo_id;
};

struct OblivQueue {
    int N;
    int ptr;
    struct QueueItem* a;
    struct QueueItem* b;
};

static struct OblivQueue* obliv_queue = NULL;

static void compare_and_swap(struct QueueItem* a, struct QueueItem* b) {
	float dista = a->dist;
	float distb = b->dist;

	if(dista < distb) {
	    int vida = a->vid, vidb = b->vid;
        int siloa = a->silo_id, silob = b->silo_id;

        a->dist = distb; a->vid = vidb; a->silo_id = silob;
        b->dist = dista; b->vid = vida; b->silo_id = siloa;
	}
}

static void obliv_bitonic(struct QueueItem *mem, int N) {
	int K = log2(N);
	int d = 1 << K;
	for (int n = 0; n < (d >> 1); n++) {
		compare_and_swap(&mem[n], &mem[d - n - 1]);
	}
	K--;
	if (K <= 0) {
		return;
	}
	for (int k = K; k > 0; k--) {
		d = 1 << k;
		for (int m = 0; m < N; m += d) {
			for (int n = 0; n < (d >> 1); n++) {
				compare_and_swap(&mem[m + n], &mem[m + (d >> 1) + n]);
			}
		}
	}
}

static void obliv_bitonic_sort(struct QueueItem *mem, int N) {
	struct QueueItem* map;
	map = (struct QueueItem *)malloc(N * sizeof(struct QueueItem));
	for (int i = 0; i < N; i++) {
		map[i].dist = mem[i].dist;
        map[i].vid = mem[i].vid;
        map[i].silo_id = mem[i].silo_id;
	}
	int K = log2(N);
	for (int k = 1; k <= K; k++) {
		int d = 1 << k;
		for (int n = 0; n < N; n += d) {
			struct QueueItem* map_ptr = &map[n];
			obliv_bitonic(map_ptr, d);
		}
	}
	for (int n = 0; n < N; n++) {
		mem[n].dist = map[n].dist;
        mem[n].vid = map[n].vid;
        mem[n].silo_id = map[n].silo_id;
	}
	free(map);
}

static void obliv_resort(struct OblivQueue *q) {
    obliv_bitonic_sort(q->a, q->N);
}

static void obliv_append(struct OblivQueue *q, float dist, int vid, int silo_id) {
    int pos = 0;
    for (int i=0; i<q->ptr; ++i) {
        if (q->a[i].dist > dist) {
            pos = i + 1;
        }
    }
    for (int i=0; i<q->ptr+1; ++i) {
        if (i < pos) {
            q->b[i].dist = q->a[i].dist;
            q->b[i].vid = q->a[i].vid;
            q->b[i].silo_id = q->a[i].silo_id;
        } else if (i == pos) {
            q->b[i].dist = dist;
            q->b[i].vid = vid;
            q->b[i].silo_id = silo_id;
        } else {
            q->b[i].dist = q->a[i-1].dist;
            q->b[i].vid = q->a[i-1].vid;
            q->b[i].silo_id = q->a[i-1].silo_id;
        }
    }
    q->ptr += 1;
    for (int i=0; i<q->ptr; ++i) {
        q->a[i].dist = q->b[i].dist;
        q->a[i].vid = q->b[i].vid;
        q->a[i].silo_id = q->b[i].silo_id;
    }
}

static void obliv_init_queue(struct OblivQueue* q, int queue_size) {
    if (q->N != queue_size) {
        if (q->a != NULL) free(q->a);
        q->a = (struct QueueItem*) malloc(queue_size * sizeof(struct QueueItem));
        if (q->b != NULL) free(q->b);
        q->b = (struct QueueItem*) malloc(queue_size * sizeof(struct QueueItem));      
    }
    q->N = queue_size;
    q->ptr = 0;
    for(int i=0; i<queue_size; i++) {
        q->a[i].dist = INF_PAIR_FIRST;
        q->a[i].vid = -1;
        q->a[i].silo_id = -1;
    }
}

static void OblivDequeue(struct OblivQueue *q) {
    q->ptr -= 1;
    q->a[q->ptr].dist = INF_PAIR_FIRST;
    q->a[q->ptr].vid = -1;
    q->a[q->ptr].silo_id = -1;
    for (int i=0; i<q->ptr; ++i) {
        q->a[i].dist = q->a[i].dist;
        q->a[i].vid = q->a[i].vid;
        q->a[i].silo_id = q->a[i].silo_id;
    }
}

static void OblivDequeue(struct OblivQueue *q, int num) {
    for (int i=0; i<num; ++i) {
        q->ptr -= 1;
        q->a[q->ptr].dist = INF_PAIR_FIRST;
        q->a[q->ptr].vid = -1;
        q->a[q->ptr].silo_id = -1;
    }
    for (int i=0; i<q->ptr; ++i) {
        q->a[i].dist = q->a[i].dist;
        q->a[i].vid = q->a[i].vid;
        q->a[i].silo_id = q->a[i].silo_id;
    }
}

static void OblivEnqueue(struct OblivQueue *q, float dist, int id, int silo_id) {
    if(q->ptr >= q->N) {
        if (dist >= q->a[q->ptr-1].dist) {
            return;
        } else {
            OblivDequeue(q);
        }
    }
    obliv_append(q, dist, id, silo_id);
    // q->a[q->ptr].dist = dist;
    // q->a[q->ptr].vid = id;
    // q->a[q->ptr].silo_id = silo_id;
    // q->ptr += 1;
    // obliv_resort(q);
}

static int GetDataIndex(struct OblivQueue *q, int id) {
    int ans = -1;
    for(int i=0; i<q->ptr; i++) {
        if (q->a[i].vid == id) {
            ans = i;
        }
    }
    return ans;
}

static int GetSiloIndex(struct OblivQueue *q, int silo_id) {
    int ans = -1;
    for(int i=0; i<q->ptr; i++) {
        if (q->a[i].silo_id == silo_id) {
            ans = i;
        }
    }
    return ans;
}

static void OblivFreeQueue(struct OblivQueue*& q) {
    if (q!=NULL && q->a!=NULL)
        free(q->a);
    q->a = NULL;
    if (q!=NULL && q->b!=NULL)
        free(q->b);
    q->b = NULL;
    if (q != NULL)
        free(q);
    q = NULL;
}

static void OblivCreateQueue(struct OblivQueue*& q, int queue_size) {
    q = (struct OblivQueue*) malloc (sizeof(struct OblivQueue));
    q->a = (struct QueueItem*) malloc(queue_size * sizeof(struct QueueItem));
    q->b = (struct QueueItem*) malloc(queue_size * sizeof(struct QueueItem));
    q->N = queue_size;
    obliv_init_queue(q, queue_size);
}

static void OblivInitQueue(struct OblivQueue* q, int queue_size) {
    obliv_init_queue(q, queue_size);
}

static void OblivQueueHead(struct OblivQueue *q, float* dist, int* vid, int* silo_id) {
    *dist = -1.0;
    *vid = -1;
    *silo_id = -1;
    for(int i=0; i<q->ptr; ++i) {
        *dist = q->a[i].dist;
        *vid = q->a[i].vid;
        *silo_id = q->a[i].silo_id;
    }
}

void ecall_OblivDequeue(void) {
    OblivDequeue(obliv_queue);
}

void ecall_OblivDequeueMore(int num) {
    OblivDequeue(obliv_queue, num);
}

void ecall_OblivCreateQueue(int queue_size) {
    OblivCreateQueue(obliv_queue, queue_size);
}

void ecall_OblivInitQueue(int queue_size) {
    OblivInitQueue(obliv_queue, queue_size);
}

int ecall_GetSiloIndex(int silo_id) {
    return GetSiloIndex(obliv_queue, silo_id);
}

int ecall_GetDataIndex(int data_id) {
    return GetDataIndex(obliv_queue, data_id);
}

void ecall_OblivFreeQueue(void) {
    OblivFreeQueue(obliv_queue);
    obliv_queue = NULL;
}

void ecall_OblivQueueHead(int* idx, int* silo_id) {
    float dist;
    OblivQueueHead(obliv_queue, &dist, idx, silo_id);
}

int ecall_Enqueue(int silo_id, int silo_beg_idx, int num,
                    const uint8_t *aes_key, const uint8_t* aes_iv,
                    const unsigned char* encrypt_data, size_t data_size) {

    int ret_status = 0;

    #ifdef LOCAL_SGX_DEBUG
    local_printf("[START] ecall_enqueue\n");
    local_printf("silo_id = %d, silo_beg_idx = %d, num = %d\n", silo_id, silo_beg_idx, num);
    #endif

    unsigned char* plain_data = NULL;
    ret_status = mbedtls_decrypt_data(aes_key, aes_iv, encrypt_data, data_size, num, plain_data);
    if (0 != ret_status) return ret_status;

    #ifdef LOCAL_SGX_DEBUG
    local_printf("[FINISH] AES decryption\n");
    {
        int batch_size = sizeof(int32_t) + sizeof(float);
        for (size_t j=0, offset=0; j<num; ++j, offset+=batch_size) {
            int vid = UnsignedVectorToInt32(plain_data, offset);
            float dis = UnsignedVectorToFloat(plain_data, offset+sizeof(int32_t));
            int idx = silo_beg_idx + j;
            local_printf("Decrypt data %d: silo_id = %d, dis = %.6f, vid = %d, idx = %d\n", j, silo_id, dis, vid, idx);
        }
    }
    #endif
    
    int batch_size = sizeof(int32_t) + sizeof(float);
    for (size_t j=0, offset=0; j<num; ++j, offset+=batch_size) {
        int vid = UnsignedVectorToInt32(plain_data, offset);
        float dis = UnsignedVectorToFloat(plain_data, offset+sizeof(int32_t));
        int idx = silo_beg_idx + j;
        OblivEnqueue(obliv_queue, dis, idx, silo_id);
    }

    #ifdef LOCAL_SGX_DEBUG
    local_printf("[FINISH] ecall_enqueue\n");
    #endif

    free(plain_data);

    return ret_status;
}

