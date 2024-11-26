#ifndef MPC_Oblivious_Queue_HPP
#define MPC_Oblivious_Queue_HPP

#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

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

void OblivDequeue(struct OblivQueue *q);
void OblivDequeue(struct OblivQueue *q, int num);
void OblivEnqueue(struct OblivQueue *q, float dist, int id, int silo_id);
void OblivQueueHead(struct OblivQueue *q, float* dist, int* id, int* silo_id);
int GetDataIndex(struct OblivQueue *q, int id);
int GetSiloIndex(struct OblivQueue *q, int silo_id);
void OblivFreeQueue(struct OblivQueue*& q);
void OblivCreateQueue(struct OblivQueue*& q, int queue_size);
void OblivInitQueue(struct OblivQueue* q, int queue_size);

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

void OblivDequeue(struct OblivQueue *q) {
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

void OblivDequeue(struct OblivQueue *q, int num) {
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

void OblivEnqueue(struct OblivQueue *q, float dist, int id, int silo_id) {
    // #ifdef LOCAL_DEBUG
    // std::cout << "[OblivEnqueue] START" << std::endl;
    // #endif
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

int GetDataIndex(struct OblivQueue *q, int id) {
    int ans = -1;
    for(int i=0; i<q->ptr; i++) {
        if (q->a[i].vid == id) {
            ans = i;
        }
    }
    return ans;
}

int GetSiloIndex(struct OblivQueue *q, int silo_id) {
    int ans = -1;
    for(int i=0; i<q->ptr; i++) {
        if (q->a[i].silo_id == silo_id) {
            ans = i;
        }
    }
    return ans;
}

void OblivFreeQueue(struct OblivQueue*& q) {
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

void OblivCreateQueue(struct OblivQueue*& q, int queue_size) {
    q = (struct OblivQueue*) malloc (sizeof(struct OblivQueue));
    q->a = (struct QueueItem*) malloc(queue_size * sizeof(struct QueueItem));
    q->b = (struct QueueItem*) malloc(queue_size * sizeof(struct QueueItem));
    q->N = queue_size;
    obliv_init_queue(q, queue_size);
}

void OblivInitQueue(struct OblivQueue* q, int queue_size) {
    obliv_init_queue(q, queue_size);
}

void OblivQueueHead(struct OblivQueue *q, float* dist, int* vid, int* silo_id) {
    *dist = -1.0;
    *vid = -1;
    *silo_id = -1;
    for(int i=0; i<q->ptr; ++i) {
        *dist = q->a[i].dist;
        *vid = q->a[i].vid;
        *silo_id = q->a[i].silo_id;
    }
}

#endif // MPC_Oblivious_Queue_HPP