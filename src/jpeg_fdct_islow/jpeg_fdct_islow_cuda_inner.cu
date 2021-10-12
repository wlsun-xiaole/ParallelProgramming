#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#define DCTSIZE 8
#define CENTERJSAMPLE 128
#define PASS1_BITS 2
#define CONST_BITS 13
#define ONE	((INT32) 1)

#define FIX_0_298631336  ((INT32)  2446)	/* FIX(0.298631336) */
#define FIX_0_390180644  ((INT32)  3196)	/* FIX(0.390180644) */
#define FIX_0_541196100  ((INT32)  4433)	/* FIX(0.541196100) */
#define FIX_0_765366865  ((INT32)  6270)	/* FIX(0.765366865) */
#define FIX_0_899976223  ((INT32)  7373)	/* FIX(0.899976223) */
#define FIX_1_175875602  ((INT32)  9633)	/* FIX(1.175875602) */
#define FIX_1_501321110  ((INT32)  12299)	/* FIX(1.501321110) */
#define FIX_1_847759065  ((INT32)  15137)	/* FIX(1.847759065) */
#define FIX_1_961570560  ((INT32)  16069)	/* FIX(1.961570560) */
#define FIX_2_053119869  ((INT32)  16819)	/* FIX(2.053119869) */
#define FIX_2_562915447  ((INT32)  20995)	/* FIX(2.562915447) */
#define FIX_3_072711026  ((INT32)  25172)	/* FIX(3.072711026) */

#define RIGHT_SHIFT(x,shft)	((x) >> (shft))
#define GETJSAMPLE(value)  ((int) (value))
#define MULTIPLY16C16(var,const)  ((var) * (const))
#define MULTIPLY(var,const)  MULTIPLY16C16(var,const)

typedef int DCTELEM;
typedef long INT32;


// #define TOTAL_FUNC_CALL 49761
#define TOTAL_FUNC_CALL 4969077
char *str_islow = "islow function";
int *elemd, *dataptrd;
int numOfThreads = 4;
int numOfBlocks = 2;
int batchSize = DCTSIZE*DCTSIZE*sizeof(int);

void readFile(int *elem){
    FILE *fp;
    char *line;
    size_t len = 0;
    int ctr = 0, iter = 0;
    int *elemptr;
    // fp = fopen("./data_small.txt", "r");
    fp = fopen("./data_large.txt", "r");

    printf("start reading file...\n");
    elemptr = elem;
    while(getline(&line, &len, fp) != -1){
        if(strstr(line, str_islow) != NULL){
            iter++;
            if(iter % 1000 == 0)
                printf("%s", line);
            ctr = DCTSIZE;
        }
        else{
            if(ctr != 0){
                ctr--;
                sscanf(line, "elemptr[0]=%d,elemptr[1]=%d,elemptr[2]=%d,elemptr[3]=%d,elemptr[4]=%d,elemptr[5]=%d,elemptr[6]=%d,elemptr[7]=%d\n",
                        &elemptr[0], &elemptr[1], &elemptr[2], &elemptr[3], &elemptr[4], &elemptr[5], &elemptr[6], &elemptr[7]);
                elemptr += DCTSIZE;     // jump 8 elements
            }
            else{
                continue;
            }
        }
    }
    printf("done reading file\n");
}

void printFile(int *elemptr, DCTELEM *dataptr){
    int *e_p, *d_p;
    for(int i = 0 ; i<TOTAL_FUNC_CALL ; i++){
        fprintf(stderr, "%dth islow function\n", i+1);
        int ctr;
        e_p = elemptr;
        for(ctr = 0 ; ctr<DCTSIZE ; ctr++){
            fprintf(stderr, "elemptr[0]=%d,", (int)e_p[0]);
            fprintf(stderr, "elemptr[1]=%d,", (int)e_p[1]);
            fprintf(stderr, "elemptr[2]=%d,", (int)e_p[2]);
            fprintf(stderr, "elemptr[3]=%d,", (int)e_p[3]);
            fprintf(stderr, "elemptr[4]=%d,", (int)e_p[4]);
            fprintf(stderr, "elemptr[5]=%d,", (int)e_p[5]);
            fprintf(stderr, "elemptr[6]=%d,", (int)e_p[6]);
            fprintf(stderr, "elemptr[7]=%d\n", (int)e_p[7]);

            e_p += DCTSIZE;
        }
        d_p = dataptr;
        for(ctr = 0 ; ctr<DCTSIZE ; ctr++){
            fprintf(stderr, "dataptr[0]=%d,", d_p[DCTSIZE*0]);
            fprintf(stderr, "dataptr[1]=%d,", d_p[DCTSIZE*1]);
            fprintf(stderr, "dataptr[2]=%d,", d_p[DCTSIZE*2]);
            fprintf(stderr, "dataptr[3]=%d,", d_p[DCTSIZE*3]);
            fprintf(stderr, "dataptr[4]=%d,", d_p[DCTSIZE*4]);
            fprintf(stderr, "dataptr[5]=%d,", d_p[DCTSIZE*5]);
            fprintf(stderr, "dataptr[6]=%d,", d_p[DCTSIZE*6]);
            fprintf(stderr, "dataptr[7]=%d\n", d_p[DCTSIZE*7]);

            d_p += 1;
        }
        elemptr += DCTSIZE*DCTSIZE;
        dataptr += DCTSIZE*DCTSIZE;
    }
}

__global__ void islow_cuda_1(int *elemptrd, int *dataptrd){
    INT32 tmp0, tmp1, tmp2, tmp3;
    INT32 tmp10, tmp11, tmp12, tmp13;
    INT32 z1;

    int ID = blockIdx.x * blockDim.x + threadIdx.x;
    elemptrd += ID * DCTSIZE;
    dataptrd += ID * DCTSIZE;

    tmp0 = GETJSAMPLE(elemptrd[0]) + GETJSAMPLE(elemptrd[7]);
    tmp1 = GETJSAMPLE(elemptrd[1]) + GETJSAMPLE(elemptrd[6]);
    tmp2 = GETJSAMPLE(elemptrd[2]) + GETJSAMPLE(elemptrd[5]);
    tmp3 = GETJSAMPLE(elemptrd[3]) + GETJSAMPLE(elemptrd[4]);

    tmp10 = tmp0 + tmp3;
    tmp12 = tmp0 - tmp3;
    tmp11 = tmp1 + tmp2;
    tmp13 = tmp1 - tmp2;

    tmp0 = GETJSAMPLE(elemptrd[0]) - GETJSAMPLE(elemptrd[7]);
    tmp1 = GETJSAMPLE(elemptrd[1]) - GETJSAMPLE(elemptrd[6]);
    tmp2 = GETJSAMPLE(elemptrd[2]) - GETJSAMPLE(elemptrd[5]);
    tmp3 = GETJSAMPLE(elemptrd[3]) - GETJSAMPLE(elemptrd[4]);

    dataptrd[0] = (DCTELEM) ((tmp10 + tmp11 - 8 * CENTERJSAMPLE) << PASS1_BITS);
    dataptrd[4] = (DCTELEM) ((tmp10 - tmp11) << PASS1_BITS);

    z1 = MULTIPLY(tmp12 + tmp13, FIX_0_541196100);       /* c6 */
    
    z1 += ONE << (CONST_BITS-PASS1_BITS-1);
    
    dataptrd[2] = (DCTELEM)
        RIGHT_SHIFT(z1 + MULTIPLY(tmp12, FIX_0_765366865), /* c2-c6 */
        CONST_BITS-PASS1_BITS);
    dataptrd[6] = (DCTELEM)
        RIGHT_SHIFT(z1 - MULTIPLY(tmp13, FIX_1_847759065), /* c2+c6 */
        CONST_BITS-PASS1_BITS);
    
    tmp12 = tmp0 + tmp2;
    tmp13 = tmp1 + tmp3;

    z1 = MULTIPLY(tmp12 + tmp13, FIX_1_175875602);       /*  c3 */
    /* Add fudge factor here for final descale. */
    z1 += ONE << (CONST_BITS-PASS1_BITS-1);

    tmp12 = MULTIPLY(tmp12, - FIX_0_390180644);          /* -c3+c5 */
    tmp13 = MULTIPLY(tmp13, - FIX_1_961570560);          /* -c3-c5 */
    tmp12 += z1;
    tmp13 += z1;

    z1 = MULTIPLY(tmp0 + tmp3, - FIX_0_899976223);       /* -c3+c7 */
    tmp0 = MULTIPLY(tmp0, FIX_1_501321110);              /*  c1+c3-c5-c7 */
    tmp3 = MULTIPLY(tmp3, FIX_0_298631336);              /* -c1+c3+c5-c7 */
    tmp0 += z1 + tmp12;
    tmp3 += z1 + tmp13;

    z1 = MULTIPLY(tmp1 + tmp2, - FIX_2_562915447);       /* -c1-c3 */
    tmp1 = MULTIPLY(tmp1, FIX_3_072711026);              /*  c1+c3+c5-c7 */
    tmp2 = MULTIPLY(tmp2, FIX_2_053119869);              /*  c1+c3-c5+c7 */
    tmp1 += z1 + tmp13;
    tmp2 += z1 + tmp12;

    dataptrd[1] = (DCTELEM) RIGHT_SHIFT(tmp0, CONST_BITS-PASS1_BITS);
    dataptrd[3] = (DCTELEM) RIGHT_SHIFT(tmp1, CONST_BITS-PASS1_BITS);
    dataptrd[5] = (DCTELEM) RIGHT_SHIFT(tmp2, CONST_BITS-PASS1_BITS);
    dataptrd[7] = (DCTELEM) RIGHT_SHIFT(tmp3, CONST_BITS-PASS1_BITS);
}

__global__ void islow_cuda_2(int *dataptrd){
    INT32 tmp0, tmp1, tmp2, tmp3;
    INT32 tmp10, tmp11, tmp12, tmp13;
    INT32 z1;

    int ID = blockIdx.x * blockDim.x + threadIdx.x;
    dataptrd += ID;

    tmp0 = dataptrd[DCTSIZE*0] + dataptrd[DCTSIZE*7];
    tmp1 = dataptrd[DCTSIZE*1] + dataptrd[DCTSIZE*6];
    tmp2 = dataptrd[DCTSIZE*2] + dataptrd[DCTSIZE*5];
    tmp3 = dataptrd[DCTSIZE*3] + dataptrd[DCTSIZE*4];

    /* Add fudge factor here for final descale. */
    tmp10 = tmp0 + tmp3 + (ONE << (PASS1_BITS-1));
    tmp12 = tmp0 - tmp3;
    tmp11 = tmp1 + tmp2;
    tmp13 = tmp1 - tmp2;

    tmp0 = dataptrd[DCTSIZE*0] - dataptrd[DCTSIZE*7];
    tmp1 = dataptrd[DCTSIZE*1] - dataptrd[DCTSIZE*6];
    tmp2 = dataptrd[DCTSIZE*2] - dataptrd[DCTSIZE*5];
    tmp3 = dataptrd[DCTSIZE*3] - dataptrd[DCTSIZE*4];

    dataptrd[DCTSIZE*0] = (DCTELEM) RIGHT_SHIFT(tmp10 + tmp11, PASS1_BITS);
    dataptrd[DCTSIZE*4] = (DCTELEM) RIGHT_SHIFT(tmp10 - tmp11, PASS1_BITS);

    z1 = MULTIPLY(tmp12 + tmp13, FIX_0_541196100);       /* c6 */
    /* Add fudge factor here for final descale. */
    z1 += ONE << (CONST_BITS+PASS1_BITS-1);

    dataptrd[DCTSIZE*2] = (DCTELEM)
    RIGHT_SHIFT(z1 + MULTIPLY(tmp12, FIX_0_765366865), /* c2-c6 */
        CONST_BITS+PASS1_BITS);
    dataptrd[DCTSIZE*6] = (DCTELEM)
    RIGHT_SHIFT(z1 - MULTIPLY(tmp13, FIX_1_847759065), /* c2+c6 */
        CONST_BITS+PASS1_BITS);

    tmp12 = tmp0 + tmp2;
    tmp13 = tmp1 + tmp3;

    z1 = MULTIPLY(tmp12 + tmp13, FIX_1_175875602);       /*  c3 */
    /* Add fudge factor here for final descale. */
    z1 += ONE << (CONST_BITS+PASS1_BITS-1);

    tmp12 = MULTIPLY(tmp12, - FIX_0_390180644);          /* -c3+c5 */
    tmp13 = MULTIPLY(tmp13, - FIX_1_961570560);          /* -c3-c5 */
    tmp12 += z1;
    tmp13 += z1;

    z1 = MULTIPLY(tmp0 + tmp3, - FIX_0_899976223);       /* -c3+c7 */
    tmp0 = MULTIPLY(tmp0, FIX_1_501321110);              /*  c1+c3-c5-c7 */
    tmp3 = MULTIPLY(tmp3, FIX_0_298631336);              /* -c1+c3+c5-c7 */
    tmp0 += z1 + tmp12;
    tmp3 += z1 + tmp13;

    z1 = MULTIPLY(tmp1 + tmp2, - FIX_2_562915447);       /* -c1-c3 */
    tmp1 = MULTIPLY(tmp1, FIX_3_072711026);              /*  c1+c3+c5-c7 */
    tmp2 = MULTIPLY(tmp2, FIX_2_053119869);              /*  c1+c3-c5+c7 */
    tmp1 += z1 + tmp13;
    tmp2 += z1 + tmp12;

    dataptrd[DCTSIZE*1] = (DCTELEM) RIGHT_SHIFT(tmp0, CONST_BITS+PASS1_BITS);
    dataptrd[DCTSIZE*3] = (DCTELEM) RIGHT_SHIFT(tmp1, CONST_BITS+PASS1_BITS);
    dataptrd[DCTSIZE*5] = (DCTELEM) RIGHT_SHIFT(tmp2, CONST_BITS+PASS1_BITS);
    dataptrd[DCTSIZE*7] = (DCTELEM) RIGHT_SHIFT(tmp3, CONST_BITS+PASS1_BITS);
}

// feed 8*8 elements at a time
void func_islow(int *elemptr, DCTELEM * dataptr){
    cudaMemcpy(elemd, elemptr, batchSize, cudaMemcpyHostToDevice);

    islow_cuda_1<<<numOfBlocks, numOfThreads>>>(elemd, dataptrd);
    islow_cuda_2<<<numOfBlocks, numOfThreads>>>(dataptrd);

    cudaMemcpy(dataptr, dataptrd, batchSize, cudaMemcpyDeviceToHost);
}



int main(void){

    freopen("./data_self_cuda.txt", "w", stderr);
    int size = TOTAL_FUNC_CALL*DCTSIZE*DCTSIZE*(sizeof(int));
    int *elem = (int *)malloc(size);
    int *dataptr = (int *)malloc(size);
    int *elem_start, *dataptr_start;

    clock_t start, end;
    double time_used = 0;

    readFile(elem);

    cudaMalloc(&elemd, batchSize);
    cudaMalloc(&dataptrd, batchSize);

    elem_start = elem;
    dataptr_start = dataptr;
    printf("start processing data\n");
    start = clock();
    for(int i = 0 ; i<TOTAL_FUNC_CALL ; i++){
        func_islow(elem, dataptr);
        elem += DCTSIZE*DCTSIZE;
        dataptr += DCTSIZE*DCTSIZE;
    }
    end = clock();
    time_used += ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("done processing data\n");
    printf("time used:%lfs\n", time_used);

    printFile(elem_start, dataptr_start);


    return 0;
}
