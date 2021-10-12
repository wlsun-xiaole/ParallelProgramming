#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<sys/time.h>
#include<omp.h>
#include<math.h>

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

#define NUM_THREADS 8
// #define TOTAL_FUNC_CALL 49761
#define TOTAL_FUNC_CALL 4969077
char *str_islow = "islow function";

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

// feed 8*8 elements at a time
void func_islow(int *elemptr, DCTELEM * dataptr){
    INT32 tmp0, tmp1, tmp2, tmp3;
    INT32 tmp10, tmp11, tmp12, tmp13;
    INT32 z1;
    // DCTELEM *dataptr = malloc(DCTSIZE*DCTSIZE*sizeof(int));
    // DCTELEM *dataptr_start;
    // dataptr_start = dataptr;    // save the starting point of dataptr
    int ctr;

    // #pragma omp parallel for private(tmp0,tmp1,tmp2,tmp3,tmp10,tmp11,tmp12,tmp13,z1) num_threads(NUM_THREADS)
    #pragma omp parallel for private(tmp0,tmp1,tmp2,tmp3,tmp10,tmp11,tmp12,tmp13,z1)
    for(ctr = 0 ; ctr<DCTSIZE ; ctr++){
        int offset = ctr*DCTSIZE;

        tmp0 = GETJSAMPLE(elemptr[offset+0]) + GETJSAMPLE(elemptr[offset+7]);
        tmp1 = GETJSAMPLE(elemptr[offset+1]) + GETJSAMPLE(elemptr[offset+6]);
        tmp2 = GETJSAMPLE(elemptr[offset+2]) + GETJSAMPLE(elemptr[offset+5]);
        tmp3 = GETJSAMPLE(elemptr[offset+3]) + GETJSAMPLE(elemptr[offset+4]);

        tmp10 = tmp0 + tmp3;
        tmp12 = tmp0 - tmp3;
        tmp11 = tmp1 + tmp2;
        tmp13 = tmp1 - tmp2;

        tmp0 = GETJSAMPLE(elemptr[offset+0]) - GETJSAMPLE(elemptr[offset+7]);
        tmp1 = GETJSAMPLE(elemptr[offset+1]) - GETJSAMPLE(elemptr[offset+6]);
        tmp2 = GETJSAMPLE(elemptr[offset+2]) - GETJSAMPLE(elemptr[offset+5]);
        tmp3 = GETJSAMPLE(elemptr[offset+3]) - GETJSAMPLE(elemptr[offset+4]);

        dataptr[offset+0] = (DCTELEM) ((tmp10 + tmp11 - 8 * CENTERJSAMPLE) << PASS1_BITS);
        dataptr[offset+4] = (DCTELEM) ((tmp10 - tmp11) << PASS1_BITS);

        z1 = MULTIPLY(tmp12 + tmp13, FIX_0_541196100);       /* c6 */
        
        z1 += ONE << (CONST_BITS-PASS1_BITS-1);
        
        dataptr[offset+2] = (DCTELEM)
            RIGHT_SHIFT(z1 + MULTIPLY(tmp12, FIX_0_765366865), /* c2-c6 */
            CONST_BITS-PASS1_BITS);
        dataptr[offset+6] = (DCTELEM)
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

        dataptr[offset+1] = (DCTELEM) RIGHT_SHIFT(tmp0, CONST_BITS-PASS1_BITS);
        dataptr[offset+3] = (DCTELEM) RIGHT_SHIFT(tmp1, CONST_BITS-PASS1_BITS);
        dataptr[offset+5] = (DCTELEM) RIGHT_SHIFT(tmp2, CONST_BITS-PASS1_BITS);
        dataptr[offset+7] = (DCTELEM) RIGHT_SHIFT(tmp3, CONST_BITS-PASS1_BITS);
        
        // elemptr += DCTSIZE;
        // dataptr += DCTSIZE;
    }

    // dataptr = dataptr_start;
    // #pragma omp parallel for private(tmp0,tmp1,tmp2,tmp3,tmp10,tmp11,tmp12,tmp13,z1) num_threads(NUM_THREADS)
    #pragma omp parallel for private(tmp0,tmp1,tmp2,tmp3,tmp10,tmp11,tmp12,tmp13,z1)
    for(ctr = 0 ; ctr < DCTSIZE ; ctr++){
        tmp0 = dataptr[ctr+DCTSIZE*0] + dataptr[ctr+DCTSIZE*7];
        tmp1 = dataptr[ctr+DCTSIZE*1] + dataptr[ctr+DCTSIZE*6];
        tmp2 = dataptr[ctr+DCTSIZE*2] + dataptr[ctr+DCTSIZE*5];
        tmp3 = dataptr[ctr+DCTSIZE*3] + dataptr[ctr+DCTSIZE*4];

        /* Add fudge factor here for final descale. */
        tmp10 = tmp0 + tmp3 + (ONE << (PASS1_BITS-1));
        tmp12 = tmp0 - tmp3;
        tmp11 = tmp1 + tmp2;
        tmp13 = tmp1 - tmp2;

        tmp0 = dataptr[ctr+DCTSIZE*0] - dataptr[ctr+DCTSIZE*7];
        tmp1 = dataptr[ctr+DCTSIZE*1] - dataptr[ctr+DCTSIZE*6];
        tmp2 = dataptr[ctr+DCTSIZE*2] - dataptr[ctr+DCTSIZE*5];
        tmp3 = dataptr[ctr+DCTSIZE*3] - dataptr[ctr+DCTSIZE*4];

        dataptr[ctr+DCTSIZE*0] = (DCTELEM) RIGHT_SHIFT(tmp10 + tmp11, PASS1_BITS);
        dataptr[ctr+DCTSIZE*4] = (DCTELEM) RIGHT_SHIFT(tmp10 - tmp11, PASS1_BITS);

        z1 = MULTIPLY(tmp12 + tmp13, FIX_0_541196100);       /* c6 */
        /* Add fudge factor here for final descale. */
        z1 += ONE << (CONST_BITS+PASS1_BITS-1);

        dataptr[ctr+DCTSIZE*2] = (DCTELEM)
        RIGHT_SHIFT(z1 + MULTIPLY(tmp12, FIX_0_765366865), /* c2-c6 */
            CONST_BITS+PASS1_BITS);
        dataptr[ctr+DCTSIZE*6] = (DCTELEM)
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

        dataptr[ctr+DCTSIZE*1] = (DCTELEM) RIGHT_SHIFT(tmp0, CONST_BITS+PASS1_BITS);
        dataptr[ctr+DCTSIZE*3] = (DCTELEM) RIGHT_SHIFT(tmp1, CONST_BITS+PASS1_BITS);
        dataptr[ctr+DCTSIZE*5] = (DCTELEM) RIGHT_SHIFT(tmp2, CONST_BITS+PASS1_BITS);
        dataptr[ctr+DCTSIZE*7] = (DCTELEM) RIGHT_SHIFT(tmp3, CONST_BITS+PASS1_BITS);

        // dataptr++;
    }
}

void process_result(double *time_data){
    double total_time = 0;
    double mean;
    double std = 0;
    for(int i = 0 ; i<TOTAL_FUNC_CALL ; i++){
        total_time += time_data[i];
    }
    mean = total_time / TOTAL_FUNC_CALL;
    for(int i = 0 ; i<TOTAL_FUNC_CALL ; i++){
        std += pow(time_data[i] - mean, 2);
    }
    std = sqrt(std / TOTAL_FUNC_CALL);
    printf("total time: %lf sec\n", total_time);
    printf("std: %lf\n", std);
}

int main(void){

    freopen("./data_self_openMP.txt", "w", stderr);
    int *elem = malloc(TOTAL_FUNC_CALL*8*8*(sizeof(int)));
    int *dataptr = malloc(TOTAL_FUNC_CALL*8*8*(sizeof(int)));
    int *elem_start, *dataptr_start;
    double time_jpeg_fdct_islow = 0;
    double time_data;

    FILE *fptr;
    fptr = fopen("time_data.txt", "w");

    struct timeval start, end;
    readFile(elem);

    printf("start processing data\n");
    elem_start = elem;
    dataptr_start = dataptr;
    for(int i = 0 ; i<TOTAL_FUNC_CALL ; i++){
        gettimeofday(&start, NULL);
        func_islow(elem, dataptr);
        elem += DCTSIZE*DCTSIZE;
        dataptr += DCTSIZE*DCTSIZE;
        gettimeofday(&end, NULL);
        time_jpeg_fdct_islow += ((double) end.tv_sec - (double) start.tv_sec) + ((double) end.tv_usec - (double) start.tv_usec)/1000000.0;
        // time_data = ((double) end.tv_sec - (double) start.tv_sec) + ((double) end.tv_usec - (double) start.tv_usec)/1000000.0;   // unit is us
        fprintf(fptr, "%lf\n", time_data * 1000000);
    }
    printf("done processing data\n");
    
    fclose(fptr);

    printf("time used:%lfs\n", time_jpeg_fdct_islow);
    // printFile(elem_start, dataptr_start);


    return 0;
}
