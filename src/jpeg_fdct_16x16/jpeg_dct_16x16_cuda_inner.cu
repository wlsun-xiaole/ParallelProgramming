#define JPEG_INTERNALS
#include "jinclude.h"
#include "jpeglib.h"
#include "jdct.h"		/* Private declarations for DCT subsystem */
#include <fcntl.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include<sys/time.h>
#ifdef DCT_ISLOW_SUPPORTED


/*
 * This module is specialized to the case DCTSIZE = 8.
 */

#if DCTSIZE != 8
  Sorry, this code only copes with 8x8 DCT blocks. /* deliberate syntax err */
#endif


/*
 * The poop on this scaling stuff is as follows:
 *
 * Each 1-D DCT step produces outputs which are a factor of sqrt(N)
 * larger than the true DCT outputs.  The final outputs are therefore
 * a factor of N larger than desired; since N=8 this can be cured by
 * a simple right shift at the end of the algorithm.  The advantage of
 * this arrangement is that we save two multiplications per 1-D DCT,
 * because the y0 and y4 outputs need not be divided by sqrt(N).
 * In the IJG code, this factor of 8 is removed by the quantization step
 * (in jcdctmgr.c), NOT in this module.
 *
 * We have to do addition and subtraction of the integer inputs, which
 * is no problem, and multiplication by fractional constants, which is
 * a problem to do in integer arithmetic.  We multiply all the constants
 * by CONST_SCALE and convert them to integer constants (thus retaining
 * CONST_BITS bits of precision in the constants).  After doing a
 * multiplication we have to divide the product by CONST_SCALE, with proper
 * rounding, to produce the correct output.  This division can be done
 * cheaply as a right shift of CONST_BITS bits.  We postpone shifting
 * as long as possible so that partial sums can be added together with
 * full fractional precision.
 *
 * The outputs of the first pass are scaled up by PASS1_BITS bits so that
 * they are represented to better-than-integral precision.  These outputs
 * require BITS_IN_JSAMPLE + PASS1_BITS + 3 bits; this fits in a 16-bit word
 * with the recommended scaling.  (For 12-bit sample data, the intermediate
 * array is INT32 anyway.)
 *
 * To avoid overflow of the 32-bit intermediate results in pass 2, we must
 * have BITS_IN_JSAMPLE + CONST_BITS + PASS1_BITS <= 26.  Error analysis
 * shows that the values given below are the most effective.
 */

#if BITS_IN_JSAMPLE == 8
#define CONST_BITS  13
#define PASS1_BITS  2
#else
#define CONST_BITS  13
#define PASS1_BITS  1		/* lose a little precision to avoid overflow */
#endif

/* Some C compilers fail to reduce "FIX(constant)" at compile time, thus
 * causing a lot of useless floating-point operations at run time.
 * To get around this we use the following pre-calculated constants.
 * If you change CONST_BITS you may want to add appropriate values.
 * (With a reasonable C compiler, you can just rely on the FIX() macro...)
 */

#if CONST_BITS == 13
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
#else
#define FIX_0_298631336  FIX(0.298631336)
#define FIX_0_390180644  FIX(0.390180644)
#define FIX_0_541196100  FIX(0.541196100)
#define FIX_0_765366865  FIX(0.765366865)
#define FIX_0_899976223  FIX(0.899976223)
#define FIX_1_175875602  FIX(1.175875602)
#define FIX_1_501321110  FIX(1.501321110)
#define FIX_1_847759065  FIX(1.847759065)
#define FIX_1_961570560  FIX(1.961570560)
#define FIX_2_053119869  FIX(2.053119869)
#define FIX_2_562915447  FIX(2.562915447)
#define FIX_3_072711026  FIX(3.072711026)
#endif
#ifdef DCT_SCALING_SUPPORTED

/* Multiply an INT32 variable by an INT32 constant to yield an INT32 result.
 * For 8-bit samples with the recommended scaling, all the variable
 * and constant values involved are no more than 16 bits wide, so a
 * 16x16->32 bit multiply can be used instead of a full 32x32 multiply.
 * For 12-bit samples, a full 32-bit multiplication will be needed.
 */

#if BITS_IN_JSAMPLE == 8
#define MULTIPLY(var,const)  MULTIPLY16C16(var,const)
#else
#define MULTIPLY(var,const)  ((var) * (const))
#endif
struct  timeval start,end;
int *devicedataptr,*deviceelem,*devicedata;
DCTELEM *dataptr,*elem;
DCTELEM *workspace;
int threadnumber = 16;
int blocknumber = 1;
__global__ void jpeg_fdct_16x16_row (DCTELEM *elemptr_local, int *dataptr_local,DCTELEM *workspace)
{
  int ID = blockIdx.x * blockDim.x + threadIdx.x;
  INT32 tmp[18];
  

    /* Pass 1: process rows.
    * Note results are scaled up by sqrt(8) compared to a true DCT;
    * furthermore, we scale the results by 2**PASS1_BITS.
    * cK represents sqrt(2) * cos(K*pi/32).
    */
    elemptr_local += ID*16;
    if(ID >= 8)
      dataptr_local=workspace+(ID-8)*8;
    else
      dataptr_local+=ID*8;

    
   
    /* Even part */
   
    for(int i = 0;i<8;i++){
      tmp[i] = elemptr_local[i]+elemptr_local[15-i]; 
    }
    for( int i =0; i<4;i++){
      tmp[i+10] = tmp[i]+tmp[7-i];
      tmp[i+14] = tmp[i]-tmp[7-i];
    }
    for(int i =0;i<8;i++){
      tmp[i] = elemptr_local[i]-elemptr_local[15-i];
      
    }
    /* Apply unsigned->signed conversion */
    
        dataptr_local[0] = (DCTELEM)
        ((tmp[10] + tmp[11] + tmp[12] + tmp[13] - 16 * CENTERJSAMPLE) << PASS1_BITS);
        dataptr_local[4] = (DCTELEM)
        DESCALE(MULTIPLY(tmp[10] - tmp[13], FIX(1.306562965)) + /* c4[16] = c2[8] */
            MULTIPLY(tmp[11] - tmp[12], FIX_0_541196100),   /* c12[16] = c6[8] */
            CONST_BITS-PASS1_BITS);

        tmp[10] = MULTIPLY(tmp[17] - tmp[15], FIX(0.275899379)) +   /* c14[16] = c7[8] */
            MULTIPLY(tmp[14] - tmp[16], FIX(1.387039845));    /* c2[16] = c1[8] */

        dataptr_local[2] = (DCTELEM)
        DESCALE(tmp[10] + MULTIPLY(tmp[15], FIX(1.451774982))   /* c6+c14 */
            + MULTIPLY(tmp[16], FIX(2.172734804)),        /* c2+c10 */
            CONST_BITS-PASS1_BITS);
        dataptr_local[6] = (DCTELEM)
        DESCALE(tmp[10] - MULTIPLY(tmp[14], FIX(0.211164243))   /* c2-c6 */
            - MULTIPLY(tmp[17], FIX(1.061594338)),        /* c10+c14 */
            CONST_BITS-PASS1_BITS);

    /* Odd part */

        tmp[11] = MULTIPLY(tmp[0] + tmp[1], FIX(1.353318001)) +         /* c3 */
            MULTIPLY(tmp[6] - tmp[7], FIX(0.410524528));          /* c13 */
        tmp[12] = MULTIPLY(tmp[0] + tmp[2], FIX(1.247225013)) +         /* c5 */
            MULTIPLY(tmp[5] + tmp[7], FIX(0.666655658));          /* c11 */
        tmp[13] = MULTIPLY(tmp[0] + tmp[3], FIX(1.093201867)) +         /* c7 */
            MULTIPLY(tmp[4] - tmp[7], FIX(0.897167586));          /* c9 */
        tmp[14] = MULTIPLY(tmp[1] + tmp[2], FIX(0.138617169)) +         /* c15 */
            MULTIPLY(tmp[6] - tmp[5], FIX(1.407403738));          /* c1 */
        tmp[15] = MULTIPLY(tmp[1] + tmp[3], - FIX(0.666655658)) +       /* -c11 */
            MULTIPLY(tmp[4] + tmp[6], - FIX(1.247225013));        /* -c5 */
        tmp[16] = MULTIPLY(tmp[2] + tmp[3], - FIX(1.353318001)) +       /* -c3 */
            MULTIPLY(tmp[5] - tmp[4], FIX(0.410524528));          /* c13 */
        tmp[10] = tmp[11] + tmp[12] + tmp[13] -
            MULTIPLY(tmp[0], FIX(2.286341144)) +                /* c7+c5+c3-c1 */
            MULTIPLY(tmp[7], FIX(0.779653625));                 /* c15+c13-c11+c9 */
        tmp[11] += tmp[14] + tmp[15] + MULTIPLY(tmp[1], FIX(0.071888074)) /* c9-c3-c15+c11 */
            - MULTIPLY(tmp[6], FIX(1.663905119));              /* c7+c13+c1-c5 */
        tmp[12] += tmp[14] + tmp[16] - MULTIPLY(tmp[2], FIX(1.125726048)) /* c7+c5+c15-c3 */
            + MULTIPLY(tmp[5], FIX(1.227391138));              /* c9-c11+c1-c13 */
        tmp[13] += tmp[15] + tmp[16] + MULTIPLY(tmp[3], FIX(1.065388962)) /* c15+c3+c11-c7 */
            + MULTIPLY(tmp[4], FIX(2.167985692));              /* c1+c13+c5-c9 */

        dataptr_local[1] = (DCTELEM) DESCALE(tmp[10], CONST_BITS-PASS1_BITS);
        dataptr_local[3] = (DCTELEM) DESCALE(tmp[11], CONST_BITS-PASS1_BITS);
        dataptr_local[5] = (DCTELEM) DESCALE(tmp[12], CONST_BITS-PASS1_BITS);
        dataptr_local[7] = (DCTELEM) DESCALE(tmp[13], CONST_BITS-PASS1_BITS);

     
    
    //if(ID%8==7 && ((ID+1)/8)%2==1)
    //    dataptr_local = workspace;	/* switch pointer to extended workspace */


}
__global__ void jpeg_fdct_16x16_col(int *dataptr_local, DCTELEM *workspace)
{
   /* Pass 2: process columns.
   * We remove the PASS1_BITS s
   * caling, but leave the results scaled up
   * by an overall factor of 8.
   * We must also scale the output by (8/16)**2 = 1/2**2.
   * cK represents sqrt(2) * cos(K*pi/32).
   */

   int ID = blockIdx.x * blockDim.x + threadIdx.x;
   INT32 tmp[18];
   DCTELEM *wsptr;

    dataptr_local += ID; 
    wsptr = workspace + ID;

    /* Even part */

    for(int i=0;i<8;i++){
      tmp[i] = dataptr_local[DCTSIZE*i]+ wsptr[DCTSIZE*(7-i)];
    }

    for(int i=0;i<4;i++){
      tmp[i+10] = tmp[i]+tmp[7-i];
      tmp[i+14] = tmp[i]-tmp[7-i];
    }

    for(int i=0;i<8;i++)
      tmp[i] = dataptr_local[DCTSIZE*i] - wsptr[DCTSIZE*(7-i)];


    dataptr_local[DCTSIZE*0] = (DCTELEM)
      DESCALE(tmp[10] + tmp[11] + tmp[12] + tmp[13], PASS1_BITS+2);
    dataptr_local[DCTSIZE*4] = (DCTELEM)
      DESCALE(MULTIPLY(tmp[10] - tmp[13], FIX(1.306562965)) + /* c4[16] = c2[8] */
	      MULTIPLY(tmp[11] - tmp[12], FIX_0_541196100),   /* c12[16] = c6[8] */
	      CONST_BITS+PASS1_BITS+2);

    tmp[10] = MULTIPLY(tmp[17] - tmp[15], FIX(0.275899379)) +   /* c14[16] = c7[8] */
	    MULTIPLY(tmp[14] - tmp[16], FIX(1.387039845));    /* c2[16] = c1[8] */

    dataptr_local[DCTSIZE*2] = (DCTELEM)
      DESCALE(tmp[10] + MULTIPLY(tmp[15], FIX(1.451774982))   /* c6+c14 */
	      + MULTIPLY(tmp[16], FIX(2.172734804)),        /* c2+10 */
	      CONST_BITS+PASS1_BITS+2);
    dataptr_local[DCTSIZE*6] = (DCTELEM)
      DESCALE(tmp[10] - MULTIPLY(tmp[14], FIX(0.211164243))   /* c2-c6 */
	      - MULTIPLY(tmp[17], FIX(1.061594338)),        /* c10+c14 */
	      CONST_BITS+PASS1_BITS+2);

    /* Odd part */

    tmp[11] = MULTIPLY(tmp[0]+ tmp[1], FIX(1.353318001)) +         /* c3 */
	    MULTIPLY(tmp[6] - tmp[7], FIX(0.410524528));          /* c13 */
    tmp[12] = MULTIPLY(tmp[0] + tmp[2], FIX(1.247225013)) +         /* c5 */
	    MULTIPLY(tmp[5] + tmp[7], FIX(0.666655658));          /* c11 */
    tmp[13] = MULTIPLY(tmp[0] + tmp[3], FIX(1.093201867)) +         /* c7 */
	    MULTIPLY(tmp[4] - tmp[7], FIX(0.897167586));          /* c9 */
    tmp[14] = MULTIPLY(tmp[1] + tmp[2], FIX(0.138617169)) +         /* c15 */
	    MULTIPLY(tmp[6] - tmp[5], FIX(1.407403738));          /* c1 */
    tmp[15] = MULTIPLY(tmp[1] + tmp[3], - FIX(0.666655658)) +       /* -c11 */
	    MULTIPLY(tmp[4] + tmp[6], - FIX(1.247225013));        /* -c5 */
    tmp[16] = MULTIPLY(tmp[2] + tmp[3], - FIX(1.353318001)) +       /* -c3 */
	    MULTIPLY(tmp[5] - tmp[4], FIX(0.410524528));          /* c13 */
    tmp[10] = tmp[11] + tmp[12] + tmp[13] -
	    MULTIPLY(tmp[0], FIX(2.286341144)) +                /* c7+c5+c3-c1 */
	    MULTIPLY(tmp[7], FIX(0.779653625));                 /* c15+c13-c11+c9 */
    tmp[11] += tmp[14] + tmp[15] + MULTIPLY(tmp[1], FIX(0.071888074)) /* c9-c3-c15+c11 */
	     - MULTIPLY(tmp[6], FIX(1.663905119));              /* c7+c13+c1-c5 */
    tmp[12] += tmp[14] + tmp[16] - MULTIPLY(tmp[2], FIX(1.125726048)) /* c7+c5+c15-c3 */
	     + MULTIPLY(tmp[5], FIX(1.227391138));              /* c9-c11+c1-c13 */
    tmp[13] += tmp[15] + tmp[16] + MULTIPLY(tmp[3], FIX(1.065388962)) /* c15+c3+c11-c7 */
	     + MULTIPLY(tmp[4], FIX(2.167985692));              /* c1+c13+c5-c9 */

    dataptr_local[DCTSIZE*1] = (DCTELEM) DESCALE(tmp[10], CONST_BITS+PASS1_BITS+2);
    dataptr_local[DCTSIZE*3] = (DCTELEM) DESCALE(tmp[11], CONST_BITS+PASS1_BITS+2);
    dataptr_local[DCTSIZE*5] = (DCTELEM) DESCALE(tmp[12], CONST_BITS+PASS1_BITS+2);
    dataptr_local[DCTSIZE*7] = (DCTELEM) DESCALE(tmp[13], CONST_BITS+PASS1_BITS+2);
    
}
int main()
{
  FILE *sc,*fds;
 
  int *elem_start, *dataptr_start;
  int iter = 0;
  char check;
  unsigned int col;
  
  double  time_used = 0.0;
  fds = fopen("sample_data_l.txt","r");
  sc = fopen("start_col_l.txt","r");
  if(fds==NULL ||  sc ==NULL){
    printf("open error\n");
    return 0;
  }
  
  do{
    fscanf(sc,"%u",&col);
    
    check = fgetc(sc);
    iter++;
  }while(check!=EOF);
  fclose(sc);
  printf("%d\n",iter);
  
  dataptr = (int*) malloc(iter*256*sizeof(int));
  elem= (int*) malloc(iter*256*sizeof(int));

  

  for(int i=0;i<iter*256;i++){
    fscanf(fds, "%d",&elem[i]);
  }
  fclose(fds);

  cudaMalloc(&devicedataptr,16*16*sizeof(int));
  cudaMalloc(&deviceelem,16*16*sizeof(int));
  cudaMalloc(&workspace,64*sizeof(int));

  dataptr_start = dataptr;
  elem_start = elem;
  
  for(int i=0;i<iter;i++){
    gettimeofday(&start,NULL);
    cudaMemcpy(deviceelem, elem ,256*sizeof(int), cudaMemcpyHostToDevice);
    jpeg_fdct_16x16_row<<<blocknumber,threadnumber>>>(deviceelem,devicedataptr,workspace);
    jpeg_fdct_16x16_col<<<blocknumber,threadnumber>>>(devicedataptr,workspace);
    cudaMemcpy(dataptr, devicedataptr, 256*sizeof(int), cudaMemcpyDeviceToHost);
    gettimeofday(&end,NULL);
    elem += 256;
    dataptr += 256;
    time_used += (double)(end.tv_sec-start.tv_sec);
  time_used += (double)(end.tv_usec-start.tv_usec)/1000000;
  }
  
  
  fprintf(stderr,"time used:%lfs\n", time_used);
  /*for(int ctr = 0;ctr<iter;ctr++){
    if(ctr >21246&&ctr<22000){
      for(int i=0;i<57;i++){
        fprintf(stderr,"%d ",(int)dataptr_start[i]);
          
      }
      printf("%d\n",ctr);
    }
    dataptr_start+=256;
  }
 */
  return 0;
  
}



#endif /* DCT_SCALING_SUPPORTED */
#endif /* DCT_ISLOW_SUPPORTED */
