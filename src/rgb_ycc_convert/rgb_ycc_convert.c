#include<stdio.h>
#include<stdlib.h>
#include<sys/time.h>

#ifdef MP_IN
#include<omp.h>
#endif 


typedef long INT32;
typedef unsigned int JDIMENSION;
#define SIZEOF(object) ((size_t) sizeof(object))
#define FAR far

typedef unsigned char JSAMPLE;

typedef JSAMPLE *JSAMPROW;
typedef JSAMPROW *JSAMPARRAY;
typedef JSAMPARRAY *JSAMPIMAGE;

#define GETJSAMPLE(value) ((int) (value))
#define MAXJSAMPLE 255
#define CENTERJSAMPLE 128


#define SCALEBITS	16	
#define CBCR_OFFSET	((INT32) CENTERJSAMPLE << SCALEBITS)
#define ONE_HALF	((INT32) 1 << (SCALEBITS-1))
#define FIX(x)		((INT32) ((x) * (1L<<SCALEBITS) + 0.5))

#define R_Y_OFF		0			
#define G_Y_OFF		(1*(MAXJSAMPLE+1))	
#define B_Y_OFF		(2*(MAXJSAMPLE+1))	
#define R_CB_OFF	(3*(MAXJSAMPLE+1))
#define G_CB_OFF	(4*(MAXJSAMPLE+1))
#define B_CB_OFF	(5*(MAXJSAMPLE+1))
#define R_CR_OFF	B_CB_OFF		
#define G_CR_OFF	(6*(MAXJSAMPLE+1))
#define B_CR_OFF	(7*(MAXJSAMPLE+1))
#define TABLE_SIZE	(8*(MAXJSAMPLE+1))

#define RGB_RED 0
#define RGB_GREEN 1
#define RGB_BLUE 2
#define RGB_PIXELSIZE 3




void rgb_ycc_start(INT32 *rgb_ycc_tab) {
	INT32 i;
	for (i = 0; i <= MAXJSAMPLE; i++) {
		rgb_ycc_tab[i+R_Y_OFF] = FIX(0.299) * i;
		rgb_ycc_tab[i+G_Y_OFF] = FIX(0.587) * i;
		rgb_ycc_tab[i+B_Y_OFF] = FIX(0.114) * i   + ONE_HALF;
		rgb_ycc_tab[i+R_CB_OFF] = (-FIX(0.168735892)) * i;
		rgb_ycc_tab[i+G_CB_OFF] = (-FIX(0.331264108)) * i;
		/* We use a rounding fudge-factor of 0.5-epsilon for Cb and Cr.
		 * This ensures that the maximum output will round to MAXJSAMPLE
		 * not MAXJSAMPLE+1, and thus that we don't have to range-limit.
		 */
		rgb_ycc_tab[i+B_CB_OFF] = FIX(0.5) * i    + CBCR_OFFSET + ONE_HALF-1;
	/*  B=>Cb and R=>Cr tables are the same
		rgb_ycc_tab[i+R_CR_OFF] = FIX(0.5) * i    + CBCR_OFFSET + ONE_HALF-1;
	*/
		rgb_ycc_tab[i+G_CR_OFF] = (-FIX(0.418687589)) * i;
		rgb_ycc_tab[i+B_CR_OFF] = (-FIX(0.081312411)) * i;
	}
}

#ifdef ORIGINAL
void rgb_ycc_convert (INT32 *rgb_ycc_tab, JSAMPARRAY input_buf, JSAMPIMAGE output_buf,
		 JDIMENSION output_row, int num_rows, JDIMENSION image_width)
{
  register INT32 * ctab = rgb_ycc_tab;
  register int r, g, b;
  register JSAMPROW inptr;
  register JSAMPROW outptr0, outptr1, outptr2;
  register JDIMENSION col;
  JDIMENSION num_cols = image_width;

  while (--num_rows >= 0) {
	
    outptr0 = output_buf[0][output_row];
    outptr1 = output_buf[1][output_row];
    outptr2 = output_buf[2][output_row];
    for (col = 0; col < num_cols; col++) {
	  inptr = input_buf[col];
      r = GETJSAMPLE(inptr[RGB_RED]);
      g = GETJSAMPLE(inptr[RGB_GREEN]);
      b = GETJSAMPLE(inptr[RGB_BLUE]);
      // If the inputs are 0..MAXJSAMPLE, the outputs of these equations
      //  must be too; we do not need an explicit range-limiting operation.
      //  Hence the value being shifted is never negative, and we don't
      //  need the general RIGHT_SHIFT macro.
       
      // Y 
      outptr0[col] = (JSAMPLE)
		((ctab[r+R_Y_OFF] + ctab[g+G_Y_OFF] + ctab[b+B_Y_OFF])
		 >> SCALEBITS);
      // Cb 
      outptr1[col] = (JSAMPLE)
		((ctab[r+R_CB_OFF] + ctab[g+G_CB_OFF] + ctab[b+B_CB_OFF])
		 >> SCALEBITS);
      // Cr 
      outptr2[col] = (JSAMPLE)
		((ctab[r+R_CR_OFF] + ctab[g+G_CR_OFF] + ctab[b+B_CR_OFF])
		 >> SCALEBITS);

	  
/*	  if(col == 1) {
		  printf("in: %d %d %d\n", inptr[RGB_RED], inptr[RGB_GREEN], inptr[RGB_BLUE]);
		  printf("out1: %d %d %d\n",output_buf[0][0][col],output_buf[1][0][col],output_buf[2][0][col]);
		  printf("out2: %d %d %d\n",outptr0[col],outptr1[col],outptr2[col]);
	  }
*/	  
	  //inptr += RGB_PIXELSIZE;
    }
  }
}
#endif


#ifdef MP_IN
void rgb_ycc_convert (INT32 *rgb_ycc_tab, JSAMPARRAY input_buf, JSAMPIMAGE output_buf,
		 JDIMENSION output_row, int num_rows, JDIMENSION image_width)
{
  register INT32 * ctab = rgb_ycc_tab;
  register int r, g, b;
  register JSAMPROW inptr;
  register JSAMPROW outptr0, outptr1, outptr2;
  register JDIMENSION col;
  JDIMENSION num_cols = image_width;

  while (--num_rows >= 0) {
	
    outptr0 = output_buf[0][output_row];
    outptr1 = output_buf[1][output_row];
    outptr2 = output_buf[2][output_row];
#pragma omp parallel firstprivate(outptr0,outptr1,outptr2)
#pragma omp for private(col,r,g,b,inptr)
    for (col = 0; col < num_cols; col++) {
	  inptr = input_buf[col];
      r = GETJSAMPLE(inptr[RGB_RED]);
      g = GETJSAMPLE(inptr[RGB_GREEN]);
      b = GETJSAMPLE(inptr[RGB_BLUE]);
      // If the inputs are 0..MAXJSAMPLE, the outputs of these equations
      //  must be too; we do not need an explicit range-limiting operation.
      //  Hence the value being shifted is never negative, and we don't
      //  need the general RIGHT_SHIFT macro.
       
      // Y 
      outptr0[col] = (JSAMPLE)
		((ctab[r+R_Y_OFF] + ctab[g+G_Y_OFF] + ctab[b+B_Y_OFF])
		 >> SCALEBITS);
      // Cb 
      outptr1[col] = (JSAMPLE)
		((ctab[r+R_CB_OFF] + ctab[g+G_CB_OFF] + ctab[b+B_CB_OFF])
		 >> SCALEBITS);
      // Cr 
      outptr2[col] = (JSAMPLE)
		((ctab[r+R_CR_OFF] + ctab[g+G_CR_OFF] + ctab[b+B_CR_OFF])
		 >> SCALEBITS);

	  
/*	  if(col == 1) {
		  printf("in: %d %d %d\n", inptr[RGB_RED], inptr[RGB_GREEN], inptr[RGB_BLUE]);
		  printf("out1: %d %d %d\n",output_buf[0][0][col],output_buf[1][0][col],output_buf[2][0][col]);
		  printf("out2: %d %d %d\n",outptr0[col],outptr1[col],outptr2[col]);
	  }
*/	  
	  //inptr += RGB_PIXELSIZE;
    }
  }
}
#endif

int main() {
	
	struct timeval st,ed;
	
	char infile[64], outfile[64], csvfile[64];
	FILE *infp, *outfp, *csvfp;
	JDIMENSION col;
	int num_rows;
	JDIMENSION num_cols;
	JSAMPARRAY input_buf;
	JSAMPIMAGE output_buf;
	JSAMPARRAY check_buf;
	INT32 *rgb_ycc_tab;
	rgb_ycc_tab = (INT32 *) malloc(TABLE_SIZE * SIZEOF(INT32));
	rgb_ycc_start(rgb_ycc_tab);

#ifdef ORIGINAL
	sprintf(csvfile,"./result/original.csv");
#endif

#ifdef MP_IN
	sprintf(csvfile,"./result/mp_in.csv");
#endif
	csvfp = fopen(csvfile, "w");
	
	for(int run_count = 0; run_count<13656 ;run_count++) {
		sprintf(infile, "./input/%d.in",run_count);
		sprintf(outfile, "./output/%d.out",run_count);
		infp = fopen(infile, "rb");
		outfp = fopen(outfile, "rb");
		if(infp==NULL || outfp ==NULL) {
			printf("count: %d\n",run_count);
			printf("%s\n",infile);
			printf("%s\n",outfile);
			break;
		}
		
		fread(&num_rows, sizeof(int), 1, infp);
		fread(&num_cols, SIZEOF(JDIMENSION), 1, infp);	
		
		if(num_rows != 1) {
			printf("num_rows: %d\n", num_rows);
			break;
		}
		
		input_buf = (JSAMPARRAY) malloc(SIZEOF(JSAMPROW)*num_cols);
		for(col = 0; col<num_cols ;col++){
			input_buf[col] = (JSAMPROW) malloc(SIZEOF(JSAMPLE)*RGB_PIXELSIZE);
			fread(&input_buf[col][RGB_RED], SIZEOF(JSAMPLE), 1, infp);
		    fread(&input_buf[col][RGB_GREEN], SIZEOF(JSAMPLE), 1, infp);
		    fread(&input_buf[col][RGB_BLUE], SIZEOF(JSAMPLE), 1, infp);
		}
		
		output_buf = (JSAMPIMAGE) malloc(SIZEOF(JSAMPARRAY)*RGB_PIXELSIZE);
		for(int i=0; i<RGB_PIXELSIZE ;i++) {
			output_buf[i] = (JSAMPARRAY) malloc(SIZEOF(JSAMPROW)*1);
			output_buf[i][0] = (JSAMPROW) malloc(SIZEOF(JSAMPLE)*num_cols);
		}
		
		check_buf = (JSAMPARRAY) malloc(SIZEOF(JSAMPROW)*RGB_PIXELSIZE);
		for(int i=0; i<RGB_PIXELSIZE ;i++) {
			check_buf[i] = (JSAMPROW) malloc(SIZEOF(JSAMPLE)*num_cols);
		}
		for(col = 0; col<num_cols ;col++){
			fread(&check_buf[0][col], SIZEOF(JSAMPLE), 1, outfp);
			fread(&check_buf[1][col], SIZEOF(JSAMPLE), 1, outfp);
			fread(&check_buf[2][col], SIZEOF(JSAMPLE), 1, outfp);
		}
		
		for(int i=0; i<10 ;i++) {
			gettimeofday(&st, NULL);
		
			rgb_ycc_convert (rgb_ycc_tab, input_buf, output_buf, 0, num_rows, num_cols);
		
			gettimeofday(&ed, NULL);
			fprintf(csvfp, "%lf,",(double) (ed.tv_usec-st.tv_usec)/1000000 + (double) (ed.tv_sec-st.tv_sec));
		}
		fprintf(csvfp, "\n");
		for(col = 0; col<num_cols ; col++) {
			if(check_buf[0][col] != output_buf[0][0][col]) {
				printf("R col: %d\n",col);
				printf("check: %d %d %d\n",check_buf[0][col],check_buf[1][col],check_buf[2][col]);
				printf("outpu: %d %d %d\n",output_buf[0][0][col],output_buf[1][0][col],output_buf[2][0][col]);
				break;
			}
			if(check_buf[1][col] != output_buf[1][0][col]) {
				printf("G col: %d\n",col);
				printf("check: %d %d %d\n",check_buf[0][col],check_buf[1][col],check_buf[2][col]);
				printf("outpu: %d %d %d\n",output_buf[0][0][col],output_buf[1][0][col],output_buf[2][0][col]);
				break;
			}
			if(check_buf[2][col] != output_buf[2][0][col]) {
				printf("B col: %d\n",col);
				printf("check: %d %d %d\n",check_buf[0][col],check_buf[1][col],check_buf[2][col]);
				printf("outpu: %d %d %d\n",output_buf[0][0][col],output_buf[1][0][col],output_buf[2][0][col]);
				break;
			}
		}
		
		
		for(col = 0; col<num_cols ;col++){
			free(input_buf[col]);
		}
		free(input_buf);
		
		
		for(int i=0; i<RGB_PIXELSIZE ;i++) {
			free(output_buf[i][0]);
			free(output_buf[i]);
			
		}
		free(output_buf);
		
		for(int i=0; i<RGB_PIXELSIZE ;i++) {
			free(check_buf[i]);
		}
		free(check_buf);
		
		fclose(infp);
		fclose(outfp);
		//printf("Done %d\n",run_count);
	}
	
#ifdef ORIGINAL
	fclose(csvfp);
#endif
	
	return 0;
}