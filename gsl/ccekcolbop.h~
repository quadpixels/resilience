#include <stdio.h>
#include <stdlib.h>

#ifndef _CCEKCOLBOP_H
#define _CCEKCOLBOP_H

#define noinline __attribute__((noinline))
#ifdef __cplusplus
extern "C" {
#endif
/* Re-writing the Po block ECC in C language. */
/* 2012-05-06 */

#define BLK_LEN 3
#define BLKSIZE (BLK_LEN*BLK_LEN)

/* How many times do we correct the ECC code itself? */
#define ECCECC 1

/* How many doubles are reserved for storing the size */
// offset to colsums and rowsums
#define TDOFFSET 1 

/* TESTING CONFIGURATIONS. */
//#define DEBUG
/* According to the configuration of May 6, tampering the ECC would make ECC's recovery ability down by 50%,
 * when ECC ECC is 1. */
//#define TAMPER_WITH_ECC_ITSELF

#define print_array(array, description, size) { \
        int printlimit = size; \
	if(printlimit > 100) printlimit=100;\
        printf("[%s]\n", description); \
        int i; for(i=0; i<printlimit; i++) { \
		if(i%7==0) printf("%04d | ", i); \
		printf("%10f ", array[i]); \
		if(((i+1)%7)==0) printf(" | %04d\n", i); \
        } \
        if(printlimit < size) printf(" ... "); \
	printf("\n"); \
}

/* Those 2 routines are slightly different from each other..... Could their tolerance also be different ?! */
noinline void do_encode(const double*, const int, const int, double*, const int);
noinline void do_encode_2(const double*, const int, const int, double*, const int); /* A slightly different implementation */
noinline unsigned int encode(const double*, const int, void**); // Should return size of alloc'd ary
noinline unsigned int decode(double* patient, const int lenPatient, const double* doc);
noinline void POECC_SUMMARY();

// And we have single precision version
noinline void do_encode_float(const float*, const int, const int, float*, const int);
noinline void do_encode_2_float(const float*, const int, const int, float*, const int); /* A slightly different implementation */
noinline void encode_float(const float*, const int, void**);
noinline void decode_float(float* patient, const int lenPatient, const float* doc);
noinline void POECC_FLOAT_SUMMARY();

#ifdef __cplusplus
}
#endif
#endif
