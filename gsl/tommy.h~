// For branch "master"

#pragma once
#ifndef TOMMY_H
#define TOMMY_H

// 2012-08-18: Added a branch "master". May merge some time.
#define GO_HOME_IN_A_TANK

#define forceinline __inline__ __attribute__((always_inline))
#define noinline __attribute__((noinline))
#include <execinfo.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

/* Also includes header files for setjmp and longjmp */
#include <string.h>
#include <sys/signal.h>
#include <signal.h>
#include <setjmp.h>

#ifdef __cplusplus
extern "C" {
#endif
extern jmp_buf buf;
extern jmp_buf buf_1;
extern jmp_buf buf_2;
#ifdef __cplusplus
}
#endif

/* To enable/disable verbose output */


#define DEBUG
#ifdef DEBUG    
#define DBG(call) { call; }
#else   
#define DBG(call) {}
#endif


/* Schifra has been eliminated because it was too inefficient */

/* CRC-32b.h has also been eliminated. A much simpler sum is OK. */

/* The SIGSEGV handler */
#define UNWIND
#ifdef __cplusplus
extern "C" {
#endif
void MY_SET_SIGSEGV_HANDLER();
void MY_REMOVE_SIGSEGV_HANDLER();
void MY_SIGSEGV_HANDLER_SUMMARY();
unsigned int GetJmpBufChecksum(jmp_buf*);
void do_backtrace();
void my_action(int sig);
#ifdef __cplusplus
}
#endif

/* The super sigsetjmp! */
// To use this macro, you must declare a variable called 'int jmpret'.
#define SUPERSETJMP(message) \
	jmpret = sigsetjmp(buf, 1);  \
	trick_me_jr(jmpret); \
	if(jmpret == 0) { \
		memcpy(&buf_1, &buf, sizeof(jmp_buf)); \
		memcpy(&buf_2, &buf, sizeof(jmp_buf)); \
		printf("[SUPERSETJMP] %s (sum=%d)\n", message, GetJmpBufChecksum(&buf)); \
	} else { \
		printf("[SUPERSETJMP] returned from %s (sum=%d)\n", message, GetJmpBufChecksum(&buf)); \
	}\

/* To enable libunwind for the handler */

#define TRIO_DOUBLE(a, b, c) ((a==b)?a:c)

/* To enable stopwatching */

#ifdef __cplusplus
extern "C" {
#endif
noinline void my_stopwatch_checkpoint(int id);
noinline void my_stopwatch_stop(int id);
noinline unsigned long my_stopwatch_get(int id);
noinline void my_stopwatch_show(int id, char* extra_msg, char writeToFile);
void STOPWATCH_INCREMENT(int id);

int Is_GSL_Vector_Equal(const gsl_vector* A, const gsl_vector* B);
int Is_GSL_Matrix_Equal(const gsl_matrix* A, const gsl_matrix* B);
noinline double GetIsEqualKnob();
noinline double ResetIsEqualKnob();
noinline double AdjustIsEqualKnob(int pwr10); // Multiply by 10^(pwr10)
double Get_GSL_Matrix_RMSD(const gsl_matrix* A, const gsl_matrix*);
double Get_GSL_Vector_RMSD(const gsl_vector* A, const gsl_vector*);
int my_dgemv(CBLAS_TRANSPOSE_t Trans, double alpha, const gsl_matrix* A, const gsl_vector* X, double beta, gsl_vector* Y);
int my_dgemm(CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, double alpha, const gsl_matrix* A, const gsl_matrix* B, double beta, gsl_matrix* C);

/* These routines are related to checking DSYRK */

/* Check the submatrix of C where C is result of A * t(A). */
int check_aat_submatrix(const gsl_matrix* C, int rowstart, int rowend, int colstart, int colend, double alpha, const gsl_matrix* A, CBLAS_TRANSPOSE_t trans);
/* Check the submatrix of C. The region to be checked is a triangle, whose vertex with the right angle is located at 
 * row vx and column vy. It points to upper-right if uplo is CblasUpper and points down if uplo is CblasLower.
 */
int check_aat_triangle(const gsl_matrix* C, CBLAS_UPLO_t uplo, int row, int col, int len, double alpha, const gsl_matrix* A, CBLAS_TRANSPOSE_t trans);
/* Check if one element is correct */
int check_aat_elem(const gsl_matrix* C, int row, int col, double alpha, const gsl_matrix* A, CBLAS_TRANSPOSE_t trans);

int Is_GSL_MM_Equal(CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB,
		    double alpha, const gsl_matrix* A, const gsl_matrix* B,
		    double beta,  const gsl_matrix* C, const gsl_matrix* C2);
int Is_GSL_DTRSV_Equal(CBLAS_UPLO_t, CBLAS_TRANSPOSE_t, CBLAS_DIAG_t, const gsl_matrix*, const gsl_vector*, const gsl_vector*);
int Is_GSL_DSYRK_Equal(const CBLAS_UPLO_t uplo, CBLAS_TRANSPOSE_t trans,
		       double alpha, const gsl_matrix* A,
		       double beta,  const gsl_matrix* C, const gsl_matrix* C2);
/* By decomposing DSYRK to a series of matrix-matrix multiplication and validate them recursively. 04-20-2012 */
int Is_GSL_DSYRK_Equal2(const CBLAS_UPLO_t uplo, CBLAS_TRANSPOSE_t trans,
		       double alpha, const gsl_matrix* A,
		       double beta,  const gsl_matrix* C, const gsl_matrix* C2);
int Is_GSL_DGEMV_Equal(CBLAS_TRANSPOSE_t Trans, double alpha, const gsl_matrix* A, const gsl_vector* X, double beta, const gsl_vector* Y, const gsl_vector* Y2);
int Is_GSL_linalg_cholesky_decomp_Equal(const gsl_matrix* A, const gsl_matrix* A_out);
void Print_Matrix(CBLAS_TRANSPOSE_t Trans, const gsl_matrix* M, const char* info);
void Print_Vector(const gsl_vector* V, const char* info);
void my_set_SIGSEGV_handler(void);
void my_remove_SIGSEGV_handler(void);
int Print_Hello(void);
unsigned long CRC_Matrix(const gsl_matrix*);
unsigned long CRC_Vector(const gsl_vector*);
 /* Added on Apr 18 */
double my_sum_vector(const gsl_vector*);
double my_sum_matrix(const gsl_matrix*);
#ifdef __cplusplus
}
#endif



/* To configure PoECC's version */
/* Changed on 2012-05-06 */
/* Better to specify a POECC_VERSION in makefile if you wish to use FT3 */
/* If no version is specified then default to version 2 */
#ifdef POECC_VERSION1
	#define MAKE_POEC_DATA_MATRIX(mat) make_PoEC_data_matrix(mat)
	#define MAKE_POEC_DATA_VECTOR(vec) make_PoEC_data_vector(vec)
	#define MY_MAT_CHK_RECOVER_POECC(sum, ecMat, mat) {\
		my_stopwatch_checkpoint(9); \
		double c = my_sum_matrix(mat); \
		if(c!=sum) {\
			my_stopwatch_checkpoint(8); \
			PoEC_restore_matrix(ecMat, mat); \
			my_stopwatch_stop(8); \
		}\
		my_stopwatch_stop(9); \
	}

	#define MY_VEC_CHK_RECOVER_POECC(sum, ecVec, vec) {\
		my_stopwatch_checkpoint(9); \
		double c = my_sum_vector(vec); \
		if(c!=sum) {\
			my_stopwatch_checkpoint(8); \
			PoEC_restore_vector(ecVec, vec); \
			my_stopwatch_stop(8); \
		}\
		my_stopwatch_stop(9); \
	}
#else
	#include "ccekcolbop.h"
 	noinline float MY_MAT_CHK_RECOVER_POECC(double sum, void* ecMat, gsl_matrix* mat); 
	noinline float MY_VEC_CHK_RECOVER_POECC(double sum, void* ecVec, gsl_vector* vec);
#endif
/* To enable stopwatching */
#include <sys/time.h>

#define SW3START my_stopwatch_checkpoint(3)
#define SW3STOP my_stopwatch_stop(3)

/* Change its definition. I'm buffled. 
 * I actually learned a lot about C programming
 * these 2 months and I want to rule out any
 * possible coding errors !!!!!! */

#define EPSILON 0.00001
// These tolerances are inclusive. E.G., if 1 is set, then a[i]/b[i]==1 is treated as equal.

#ifdef DTRSV_EXPERIMENT
#define FT_TOLERANCE 0.01
#else
#define FT_TOLERANCE       0.00001
#endif
/*
#define FT_TOLERANCE_LOWER 0.999999999999
#define FT_TOLERANCE_UPPER 1.000000000001
*/

#define NUM_OF_RERUN 35	// Up to how many re-runs ? If the number is exceeded the program would return the last result
#define NUM_OF_SIGSEGV 500 // Up to how many SIGSEGV's ? If the number is exceeded the program would be aborted.
			   // 04-21-2012 Increased for running LASSO @ FR 1e-07...
			   // 08-18-2012 To pass the tough 4-MM test, we need some more seg-fault tolerance.

struct matrix_ECC {
	int size1;
	int size2;
	double* rowsum;
	double* colsum;
	double allsum;
};

static void my_matrix_fwrite(const char*, gsl_matrix*);
static void my_matrix_fread(const char*, gsl_matrix*);

#define MY_MATRIX_FWRITE(filename, m) {\
	my_stopwatch_checkpoint(7); \
	FILE* f = fopen(filename, "wb");\
	if(f==NULL) { printf("[my_matrix_fwrite]Error opening file %d\n", filename); }\
	gsl_matrix_fwrite(f, m);\
	fclose(f);\
	my_stopwatch_stop(7); \
}

#define MY_MATRIX_FREAD(filename, m) {\
	FILE* f = fopen(filename, "r");\
	if(f==NULL) { printf("[my_matrix_fread]Error opening file %d\n", filename); }\
	gsl_matrix_fread(f, m);\
	fclose(f);\
}

#define MY_MATRIX_CHK_RECOVER_FILE(crc, filename, m) {\
	my_stopwatch_checkpoint(9); \
	unsigned long c = CRC_Matrix(m);\
	if(c!=crc) {\
		my_stopwatch_checkpoint(8); \
		MY_MATRIX_FREAD(filename, m);\
		my_stopwatch_stop(8); \
	}\
	my_stopwatch_stop(9); \
}

#define MY_MAT_CHKSUM_RECOVER_FILE(sum, filename, m) {\
	my_stopwatch_checkpoint(9); \
	double c = my_sum_matrix(m); \
	if(c != sum) {\
		my_stopwatch_checkpoint(8); \
		MY_MATRIX_FREAD(filename, m); \
		my_stopwatch_stop(8); \
	}\
	my_stopwatch_stop(9); \
}

#define MY_VEC_CHKSUM_RECOVER_FILE(sum, filename, v) {\
	my_stopwatch_checkpoint(9); \
	double c = my_sum_vector(v); \
	if(c != sum) {\
		my_stopwatch_checkpoint(8); \
		MY_VECTOR_FREAD(filename, v); \
		my_stopwatch_stop(8); \
	}\
	my_stopwatch_stop(9); \
}

#define MY_VEC_CHK_RECOVER_FILE(crc, filename, v) {\
	my_stopwatch_checkpoint(9); \
	unsigned long c = CRC_Vector(v);\
	if(c!=crc) {\
		my_stopwatch_checkpoint(8); \
		MY_VECTOR_FREAD(filename, v);\
		my_stopwatch_stop(8); \
	}\
	my_stopwatch_stop(9); \
}

#define MY_MAT_CHK_RECOVER_ECC(crc, ecMat, mat) {\
	my_stopwatch_checkpoint(9); \
	unsigned long c = CRC_Matrix((gsl_matrix*)mat);\
	if(c!=crc) {\
		my_stopwatch_checkpoint(8); \
		mat=(gsl_matrix*)EC_restore_matrix((void*)ecMat);\
		my_stopwatch_stop(8); \
	}\
	my_stopwatch_stop(9); \
}

#define MY_VEC_CHK_RECOVER_ECC(crc, ecVec, mat) {\
	my_stopwatch_checkpoint(9); \
	unsigned long c = CRC_Vector((gsl_vector*)mat);\
	if(c!=crc) {\
		my_stopwatch_checkpoint(8); \
		mat=(gsl_vector*)EC_restore_vector((void*)ecVec);\
		my_stopwatch_stop(8); \
	}\
	my_stopwatch_stop(9); \
}

#define MY_VECTOR_FWRITE(filename, v) {\
	my_stopwatch_checkpoint(7); \
	FILE* f = fopen(filename, "wb");\
	if(f==NULL) { printf("[my_vector_fwrite]Error opening file %d\n", filename); }\
	gsl_vector_fwrite(f, v);\
	fclose(f);\
	my_stopwatch_stop(7); \
}

#define MY_VECTOR_FREAD(filename, v) {\
	FILE* f = fopen(filename, "r");\
	if(f==NULL) { printf("[my_vector_fread]Error opening file %d\n", filename); }\
	gsl_vector_fread(f, v);\
	fclose(f);\
}

#define MY_FILENAME(dest, length, basepattern) {\
	int bLen = strlen(basepattern);\
	strcpy(dest, basepattern);\
	int i; for(i=bLen; i<length-1; i++) {\
		char ch = (char)(rand() % 26 + 65);\
		dest[i] = ch;\
	}\
	dest[length-1]=0;\
}

#define MY_REMOVE_FILE(filename) {\
	my_stopwatch_checkpoint(7); \
	printf("[MY_REMOVE_FILE]Removing file %s\n", filename);\
	if(remove(filename)!=0) printf("Error removing file %s\n", filename);\
	my_stopwatch_stop(7); \
}

// This is the FT'ed MM routine, using files for temporary restoration
#define GSL_BLAS_DGEMM_FT(TransA, TransB, alpha, matA, matB, beta, matC, kk)		\
{\
	unsigned long crcA, crcB, crcC;\
	crcA=CRC_Matrix(matA); crcB=CRC_Matrix(matB); crcC=CRC_Matrix(matC); \
	nonEqualCount=0;\
	MY_SET_SIGSEGV_HANDLER();\
	gsl_matrix* matC_bak = gsl_matrix_alloc(matC->size1, matC->size2);	\
	char fnA[20], fnB[20], fnC[20]; \
	srand(time(NULL)*1000+getpid());\
	MY_FILENAME(fnA, 20, "tmp_matA");\
	MY_FILENAME(fnB, 20, "tmp_matB");\
	MY_FILENAME(fnC, 20, "tmp_matC");\
	MY_MATRIX_FWRITE(fnA, matA); \
	MY_MATRIX_FWRITE(fnB, matB); \
	MY_MATRIX_FWRITE(fnC, matC); \
	int jmpret=0, isEqual;		\
	jmpret = sigsetjmp(buf, 1);\
	kk: \
	printf("[DGEMM_FT]Recovering matrices from file\n");\
	MY_MATRIX_CHK_RECOVER_FILE(crcA, fnA, matA); \
	MY_MATRIX_CHK_RECOVER_FILE(crcB, fnB, matB); \
	MY_MATRIX_CHK_RECOVER_FILE(crcC, fnC, matC); \
	gsl_matrix_memcpy(matC_bak, matC);	\
	printf("[DGEMM_FT]Normal call to dgemm.. nonEqualCount=%d\n", nonEqualCount);	\
	SW3START; \
	int ret = gsl_blas_dgemm(TransA, TransB, alpha, matA, matB, beta, matC);	\
	SW3STOP; \
	isEqual = Is_GSL_MM_Equal(TransA, TransB, alpha, matA, matB, beta, matC_bak, matC);	\
	if(isEqual==1) printf("[DGEMM_FT]Result: Equal\n");	\
	else { \
		printf("[DGEMM_FT]Result: NOT Equal\n");	\
		my_stopwatch_checkpoint(6); \
		nonEqualCount = nonEqualCount + 1; \
		if(nonEqualCount < NUM_OF_RERUN) { \
			printf("[DGEMM_FT]Restart calculation.\n");\
			goto kk; \
		}\
	}\
	my_stopwatch_stop(6); \
	printf("[DGEMM_FT]Releasing memory for mat_bak's and deleting temp files\n");\
	MY_REMOVE_FILE(fnA);\
	MY_REMOVE_FILE(fnB);\
	MY_REMOVE_FILE(fnC);\
	gsl_matrix_free(matC_bak);	\
	printf("Released.\n");\
	MY_REMOVE_SIGSEGV_HANDLER();\
}

/* 04-03-2012: Difference between FT2 and FT: using Reed-Solomon Error Correcting Code instead of file backup */
#define GSL_BLAS_DGEMM_FT2(TransA, TransB, alpha, matA, matB, beta, matC, kk)		\
{\
	unsigned long crcA, crcB, crcC;\
	crcA=CRC_Matrix(matA); crcB=CRC_Matrix(matB); crcC=CRC_Matrix(matC);\
	void* ecMatA, *ecMatB, *ecMatC;\
	nonEqualCount=0;\
	MY_SET_SIGSEGV_HANDLER();\
	ecMatA = make_EC_data_matrix(matA);\
	ecMatB = make_EC_data_matrix(matB);\
	ecMatC = make_EC_data_matrix(matC);\
	int jmpret=0, isEqual;		\
	jmpret = sigsetjmp(buf, 1);\
	kk: \
	printf("[DGEMM_FT2]Recovering matrices from backup error correcting code\n");\
	MY_MAT_CHK_RECOVER_ECC(crcA, ecMatA, matA);\
	MY_MAT_CHK_RECOVER_ECC(crcB, ecMatB, matB);\
	MY_MAT_CHK_RECOVER_ECC(crcC, ecMatC, matC);\
	gsl_matrix* matC_bak = gsl_matrix_alloc(matC->size1, matC->size2);	\
	gsl_matrix_memcpy(matC_bak, matC);	\
	printf("[DGEMM_FT2]Normal call to dgemm.. nonEqualCount=%d\n", nonEqualCount);	\
	SW3START; \
	int ret = gsl_blas_dgemm(TransA, TransB, alpha, matA, matB, beta, matC);	\
	SW3STOP; \
	printf("[DGEMM_FT2] Check Results..\n");	\
	isEqual = Is_GSL_MM_Equal(TransA, TransB, alpha, matA, matB, beta, matC_bak, matC);	\
	if(isEqual==1) printf("[DGEMM_FT]Result: Equal\n");	\
	else { \
		my_stopwatch_checkpoint(6); \
		printf("[DGEMM_FT]Result: NOT Equal\n");	\
		nonEqualCount = nonEqualCount + 1; \
		if(nonEqualCount < NUM_OF_RERUN) { \
			printf("[DGEMM_FT]Restart calculation.\n");\
			goto kk; \
		}\
	}\
	my_stopwatch_stop(6); \
	printf("[DGEMM_FT]Releasing memory for mat_bak's and deleting temp files\n");\
	gsl_matrix_free(matC_bak);	\
	printf("Released.\n");\
	MY_REMOVE_SIGSEGV_HANDLER();\
}

noinline
void GSL_BLAS_DGEMM_FT3(CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, 
			double alpha, const gsl_matrix* matA, const gsl_matrix* matB, 
			double beta, gsl_matrix* matC);

#ifdef CRC32
#define GSL_BLAS_DGEMV_FT(Trans, alpha, matA, vecX, beta, vecY, kk)	\
{\
	unsigned long crcA, crcX, crcY;\
	crcA=CRC_Matrix(matA); crcX=CRC_Vector(vecX); crcY=CRC_Vector(vecY);\
	nonEqualCount=0;\
	MY_SET_SIGSEGV_HANDLER();\
	gsl_vector* vecY_bak = gsl_vector_alloc(vecY->size);\
	char fnX[20], fnY[20], fnA[20];\
	srand(time(NULL)*1000 + getpid());\
	MY_FILENAME(fnX, 20, "tmp_vecX");\
	MY_FILENAME(fnY, 20, "tmp_vecY");\
	MY_FILENAME(fnA, 20, "tmp_matA");\
	MY_VECTOR_FWRITE(fnX, vecX);\
	MY_VECTOR_FWRITE(fnY, vecY);\
	MY_MATRIX_FWRITE(fnA, matA);\
	int jmpret=0, isEqual;\
	jmpret = sigsetjmp(buf, 1);\
	kk: \
	printf("[DGEMV_FT]Recovering matrix and vectors from file\n");\
	MY_VEC_CHK_RECOVER_FILE(crcX, fnX, vecX);\
	MY_VEC_CHK_RECOVER_FILE(crcY, fnY, vecY);\
	MY_MATRIX_CHK_RECOVER_FILE(crcA, fnA, matA);\
	gsl_vector_memcpy(vecY_bak, vecY);\
	printf("[DGEMV_FT]Normal call to dgemv.. nonEqualCount=%d\n", nonEqualCount);\
	SW3START; \
	int ret = gsl_blas_dgemv(Trans, alpha, matA, vecX, beta, vecY);\
	SW3STOP; \
	isEqual = Is_GSL_DGEMV_Equal(Trans, alpha, matA, vecX, beta, vecY_bak, vecY);\
	if(isEqual==1) printf("[DGEMV_FT]Result: Equal\n");\
	else {\
		my_stopwatch_checkpoint(6); \
		printf("[DGEMV_FT]Result: NOT Equal\n");	\
		nonEqualCount = nonEqualCount + 1;\
		if(nonEqualCount < NUM_OF_RERUN) {\
			printf("[DGEMV_FT]Restart calculation.\n");\
			goto kk; \
		}\
	}\
	my_stopwatch_stop(6); \
	printf("[DGEMV_FT]Releasing memory for vec_bak's\n");\
	MY_REMOVE_FILE(fnX);\
	MY_REMOVE_FILE(fnY);\
	MY_REMOVE_FILE(fnA);\
	gsl_vector_free(vecY_bak);\
	printf("Released.\n");\
	MY_REMOVE_SIGSEGV_HANDLER();\
}
#else
#define GSL_BLAS_DGEMV_FT(Trans, alpha, matA, vecX, beta, vecY, kk)	\
{\
	double sumA, sumX, sumY; \
	sumA=my_sum_matrix(matA); sumX=my_sum_vector(vecX); sumY=my_sum_vector(vecY); \
	nonEqualCount=0;\
	MY_SET_SIGSEGV_HANDLER();\
	gsl_vector* vecY_bak = gsl_vector_alloc(vecY->size);\
	char fnX[20], fnY[20], fnA[20];\
	srand(time(NULL)*1000 + getpid());\
	MY_FILENAME(fnX, 20, "tmp_vecX");\
	MY_FILENAME(fnY, 20, "tmp_vecY");\
	MY_FILENAME(fnA, 20, "tmp_matA");\
	MY_VECTOR_FWRITE(fnX, vecX);\
	MY_VECTOR_FWRITE(fnY, vecY);\
	MY_MATRIX_FWRITE(fnA, matA);\
	int jmpret=0, isEqual;\
	jmpret = sigsetjmp(buf, 1);\
	kk: \
	printf("[DGEMV_FT]Recovering matrix and vectors from file\n");\
	MY_MAT_CHKSUM_RECOVER_FILE(sumA, fnA, matA); \
	MY_VEC_CHKSUM_RECOVER_FILE(sumX, fnX, vecX); \
	MY_VEC_CHKSUM_RECOVER_FILE(sumY, fnY, vecY); \
	gsl_vector_memcpy(vecY_bak, vecY);\
	printf("[DGEMV_FT]Normal call to dgemv.. nonEqualCount=%d\n", nonEqualCount);\
	SW3START; \
	int ret = gsl_blas_dgemv(Trans, alpha, matA, vecX, beta, vecY);\
	SW3STOP; \
	isEqual = Is_GSL_DGEMV_Equal(Trans, alpha, matA, vecX, beta, vecY_bak, vecY);\
	if(isEqual==1) printf("[DGEMV_FT]Result: Equal\n");\
	else {\
		my_stopwatch_checkpoint(6); \
		printf("[DGEMV_FT]Result: NOT Equal\n");	\
		nonEqualCount = nonEqualCount + 1;\
		if(nonEqualCount < NUM_OF_RERUN) {\
			printf("[DGEMV_FT]Restart calculation.\n");\
			goto kk; \
		}\
	}\
	my_stopwatch_stop(6); \
	printf("[DGEMV_FT]Releasing memory for vec_bak's\n");\
	MY_REMOVE_FILE(fnX);\
	MY_REMOVE_FILE(fnY);\
	MY_REMOVE_FILE(fnA);\
	gsl_vector_free(vecY_bak);\
	printf("Released.\n");\
	MY_REMOVE_SIGSEGV_HANDLER();\
}
#endif

#define GSL_BLAS_DGEMV_FT2(Trans, alpha, matA, vecX, beta, vecY, kk)\
{\
	unsigned long crcA, crcX, crcY;\
	crcA=CRC_Matrix(matA); crcX=CRC_Vector(vecX); crcY=CRC_Vector(vecY);\
	void* ecMatA, *ecVecX, *ecVecY;\
	nonEqualCount=0;\
	MY_SET_SIGSEGV_HANDLER();\
	ecMatA = make_EC_data_matrix(matA); \
	ecVecX = make_EC_data_vector(vecX); \
	ecVecY = make_EC_data_vector(vecY); \
	int jmpret=0, isEqual;\
	jmpret = sigsetjmp(buf, 1);\
	kk: \
	printf("[DGEMV_FT2]Recovering matrix and vectors from file\n");\
	MY_MAT_CHK_RECOVER_ECC(crcA, ecMatA, matA);\
	MY_VEC_CHK_RECOVER_ECC(crcX, ecVecX, vecX);\
	MY_VEC_CHK_RECOVER_ECC(crcY, ecVecY, vecY);\
	gsl_vector* vecY_bak = gsl_vector_alloc(vecY->size);\
	gsl_vector_memcpy(vecY_bak, vecY);\
	printf("[DGEMV_FT2]Normal call to dgemv.. nonEqualCount=%d\n", nonEqualCount);\
	SW3START; \
	int ret = gsl_blas_dgemv(Trans, alpha, matA, vecX, beta, vecY);\
	SW3STOP; \
	isEqual = Is_GSL_DGEMV_Equal(Trans, alpha, matA, vecX, beta, vecY_bak, vecY);\
	if(isEqual==1) printf("[DGEMV_FT2]Result: Equal\n");\
	else {\
		my_stopwatch_checkpoint(6); \
		printf("[DGEMV_FT2]Result: NOT Equal\n");\
		nonEqualCount = nonEqualCount + 1;\
		if(nonEqualCount < NUM_OF_RERUN) {\
			printf("[DGEMV_FT2]Restart calculation.\n");\
			goto kk; \
		}\
	}\
	my_stopwatch_stop(6); \
	printf("[DGEMV_FT2]Releasing memory for vec_bak's\n");\
	gsl_vector_free(vecY_bak);\
	printf("Released.\n");\
	MY_REMOVE_SIGSEGV_HANDLER();\
}

noinline
void GSL_BLAS_DGEMV_FT3(CBLAS_TRANSPOSE_t Trans, double alpha, 
			const gsl_matrix* matA, const gsl_vector* vecX, 
			double beta, gsl_vector* vecY);

#define GSL_BLAS_DSYRK_FT(uplo, trans, alpha, A, beta, C, kk) \
{\
	unsigned long crcA, crcC;\
	crcA=CRC_Matrix(A); crcC=CRC_Matrix(C);\
	nonEqualCount=0;\
	MY_SET_SIGSEGV_HANDLER();\
	gsl_matrix* C_bak = gsl_matrix_alloc(C->size1, C->size2);\
	char fnA[20], fnC[20];\
	srand(time(NULL)*1000 + getpid());\
	MY_FILENAME(fnA, 20, "tmp_matA");\
	MY_FILENAME(fnC, 20, "tmp_matC");\
	MY_MATRIX_FWRITE(fnA, A);\
	MY_MATRIX_FWRITE(fnC, C);\
	int jmpret=0, isEqual;\
	jmpret = sigsetjmp(buf, 1);\
	kk: \
	printf("[DSYRK_FT]Recovering matrices from file\n");\
	MY_MATRIX_CHK_RECOVER_FILE(crcA, fnA, A);\
	MY_MATRIX_CHK_RECOVER_FILE(crcC, fnC, C);\
	gsl_matrix_memcpy(C_bak, C);\
	printf("[DSYRK_FT]Normal call to dsyrk.. nonEqualCount=%d\n", nonEqualCount);\
	SW3START; \
	int ret = gsl_blas_dsyrk(uplo, trans, alpha, A, beta, C);\
	SW3STOP; \
	isEqual = Is_GSL_DSYRK_Equal2(uplo, trans, alpha, A, beta, C_bak, C);\
	if(isEqual==1) printf("[DSYRK_FT]Result: Equal\n");\
	else {\
		my_stopwatch_checkpoint(6); \
		printf("[DSYRK_FT]Result: NOT Equal\n");\
		nonEqualCount = nonEqualCount + 1;\
		if(nonEqualCount < NUM_OF_RERUN) {\
			printf("[DSYRK_FT]Restart calculation.\n");\
			goto kk; \
		}\
	}\
	my_stopwatch_stop(6); \
	printf("[DSYRK_FT]Releasing memory for C_bak\n");\
	MY_REMOVE_FILE(fnA);\
	MY_REMOVE_FILE(fnC);\
	gsl_matrix_free(C_bak);\
	printf("Released.\n");\
	MY_REMOVE_SIGSEGV_HANDLER();\
}

#define GSL_BLAS_DSYRK_FT2(uplo, trans, alpha, A, beta, C, kk) \
{\
	unsigned long crcA, crcC;\
	crcA=CRC_Matrix(A);\
	crcC=CRC_Matrix(C);\
	void *ecMatA, *ecMatC;\
	nonEqualCount=0;\
	MY_SET_SIGSEGV_HANDLER();\
	ecMatA = (void*)make_EC_data_matrix(A);\
	ecMatC = make_EC_data_matrix(C);\
	int jmpret=0, isEqual;\
	jmpret = sigsetjmp(buf, 1);\
	kk: \
	printf("[DSYRK_FT2]Recovering matrices from error correction data\n");\
	MY_MAT_CHK_RECOVER_ECC(crcA, ecMatA, A);\
	MY_MAT_CHK_RECOVER_ECC(crcC, ecMatC, C);\
	gsl_matrix* C_bak = gsl_matrix_alloc(C->size1, C->size2);\
	gsl_matrix_memcpy(C_bak, C);\
	printf("[DSYRK_FT2]Normal call to dsyrk.. nonEqualCount=%d\n", nonEqualCount);\
	SW3START; \
	int ret = gsl_blas_dsyrk(uplo, trans, alpha, A, beta, C);\
	SW3STOP; \
	isEqual = Is_GSL_DSYRK_Equal2(uplo, trans, alpha, A, beta, C_bak, C);\
	if(isEqual==1) printf("[DSYRK_FT2]Result: Equal\n");\
	else {\
		my_stopwatch_checkpoint(6); \
		printf("[DSYRK_FT2]Result: NOT Equal\n");\
		nonEqualCount = nonEqualCount + 1;\
		if(nonEqualCount < NUM_OF_RERUN) {\
			printf("[DSYRK_FT2]Restart calculation.\n");\
			goto kk; \
		}\
	}\
	my_stopwatch_stop(6); \
	printf("[DSYRK_FT2]Releasing memory for C_bak\n");\
	gsl_matrix_free(C_bak);\
	printf("Released.\n");\
	MY_REMOVE_SIGSEGV_HANDLER();\
}

noinline
void GSL_BLAS_DSYRK_FT3(CBLAS_UPLO_t uplo, CBLAS_TRANSPOSE_t trans, 
	           double alpha, const gsl_matrix* A,
		   double beta, gsl_matrix* C);

#define GSL_BLAS_DTRSV_FT(uplo, TransA, Diag, A, X, kk)\
{\
	unsigned long crcA, crcX;\
	crcA=CRC_Matrix(A); crcX=CRC_Vector(X);\
	nonEqualCount=0;\
	MY_SET_SIGSEGV_HANDLER();\
	gsl_vector* X_bak = gsl_vector_alloc(X->size);\
	char fnA[20], fnX[20];\
	srand(time(NULL)*1000 + getpid());\
	MY_FILENAME(fnA, 20, "tmp_matA");\
	MY_FILENAME(fnX, 20, "tmp_vecX");\
	MY_MATRIX_FWRITE(fnA, A);\
	MY_VECTOR_FWRITE(fnX, X);\
	int jmpret=0, isEqual;\
	jmpret = sigsetjmp(buf, 1);\
	kk: \
	printf("[DTRSV_FT]Recovering matrix and vector from file\n");\
	MY_MATRIX_CHK_RECOVER_FILE(crcA, fnA, A);\
	MY_VEC_CHK_RECOVER_FILE(crcX, fnX, X);\
	gsl_vector_memcpy(X_bak, X);\
	printf("[DTRSV_FT]Normal call to dsyrk.. nonEqualCount=%d\n", nonEqualCount);\
	SW3START; \
	int ret = gsl_blas_dtrsv(uplo, TransA, Diag, A, X);\
	SW3STOP; \
	isEqual = Is_GSL_DTRSV_Equal(uplo, TransA, Diag, A, X_bak, X);\
	if(isEqual==1) printf("[DTRSV_FT]Result: Equal\n");\
	else {\
		my_stopwatch_checkpoint(6); \
		printf("[DTRSV_FT]Result: NOT Equal\n");\
		nonEqualCount = nonEqualCount + 1;\
		if(nonEqualCount < NUM_OF_RERUN) {\
			printf("[DTRSV_FT]Restart calculation.\n");\
			goto kk;\
		}\
	}\
	my_stopwatch_stop(6); \
	printf("[DTRSV_FT]Releasing memory for X_bak\n");\
	MY_REMOVE_FILE(fnA);\
	MY_REMOVE_FILE(fnX);\
	gsl_vector_free(X_bak);\
	printf("Released.\n");\
	MY_REMOVE_SIGSEGV_HANDLER();\
}

/* FT2 uses reed-solomon error correction code */
#define GSL_BLAS_DTRSV_FT2(uplo, TransA, Diag, A, X, kk)\
{\
	unsigned long crcA, crcX;\
	crcA=CRC_Matrix(A); crcX=CRC_Vector(X);\
	void *ecMatA, *ecVecX;\
	nonEqualCount=0;\
	MY_SET_SIGSEGV_HANDLER();\
	ecMatA = make_EC_data_matrix(A);\
	ecVecX = make_EC_data_vector(X);\
	int jmpret=0, isEqual;\
	jmpret = sigsetjmp(buf, 1);\
	kk: \
	printf("[DTRSV_FT]Recovering matrix and vector from file\n");\
	MY_MAT_CHK_RECOVER_ECC(crcA, ecMatA, A);\
	MY_VEC_CHK_RECOVER_ECC(crcX, ecVecX, X);\
	gsl_vector* X_bak = gsl_vector_alloc(X->size);\
	gsl_vector_memcpy(X_bak, X);\
	printf("[DTRSV_FT]Normal call to dsyrk.. nonEqualCount=%d\n", nonEqualCount);\
	SW3START; \
	int ret = gsl_blas_dtrsv(uplo, TransA, Diag, A, X);\
	SW3STOP; \
	isEqual = Is_GSL_DTRSV_Equal(uplo, TransA, Diag, A, X_bak, X);\
	if(isEqual==1) printf("[DTRSV_FT]Result: Equal\n");\
	else {\
		my_stopwatch_checkpoint(6); \
		printf("[DTRSV_FT]Result: NOT Equal\n");\
		nonEqualCount = nonEqualCount + 1;\
		if(nonEqualCount < NUM_OF_RERUN) {\
			printf("[DTRSV_FT]Restart calculation.\n");\
			goto kk; \
		}\
	}\
	my_stopwatch_stop(6); \
	printf("[DTRSV_FT]Releasing memory for X_bak\n");\
	gsl_vector_free(X_bak);\
	printf("Released.\n");\
	MY_REMOVE_SIGSEGV_HANDLER();\
}

/* FT3 uses po error correction code and sum-of-all checksum */
noinline
void GSL_BLAS_DTRSV_FT3(CBLAS_UPLO_t uplo, CBLAS_TRANSPOSE_t TransA,
			CBLAS_DIAG_t Diag, const gsl_matrix* A, 
			gsl_vector* X);

/* You have to supply an unique label to each call */
noinline 
void GSL_LINALG_CHOLESKY_DECOMP_FT3(gsl_matrix* A);

#endif
