// BRANCH "Master" --------> This code lets you go home in a TANK!
// 08-12: Re-write some of the FT routines. Hope it will be fault-tolerant.
#include "tommy.h"
#include "triplicate.h"
#include <sys/signal.h>
#include <signal.h>
#include "real.h"

#define FT_ROUTINES_VULN
#ifdef FT_ROUTINES_VULN
	#ifndef FT3TEST
	#define FTV_REAL_TRY(label) REAL_TRY(label)
	#define FTV_REAL_CATCH(label) REAL_CATCH(label)
	#define FTV_REAL_END(label) REAL_END(label)
	#else
		#define FTV_REAL_TRY(label) ;
		#define FTV_REAL_CATCH(label) ;
		#define FTV_REAL_END(label) ;

	#endif
#else
	#define FTV_REAL_TRY(label) ;
	#define FTV_REAL_CATCH(label) ;
	#define FTV_REAL_END(label) ;
#endif

// Added on Sep 06
#define FT_CHKR
#ifdef FT_CHKR
	#define PROTECT_IDX_I  \
	((i!=i1-1) && (i1-1==i2-2) && (i=i1-1)), \
	((i1-1!=i2-2) && (i2-2==i) && (i1=i1+1)), \
	((i2-2!=i) && (i==i1-1) && (i2=i+2))
	#define PROTECT_IDX_J  \
	((j!=j1-1) && (j1-1==j2-2) && (j=j1-1)), \
	((j1-1!=j2-2) && (j2-2==j) && (j1=j+1)), \
	((j2-2!=j) && (j==j1-1) && (j2=j+2))
#endif

#define noinline __attribute__((noinline))

int nonEqualCount = 0;
// int num_total_retries = 0; // Moved to tommy_handler.c
extern int num_total_retries;
double is_equal_knob = FT_TOLERANCE;
void check_nan(double x, char* k) { if(/*isnan(x)*/x!=x) { printf("%s is nan\n", k); }}

noinline
double my_sum_vector_actual(const gsl_vector* v) {
	double ret = 0;
	DBG(printf("[my_sum_vector] v=%lx, v->data=%lx, size=%ld, stride=%ld\n", 
		(unsigned long)v, v->data, v->size, v->stride));
	FTV_REAL_TRY(0) {
		#ifdef FT_CHKRxx
		int i, i1, i2; for(i=0; i<v->size; i++, i1++, i2++,
			PROTECT_IDX_I
			) {
		#else
		int i; for(i=0; i<v->size; i++) {
		#endif
			DBG(printf("i=%d ", i));
			// The same as my_sum_matrix_actual; the GSL's routine may not
			// play well here....
//			ret = ret + gsl_vector_get(v, i);
			ret = ret + v->data[v->stride * i];
		}
	} FTV_REAL_CATCH(0) {} FTV_REAL_END(0); 
	return ret;
}

noinline
double my_sum_vector(const gsl_vector* v) {
	MY_SET_SIGSEGV_HANDLER();
	int jmpret;
	SUPERSETJMP("my_sum_vector");
	return my_sum_vector_actual(v);
}

noinline
double my_sum_matrix_actual(const gsl_matrix* m) {
	DBG(printf("[my_sum_matrix_actual] m=%lx(%ldx%ld), m->tda=%ld, m->data=%x\n", (long)m, m->size1, m->size2, m->tda, (void*)(m->data)));
	double ret = 0;
	FTV_REAL_TRY(0) {
		int i,j; for(i=0; i<m->size1; i++) {
			for(j=0; j<m->size2; j++) {
//				ret = ret + gsl_matrix_get(m, i, j);
// Since we have our own fault checker, must stop GSL's checker from being invoked.
				ret = ret + m->data[i * m->tda + j];
			}
		}
	} FTV_REAL_CATCH(0) {} FTV_REAL_END(0);
	return ret;
}

noinline 
double my_sum_matrix(const gsl_matrix* m) {
	MY_SET_SIGSEGV_HANDLER();
	int jmpret;
	SUPERSETJMP("my_sum_matrix");
	DBG(printf("[my_sum_matrix] m=%lx", (long)m));
	DBG(printf(", m->tda=%ld, m->data=%lx\n", m->tda, (void*)(m->data)));
	return my_sum_matrix_actual(m);
}

noinline
int my_dgemm_actual(CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, double alpha, const gsl_matrix* A, const gsl_matrix* B, double beta, gsl_matrix* C) {
		int m, n, p; // A is m*n, B is n*p, C is m*p
		int i, j, k;
		int szA1, szA2, szB1, szB2, szC1, szC2;
		if(TransA == CblasNoTrans) { szA1=A->size1; szA2=A->size2; }
		else if(TransA == CblasTrans) { szA1=A->size2; szA2=A->size1; }
		if(TransB == CblasNoTrans) { szB1=B->size1; szB2=B->size2; }
		else if(TransB == CblasTrans) { szB1=B->size2; szB2=B->size1; }
	FTV_REAL_TRY(0) {
		szC1 = C->size1; szC2 = C->size2;
	} FTV_REAL_CATCH(0) {} FTV_REAL_END(0);
		if(szA2!=szB1) { DBG(printf("[my_dgemm] Error: A's column count != B's row count.\n")); return -1; }
		if(szA1!=szC1) { DBG(printf("[my_dgemm] Error: A's row count != C's row count.\n"));    return -1; }
		if(szB2!=szC2) { DBG(printf("[my_dgemm] Error: B's column count != C's column count.\n")); return -1; }
		m=szA1; n=szA2; p=szB2;
		double temp, elemA, elemB;
		for(i=0; i<m; i++) {
			for(j=0; j<p; j++) {
				temp=gsl_matrix_get(C, i, j) * beta;
				for(k=0; k<n; k++) {
					if(TransA == CblasNoTrans) elemA=gsl_matrix_get(A, i, k);
					else if(TransA==CblasTrans)elemA=gsl_matrix_get(A, k, i);
					if(TransB == CblasNoTrans) elemB=gsl_matrix_get(B, k, j);
					else if(TransB==CblasTrans)elemB=gsl_matrix_get(B, j, k);
	FTV_REAL_TRY(1) {
					temp=temp + elemA*elemB*alpha;
	} FTV_REAL_CATCH(1) {} FTV_REAL_END(1);
				}
				gsl_matrix_set(C, i, j, temp);
			}
		}
		return 0;
}

noinline int my_dgemm(CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, double alpha, const gsl_matrix* A, const gsl_matrix* B, double beta, gsl_matrix* C) {
	MY_SET_SIGSEGV_HANDLER();
	int jmpret = sigsetjmp(buf, 1);
	trick_me_jr(jmpret);
	return my_dgemm_actual(TransA, TransB, alpha, A, B, beta, C);
}

noinline 
int my_dgemv_actual(CBLAS_TRANSPOSE_t Trans, double alpha, const gsl_matrix* A, const gsl_vector* X, double beta, gsl_vector* Y) {
		if(Trans==CblasNoTrans) {
			DBG(printf("A=(%lu x %lu), X=(%lu x 1), Y=(%lu x 1)\n", A->size1, A->size2, X->size, Y->size));
			if(A->size2 != X->size) { DBG(printf("my_dgemv: A's width is not equal to X's length.\n")); return -1; }
			if(A->size1 != Y->size) { DBG(printf("my_dgemv: A's height is not equal to Y's length.\n")); return -1; }
		} else if(Trans==CblasTrans) {
			DBG(printf("A=(%lu x %lu), X=(%lu x 1), Y=(%lu x 1)\n", A->size1, A->size2, X->size, Y->size));
			if(A->size1 != X->size) { DBG(printf("my_dgemv: A_transposed's width is not equal to X's length.\n")); return -1; }
			if(A->size2 != Y->size) { DBG(printf("my_dgemv; A_transposed's height is not equal to Y's length.\n")); return -1; }
		}
		int i,j;
		double tmp=0.0, elemA, elemX, elemY;

		if(Trans==CblasNoTrans) {
			for(i=0; i<Y->size; i++) {
				tmp = beta * gsl_vector_get(Y, i);
				for(j=0; j<X->size; j++) {
					tmp += gsl_vector_get(X, j) * gsl_matrix_get(A, i, j) * alpha; }
				gsl_vector_set(Y, i, tmp);
			}
		} else { // Trans == CblasTrans
			for(i=0; i<Y->size; i++) {
				tmp = beta * gsl_vector_get(Y, i);
				for(j=0; j<X->size; j++) 
					tmp += gsl_vector_get(X, j) * gsl_matrix_get(A, j, i) * alpha;
				gsl_vector_set(Y, i, tmp);
			}
		}
	return 1;
}

noinline
int my_dgemv(CBLAS_TRANSPOSE_t Trans, double alpha, const gsl_matrix* A, const gsl_vector* X, double beta, gsl_vector* Y) {
	MY_SET_SIGSEGV_HANDLER();
	int jmpret = sigsetjmp(buf, 1);
	trick_me_jr(jmpret);
	return my_dgemv_actual(Trans, alpha, A, X, beta, Y);
}

noinline 
int my_dgemv_upperlower_actual(CBLAS_UPLO_t uplo, CBLAS_TRANSPOSE_t Trans, double alpha, const gsl_matrix* A, const gsl_vector* X, double beta, gsl_vector* Y) {
	if(A->size1 != A->size2) { DBG(printf("my_dgemv_upperlower: A is not square!\n")); return -1; }
	if(Trans==CblasNoTrans) {
		if(A->size2 != X->size) { DBG(printf("my_dgemv: A's width is not equal to X's length.\n")); return -1; }
		if(A->size1 != Y->size) { DBG(printf("my_dgemv: A's height is not equal to Y's length.\n")); return -1; } } else if(Trans==CblasTrans) {
		if(A->size1 != X->size) { DBG(printf("my_dgemv: A_transposed's width is not equal to X's length.\n")); return -1; }
		if(A->size2 != Y->size) { DBG(printf("my_dgemv; A_transposed's height is not equal to Y's length.\n")); return -1; }
	}
	int i,j;
	double tmp=0, tmp1;
	for(i=0; i<Y->size; i++) {
		tmp = beta * gsl_vector_get(Y, i);
		int start=0, end=X->size-1;
		if(uplo==CblasUpper) start=i;
		else if(uplo==CblasLower) end=i;
		for(j=0; j<end; j++) {
			if(Trans==CblasNoTrans) tmp1 = alpha * gsl_vector_get(X, j) * gsl_matrix_get(A, i, j);
			else if(Trans==CblasTrans) tmp1 = alpha * gsl_vector_get(X, j) * gsl_matrix_get(A, j, i);
FTV_REAL_TRY(0) {
			tmp = tmp + tmp1;
} FTV_REAL_CATCH(0) {} FTV_REAL_END(0);
		}
		gsl_vector_set(Y, i, tmp);
	}
	return 1;
}

noinline
int my_dgemv_upperlower(CBLAS_UPLO_t uplo, CBLAS_TRANSPOSE_t Trans, double alpha, const gsl_matrix* A, const gsl_vector* X, double beta, gsl_vector* Y) {
	MY_SET_SIGSEGV_HANDLER();
	int jmpret = sigsetjmp(buf, 1);
	trick_me_jr(jmpret);
	return my_dgemv_upperlower_actual(uplo, Trans, alpha, A, X, beta, Y);
}

void my_remove_SIGSEGV_handler() {
	DBG(printf(" >> Unsetted SIGSEGV handler.\n"));
	struct sigaction sa;
	memset(&sa, 0, sizeof(struct sigaction));
	sa.sa_handler = SIG_DFL;
	sigemptyset(&sa.sa_mask);
	sigaction(SIGSEGV, &sa, NULL);
	sigaction(SIGBUS,  &sa, NULL);
	DBG(printf(" >> Unsetted.\n"));
}

void my_set_SIGSEGV_handler(void) {
	DBG(printf(" >> My SIGSEGV Handler 1 - rerun routine engaged!\n"));
	struct sigaction sa;
	memset(&sa, 0, sizeof(struct sigaction));
	sa.sa_handler = my_action;
	sigemptyset(&sa.sa_mask);
	sigaction(SIGSEGV, &sa, NULL);
	sigaction(SIGBUS,  &sa, NULL);
}

double Matrix_Sum(const gsl_matrix* M) {
	double s=0; int i, j, rows=M->size1, cols=M->size2;
	for(i=0; i<rows; i++) {
		for(j=0; j<cols; j++) {
			s = s + gsl_matrix_get(M, i, j);
		}
	}
	return s;
}

void my_copy_double_to_matrix(double* arr, gsl_matrix* M) {
	int i,j,k;
	for(i=0; i<M->size1; i++) {
		for(j=0; j<M->size2; j++) {
			k=i*(M->size2) + j;
			gsl_matrix_set(M, i, j, arr[k]);
		}
	}
}


void Print_Matrix(CBLAS_TRANSPOSE_t Trans, const gsl_matrix* M, const char* info) {
	int i,j;
	DBG(printf("Print_Matrix"));
	if(Trans==CblasNoTrans) { printf(">> Matrix [%s] (%d x %d)\n", info, (int)M->size1, (int)M->size2); 
		for(i=0; i<M->size1; i++) {
			for(j=0; j<M->size2; j++) { printf("%g ", gsl_matrix_get(M, i, j)); }
			printf("\n");
		}
	}
	else { printf(">> Matrix [%s] (%d x %d)\n", info, (int)M->size2, (int)M->size1); 
		for(j=0; j<M->size2; j++) {
			for(i=0; i<M->size1; i++) { printf("%g ", gsl_matrix_get(M, i, j)); }
			printf("\n");
		}
	}
}

void Print_Vector(const gsl_vector* V, const char* info) {
	printf(">> Vector [%s] (%lu x 1)\n", info, V->size); 
	int i;
	for(i=0; i<V->size; i++) printf("%g ", gsl_vector_get(V,i));
	printf("\n");
	return;
}

#define STOPWATCH_INCREMENT(id) \
	if(g_ptv[id]==NULL) g_ptv[id]=&g_tvTick[id]; \
	struct timeval* ptv_late  = g_ptv[id]; \
	struct timeval* ptv_early = (g_ptv[id]==&g_tvTick[id]) ? &g_tvTock[id] : &g_tvTick[id]; \
	g_ptv[id] = ptv_early; \
	gettimeofday(ptv_late, NULL); \
	if(ptv_early->tv_usec != 0 && ptv_early->tv_sec != 0){ \
		unsigned long deltaMicros = 1000000*(ptv_late->tv_sec - ptv_early->tv_sec) + ptv_late->tv_usec - ptv_early->tv_usec; \
		g_time[id] += deltaMicros; \
	}


int Is_GSL_Vector_Equal(const gsl_vector* A, const gsl_vector* B) /* Why? gsl_vector_equal is not available. */
{
	if(A->size!=B->size) { DBG(printf("Vector A and B size not equal.\n")); return 0;}
	double a, b;
	int i; for(i=0; i<A->size; i++) {
		a=gsl_vector_get(A, i); b=gsl_vector_get(B, i);
		if(a==b) continue; // 04-19-12 Could this be the performance culprit ...
		double ab = a-b;
		double r = (double)ab/(double)b; if(r<0) r=-r;
		if(r > is_equal_knob) { // Changed on 04-20-2012
			DBG(printf(" A[%d]!=B[%d], %f and %f, (a-b)/b=%f, a-b=%f, knob=%g\n", i, i, a, b, r, ab, is_equal_knob));
			return 0; }
	}
	DBG(printf(" A[0]=%f, B[0]=%f\n", A->data[0], B->data[0]));
	return 1;
}

int Is_GSL_Matrix_Equal(const gsl_matrix* A, const gsl_matrix* B)
{
	int i, j;
	double aa, bb;
	DBG(printf("[Is_Matrix_Equal] A:%lu x %lu; B:%lu x %lu\n", A->size1, A->size2, B->size1, B->size2));
	if(!((A->size1 == B->size1)&&(A->size2 == B->size2))) { 
		DBG(printf("[Is_Matrix_Equal] Matrices A and B are of different dimensions!\n"));
		DBG(printf("                  A:%lu x %lu; B:%lu x %lu\n", A->size1, A->size2, B->size1, B->size2));
		return 0;
	}
	for(i=0; i<A->size1; i++) {
		for(j=0; j<A->size2; j++) {
			aa = gsl_matrix_get(A, i, j);
			bb = gsl_matrix_get(B, i, j);
			if(aa==bb) continue; // 04-19-12 Could this be the performance culprit ...
			double ab = aa-bb;
			double r = (double)ab/(double)bb; if(r<0) r=r*-1;
			if(r > is_equal_knob) { 
			DBG(printf(" A[%d,%d]!=B[%d,%d], %g and %g, (a-b)/b=%f, knob=%g\n", i, j, i, j, aa, bb, r, is_equal_knob));
			return 0; }
		}
	}
	return 1;
}

double Get_GSL_Matrix_RMSD(const gsl_matrix* A, const gsl_matrix* B) {
	int s1, s2;
	s1 = A->size1; s2 = A->size2; // Pretend A and B are of same size
	if(!((A->size1 == B->size1) && (A->size2 == B->size2))) {
		DBG(printf("[Get_GSL_Matrix_RMSD] A's size (%ldx%ld) != B's size(%ldx%ld).\n",
		A->size1, A->size2, B->size1, B->size2));
		s1 = (A->size1 > B->size1 ? B->size1 : A->size1);
		s2 = (A->size2 > B->size2 ? B->size2 : A->size2);
		DBG(printf("                      Assume size is %dx%d.\n", s1, s2));
	}
	
	double rsd = 0;
	int i, j; for(i=0; i<s1; i++) {
		for(j=0; j<s2; j++) {
			double a_elem = gsl_matrix_get(A, i, j);
			double b_elem = gsl_matrix_get(B, i, j);
			a_elem -= b_elem;
			rsd += (a_elem * a_elem);
		}
	}
	rsd *= (1.0 / s1 / s2);
	return rsd;
}

double Get_GSL_Vector_RMSD(const gsl_vector* X, const gsl_vector* Y) {
	int s;
	s = X->size;
	if(X->size != Y->size) {
		DBG(printf("[Get_GSL_Vector_RMSD] X's size (%ld) != Y's size (%ld).\n",
			X->size, Y->size));
		s = (X->size > Y->size) ? Y->size : X->size;
		DBG(printf("                      Assume size is %d.\n", s));
	}
	double rsd = 0;
	int i;
	for(i=0; i<s; i++) {	
		double x_elem = gsl_vector_get(X, i);
		double y_elem = gsl_vector_get(Y, i);
		x_elem -= y_elem;
		rsd += x_elem * x_elem;
	}
	rsd *= (1.0 / s);
	return rsd;
}

/* Validating if C2=A*B+C. Alpha and Beta both need be set to 1. */
/* Returns 1 for equal and 0 for not equal. */
noinline 
int Is_GSL_MM_Equal_actual(CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB,
		    double alpha, const gsl_matrix* A, const gsl_matrix* B,
		    double beta,  const gsl_matrix* C, const gsl_matrix* C2) {
	int result;
	my_stopwatch_checkpoint(4);
FTV_REAL_TRY(0) {
	int A_rc = (TransA == CblasNoTrans) ? A->size1 : A->size2;
	int A_cc = (TransA == CblasNoTrans) ? A->size2 : A->size1;
	int B_rc = (TransB == CblasNoTrans) ? B->size1 : B->size2;
	int C_rc = C->size1;
	int C_cc = C->size2;
	if(!((A_cc==B_rc)&&(A_rc==C_rc))) {
		DBG(printf("Size a != size b!\n")); 
		my_stopwatch_stop(4);
		return 0; }
	srand(time(0));
	gsl_vector* r = gsl_vector_alloc(C_cc);
	gsl_vector* C2r, *Br, *ABr;
	C2r = gsl_vector_alloc(A_rc); gsl_vector_set_zero(C2r);
	Br = gsl_vector_alloc(B_rc); gsl_vector_set_zero(Br);
	ABr = gsl_vector_alloc(A_rc); gsl_vector_set_zero(ABr);
	DBG(printf("sizeof r, C2r, Br, ABr: %ld %ld %ld %ld\n", r->size, C2r->size, Br->size, ABr->size));
	int i;
	for(i=0; i<r->size; i++) gsl_vector_set(r, i, (double)/*rand()*/i/100);
	my_dgemv(TransB, alpha, B, r, 0, Br);  // Br <- B * r
	my_dgemv(TransA, 1, A, Br, 0, ABr);    // ABr <- A * B * r
	my_dgemv(CblasNoTrans, beta, C, r, 1, ABr); // ABr <- ABr + C * r
	my_dgemv(CblasNoTrans, 1, C2, r, 0, C2r);
	result = Is_GSL_Vector_Equal(C2r, ABr);
	gsl_vector_free(C2r);
	gsl_vector_free(Br);
	gsl_vector_free(ABr);
	gsl_vector_free(r);
	my_stopwatch_stop(4);
	if(result==0) { DBG(printf("## DGEMM not equal\n")); }
	else { DBG(printf("[Is_GSL_MM_Equal] DGEMM is equal.\n")); }
} FTV_REAL_CATCH(0) {} FTV_REAL_END(0);
	return result;
}

noinline
int Is_GSL_MM_Equal(CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB,
		    double alpha, const gsl_matrix* A, const gsl_matrix* B,
		    double beta,  const gsl_matrix* C, const gsl_matrix* C2) {
	MY_SET_SIGSEGV_HANDLER();
	int jmpret;
	SUPERSETJMP("Is_GSL_MM_Equal");
	trick_me_jr(jmpret);
	return Is_GSL_MM_Equal_actual(TransA, TransB, alpha, A, B, beta, C, C2);
}

/* Validating if C2=A*t(A)+B. Same as MM. */
noinline 
int Is_GSL_DSYRK_Equal_actual(const CBLAS_UPLO_t uplo, CBLAS_TRANSPOSE_t trans,
		       double alpha, const gsl_matrix* A,
		       double beta,  const gsl_matrix* C, const gsl_matrix* C2) {
	my_stopwatch_checkpoint(4);
	int i, j, len;
	if(C->size1 != C->size2) { printf("[Is_GSL_DSYRK_Equal] Matrix C is not square!\n"); 
		my_stopwatch_stop(4);
		return 0; }
	len = C->size1;
	DBG(printf("[Is_GSL_DSYRK_Equal]\n"));
	gsl_matrix* tA = gsl_matrix_alloc(A->size2, A->size1);
	gsl_matrix* uploC = gsl_matrix_alloc(C->size1, C->size2);
	if(uplo == CblasUpper) { /* Fixed on 2012-08-21 */
		for(i=0; i<len; i++) {
			for(j=i; j<len; j++)
				gsl_matrix_set(uploC, i, j, gsl_matrix_get(C, i, j)*beta);
		}
		for(i=1; i<len; i++) {
			for(j=0; j<i; j++)
				gsl_matrix_set(uploC, i, j, gsl_matrix_get(C, i, j));
		}
	} else {
		for(i=0; i<len; i++) {
			for(j=i+1; j<len; j++) 
				gsl_matrix_set(uploC, i, j, gsl_matrix_get(C, i, j));
		}
		for(i=0; i<len; i++) {
			for(j=0; j<=i; j++)
				gsl_matrix_set(uploC, i, j, gsl_matrix_get(C, i, j)*beta);
		}
	}
	gsl_matrix_transpose_memcpy(tA, A);

	gsl_matrix* AAt = gsl_matrix_alloc(len, len);
FTV_REAL_TRY(0) {
	if(trans==CblasNoTrans) {
//		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, tA, 0.0, AAt); // I think this causes problems
		my_dgemm(CblasNoTrans, CblasNoTrans, alpha, A, tA, 0.0, AAt);
	} else {
//		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tA, A, 0.0, AAt);
		my_dgemm(CblasNoTrans, CblasNoTrans, alpha, tA, A, 0.0, AAt);
	}
} FTV_REAL_CATCH(0) {} FTV_REAL_END(0);
	if(uplo == CblasUpper) {
		for(i=1; i<len; i++) {
			for(j=0; j<i; j++) {
//				gsl_matrix_set(uploC, i, j, 0);
				gsl_matrix_set(AAt,   i, j, 0);
			}
		}
	} else if(uplo == CblasLower) {
		for(i=0; i<len; i++) {
			for(j=i+1; j<len; j++) {
//				gsl_matrix_set(uploC, i, j, 0);
				gsl_matrix_set(AAt,   i, j, 0);
			}
		}
	}
	gsl_matrix_add(uploC, AAt);
	int result = Is_GSL_Matrix_Equal(uploC, C2);
	gsl_matrix_free(uploC);
	gsl_matrix_free(AAt);
	my_stopwatch_stop(4);
	if(result==0) { DBG(printf("## DSYRK not equal\n")); }
	return result;
}

noinline
int Is_GSL_DSYRK_Equal(const CBLAS_UPLO_t uplo, CBLAS_TRANSPOSE_t trans,
		       double alpha, const gsl_matrix* A,
		       double beta,  const gsl_matrix* C, const gsl_matrix* C2) {
	MY_SET_SIGSEGV_HANDLER();
	return Is_GSL_DSYRK_Equal_actual(uplo, trans, alpha, A, beta, C, C2);
}

noinline
int Is_GSL_DSYRK_Equal2_actual(const CBLAS_UPLO_t uplo, CBLAS_TRANSPOSE_t trans,
		       double alpha, const gsl_matrix* A,
		       double beta,  const gsl_matrix* C, const gsl_matrix* C2)
{
	DBG(printf("[Is_GSL_DSYRK_Equal2]\n"));
	my_stopwatch_checkpoint(4);
	if(C2->size1 != C2->size2) { DBG(printf("[Is_GSL_DSYRK_Equal2] C2 is not square. \n")); return 0; }
	int n = C2->size1, i, j;
	if(trans==CblasNoTrans) {
		if(!((C->size1==n) && (C->size2==n) && (A->size1==n))) { DBG(printf("[Is_GSL_DSYRK_Equal2] Dimensions do not agree.\n")); return 0; }
	} else {
		if(!((C->size1==n) && (C->size2==n) && (A->size2==n))) { DBG(printf("[Is_GSL_DSYRK_Equal2] Dimensions do not agree.\n")); return 0; }
	}

	gsl_matrix* c;
	if(beta!=0) {
		c = gsl_matrix_alloc(n, n);
		for(i=0; i<n; i++) {
			for(j=0; j<n; j++) {gsl_matrix_set(c, i, j, gsl_matrix_get(C2, i, j) - beta * gsl_matrix_get(C, i, j));}
		}
	}
	int ret;
//FTV_REAL_TRY(0) {
	if(beta != 0) {
		if(uplo==CblasUpper) ret = check_aat_triangle(c, uplo, 0, n-1, n, alpha, A, trans);
		if(uplo==CblasLower) ret = check_aat_triangle(c, uplo, n-1, 0, n, alpha, A, trans);
		gsl_matrix_free(c);
	} else {
		if(uplo==CblasUpper) ret = check_aat_triangle(C2, uplo, 0, n-1, n, alpha, A, trans);
		if(uplo==CblasLower) ret = check_aat_triangle(C2, uplo, n-1, 0, n, alpha, A, trans);	
	}
//} FTV_REAL_CATCH(0) {} FTV_REAL_END(0);
	my_stopwatch_stop(4);
	if(ret==0) { DBG(printf("## DSYRK not equal\n")); }
	return ret;
}

/* Check if Y2 = A * X */
noinline 
int Is_GSL_DGEMV_Equal_actual(CBLAS_TRANSPOSE_t Trans, double alpha, const gsl_matrix* A, const gsl_vector* X, double beta, const gsl_vector* Y, const gsl_vector* Y2) {
	/* Y2 is the output alphaAX+betaY */
	my_stopwatch_checkpoint(4);
	DBG(printf(" >> \n"));
	DBG(if(X==NULL) printf("X is null\n"););
	DBG(if(Y2==NULL) printf("Y2 is null\n"););
	if(Y->size != Y2->size) { DBG(printf("[Is_GSL_DGEMV_Equal]Error: X's size != Y's size\n"));
		my_stopwatch_stop(4);
		return 0; 
	}
	int i, i1, i2, j, j1, j2, result;
	double sum1_1=0, sum1=0, sum2=0;
	double abserr, relerr;
FTV_REAL_TRY(0) {
	if(beta!=0) {
		#ifndef FT_CHKR
		for(i=0; i<Y->size; i++) {
		#else
		for(i=0, i1=1, i2=2; i<Y->size; i++, i1++, i2++,
			PROTECT_IDX_I) {
		#endif
//			sum1 = sum1 + gsl_vector_get(Y, i);
			sum1 = sum1 + Y->data[i];
		}
		
		sum1 = sum1 * beta;
	}
	#ifndef FT_CHKR
	for(i=0; i<X->size; i++) {
	#else
	for(i=0, i1=1, i2=2; i<X->size; i++, i1++, i2++,
		PROTECT_IDX_I) {
	#endif
		double currCS = 0; /* current column Sum */
		int jj = (Trans==CblasNoTrans) ? A->size1 : A->size2;
		#ifndef FT_CHKR
		for(j=0; j<jj; j++) {
		#else
		for(j=0, j1=1, j2=2; j<jj; j++, j1++, j2++,
			PROTECT_IDX_J) {
		#endif
			if(Trans==CblasNoTrans) { 
				currCS = currCS + A->data[j*A->tda+i];
				//gsl_matrix_get(A, j, i);
			} else {
				currCS = currCS + A->data[i*A->tda+j]; 
				//gsl_matrix_get(A, i, j);
			}
		}
		sum1_1 = sum1_1 + currCS * X->data[i];//gsl_vector_get(X, i);
	}
	sum1_1 *= alpha;
	sum1 += sum1_1;
	#ifndef FT_CHKR
	for(i=0; i<Y->size; i++) {
	#else
	for(i=0,i1=1,i2=2; i<Y->size; i++,i1++,i2++,
		PROTECT_IDX_I) {
	#endif
		sum2 = sum2 + Y2->data[i];//gsl_vector_get(Y2, i); 
	}
	abserr = sum1 - sum2; if(abserr < 0) abserr = -abserr;
	relerr = abserr / sum2; if(relerr < 0) relerr = -relerr;
	const double tol = FT_TOLERANCE / 100;
	{ DBG(printf("[Is_GSL_DGEMV_Equal] DGEMV Equal, abs.err=%g, rel.err=%g (tol=%g)\n", abserr, relerr, tol)); }
	if(abserr > tol && relerr > tol) {
		result = 0; /* Precision problem!!! */
	}
	else result = 1;
	my_stopwatch_stop(4);
	if(result==0) { DBG(printf("[Is_GSL_DGEMV_Equal] DGEMV not equal\n")); }
} FTV_REAL_CATCH(0) {} FTV_REAL_END(0);
	return result;
}

noinline
int Is_GSL_DSYRK_Equal2(const CBLAS_UPLO_t uplo, CBLAS_TRANSPOSE_t trans,
		       double alpha, const gsl_matrix* A,
		       double beta,  const gsl_matrix* C, const gsl_matrix* C2)
{
	MY_SET_SIGSEGV_HANDLER();
	int jmpret;
	SUPERSETJMP("Is_GSL_DSYRK_Equal2");
	trick_me_jr(jmpret);
	return Is_GSL_DSYRK_Equal2_actual(uplo, trans, alpha, A, beta, C, C2);
}

noinline
int Is_GSL_DGEMV_Equal(CBLAS_TRANSPOSE_t Trans, double alpha, const gsl_matrix* A, const gsl_vector* X, double beta, const gsl_vector* Y, const gsl_vector* Y2) {
	MY_SET_SIGSEGV_HANDLER();
	return Is_GSL_DGEMV_Equal_actual(Trans, alpha, A, X, beta, Y, Y2);
}


/* Check dtrsv routine, should be inv(op(A)) * X2 == X */
int Is_GSL_DTRSV_Equal(CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t diag, const gsl_matrix* A, const gsl_vector* X, const gsl_vector* X2) {
	my_stopwatch_checkpoint(4);
	if(A->size1 != A->size2) { printf("[Is_GSL_DTRSV_Equal] Matrix A is not square.\n"); return -1; }
	if(A->size1 != X->size)  { printf("[Is_GSL_DTRSV_Equal] Sizes of A and X don't match.\n"); return -1; }
	if(X->size  != X2->size) { printf("[Is_GSL_DTRSV_Equal] Sizes of X and X2 don't match.\n");return -1; }
	gsl_matrix* AA = gsl_matrix_alloc(A->size1, A->size2);
	gsl_matrix_memcpy(AA, A);
	int i,j;
	if(Uplo==CblasUpper) {
		for(i=1; i<A->size1; i++) {
			for(j=0; j<i; j++) gsl_matrix_set(AA, i, j, 0.0);
		}
	} else if(Uplo==CblasLower) {
		for(i=0; i<A->size1-1; i++) {
			for(j=i+1; j<A->size1; j++) gsl_matrix_set(AA, i, j, 0.0);
		}
	}
	if(diag==CblasUnit) { for(i=0; i<A->size1; i++) gsl_matrix_set(AA, i, i, 1.0); }
	gsl_vector* XX = gsl_vector_alloc(X->size);
	gsl_vector_set_zero(XX);
//	Print_Vector(X,  "[DTRSV Checker]X");
//	Print_Vector(X2, "[DTRSV Checker]X2");
	my_dgemv(TransA, 1.0, AA, X2, 0.0, XX);
//	Print_Vector(XX, "[DTRSV Checker]Result of MV multiplication");
	int result = Is_GSL_Vector_Equal(XX, X);
	gsl_vector_free(XX); gsl_matrix_free(AA);
	my_stopwatch_stop(4);
	if(result==0) printf("## DTRSV not equal\n");
	return result;
}

noinline 
int Is_GSL_linalg_cholesky_decomp_Equal_actual(const gsl_matrix* A, const gsl_matrix* A_out) {
	if(A->size1 != A->size2) { printf("[Is_GSL_linalg_cholesky_decomp_Equal] A is not square -- not to say symmetric.\n"); return 0; }
	if(A_out->size1 != A_out->size2) { printf("[Is_GSL_linalg_cholesky_decomp_Equal] A_out is not square -- not to say symmetric.\n"); return 0; }
	if(A->size1 != A_out->size1) { printf("[Is_GSL_linalg_cholesky_decomp_Equal] A and A_out's dimensions do not agree.\n"); return 0; }
	int n = A->size1, i, j, k;

//	my_dgemm(CblasNoTrans, CblasTrans, 1.0, L, L, 0.0, LU); // This is too slow. May take up to 100x time of decomposition!
//	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, L, L, 0.0, LU);  // Confirmed. Executing this on my laptop with an OPTIMIZED gsl library reduces overhead.
								       // 20x faster than unoptimized -- but still not fast enough.
	int ret = 1;
FTV_REAL_TRY(0) {
	for(i=0; i<n; i++) {
		for(j=0; j<=i; j++) {
			double tmp = 0;
			double tmp0 = gsl_matrix_get(A, i, j);
			for(k=0; k<=j; k++) {
				tmp = tmp + gsl_matrix_get(A_out, i, k)*gsl_matrix_get(A_out, j, k);
			}
			if(tmp == tmp0) continue;
			else {
				tmp=tmp-tmp0; double r = tmp/tmp0; if(r<0) r=r*-1.0;
				if(r > FT_TOLERANCE) { ret = 0; goto re; }
			}
		}
	}
	re:
	if(ret==0) printf("## Cholesky_Decomp not equal\n");
} FTV_REAL_CATCH(0) {} FTV_REAL_END(0);
	return ret;
}

noinline
int Is_GSL_linalg_cholesky_decomp_Equal(const gsl_matrix* A, const gsl_matrix* A_out) {
	MY_SET_SIGSEGV_HANDLER();
	int jmpret;
	SUPERSETJMP("Is_GSL_linalg_cholesky_decomp_Equal");
	return Is_GSL_linalg_cholesky_decomp_Equal_actual(A, A_out);
}

noinline 
int check_aat_triangle_actual(const gsl_matrix* C, CBLAS_UPLO_t uplo, int row, int col, int len, double alpha, const gsl_matrix* A, CBLAS_TRANSPOSE_t trans) {
//FTV_REAL_TRY(0) {
	DBG(printf("check_aat_triangle [%d, %d] len=%d\n", row, col, len));
	/* Let's say, we have a "body" and two "flanks" */
	if(C->size1 != C->size2) { printf("[check_aat_triangle_actual] !! C is not square.\n"); return 0; }
	int body_rs, body_re, body_cs, body_ce, body_len;
	int uflank_row, uflank_col, dflank_row, dflank_col, flank_len;
	if(len <= 2) {
		if(len==1) {
			int ret = check_aat_elem(C, row, col, alpha, A, trans);
			return ret;
		} else {
			int ret;
			ret = check_aat_elem(C, row, col, alpha, A, trans); if(ret==0) return 0;
			if(uplo==CblasUpper) {
				ret = check_aat_elem(C, row, col-1, alpha, A, trans); if(ret==0) return 0;
				ret = check_aat_elem(C, row+1, col, alpha, A, trans); if(ret==0) return 0;
				return 1;
			} else {
				ret = check_aat_elem(C, row, col+1, alpha, A, trans); if(ret==0) return 0;
				ret = check_aat_elem(C, row-1, col, alpha, A, trans); if(ret==0) return 0;
				return 1;
			}
		}
	} else {
		body_len = (int)((len+1) / 2); flank_len = (int)(len / 2);
		int ret;
		if(uplo==CblasUpper) {
			body_rs = row; body_re = row + body_len - 1;
			body_cs = col - body_len + 1; body_ce = col;
			uflank_row = row; uflank_col = body_cs - 1;
			dflank_row = body_re + 1; dflank_col = body_ce;
		} else {
			body_rs = row - body_len + 1; body_re = row;
			body_cs = col; body_ce = col + body_len - 1;
			uflank_row = body_rs - 1; uflank_col = col;
			dflank_row = body_re; dflank_col = body_ce + 1;
		}
		DBG(printf("body rows: %d ~ %d, cols: %d ~ %d, flanks at: (%d, %d) and (%d, %d)\n", body_rs, body_re, body_cs, body_ce, uflank_row, uflank_col, dflank_row, dflank_col));
		ret = check_aat_submatrix(C, body_rs, body_re, body_cs, body_ce, alpha, A, trans); if(ret==0) { DBG(printf("@")); return 0; }
		ret = check_aat_triangle(C, uplo, uflank_row, uflank_col, flank_len, alpha, A, trans); if(ret==0) { DBG(printf("#")); return 0; }
		ret = check_aat_triangle(C, uplo, dflank_row, dflank_col, flank_len, alpha, A, trans); if(ret==0) { DBG(printf("$")); return 0; }
		return 1;
	}
//} FTV_REAL_CATCH(0) {} FTV_REAL_END(0);
	return 0;
}

noinline
int check_aat_triangle(const gsl_matrix* C, CBLAS_UPLO_t uplo, int row, int col, int len, double alpha, const gsl_matrix* A, CBLAS_TRANSPOSE_t trans) {
	MY_SET_SIGSEGV_HANDLER();
	return check_aat_triangle_actual(C, uplo, row, col, len, alpha, A, trans);
}

#define AAT_CHECK_DIM \
	int n = C->size1; \
	if(C->size1 != C->size2) { printf("C is not square.\n"); return 0; }

noinline 
int check_aat_elem_actual(const gsl_matrix* C, int row, int col, double alpha, const gsl_matrix* A, CBLAS_TRANSPOSE_t trans) {
FTV_REAL_TRY(0) {
	AAT_CHECK_DIM; int i;
	double temp=0, tmp0 = gsl_matrix_get(C, row, col);
	if(trans==CblasNoTrans) {
		n=A->size2;
		for(i=0; i<A->size2; i++)  { temp = temp + gsl_matrix_get(A, row, i) * gsl_matrix_get(A, col, i) * alpha;}
	} else {
		n=A->size1;
		for(i=0; i<A->size1; i++) temp = temp + gsl_matrix_get(A, i, row) * gsl_matrix_get(A, i, col) * alpha;
	}

	if(temp != tmp0) {
		double r = (temp-tmp0)/tmp0; if(r<0.0) r=r*-1.0;
		if(r > FT_TOLERANCE) {
			printf("check_aat_elem(%d, %d)=0\n", row, col);
			printf("%f vs %f\n", temp, tmp0);
			return 0; }
	}
} FTV_REAL_CATCH(0) {} FTV_REAL_END(0);
	return 1;
}

int check_aat_elem(const gsl_matrix* C, int row, int col, double alpha, const gsl_matrix* A, CBLAS_TRANSPOSE_t trans) {
	MY_SET_SIGSEGV_HANDLER();
	return check_aat_elem_actual(C, row, col, alpha, A, trans); 
}

/* This is the routine to be integrated into tommy.h */
noinline 
int check_aat_submatrix_actual(const gsl_matrix* C, int rowstart, int rowend, int colstart, int colend, double alpha, const gsl_matrix* A, CBLAS_TRANSPOSE_t trans) {
	int ret;
FTV_REAL_TRY(0) {
	int n = C->size1;
	if(C->size1 != C->size2) { printf("C is not square.\n"); return 0; }
	gsl_vector* v = gsl_vector_alloc(n);
	int m = (trans==CblasNoTrans) ? A->size2 : A->size1;
	int i, j; for(i=0; i<n; i++) { gsl_vector_set(v, i, (double)rand()); }
	gsl_vector* aatv = gsl_vector_alloc(n), *cv=gsl_vector_alloc(n), *atv=gsl_vector_alloc(m);
	gsl_vector_set_zero(aatv); gsl_vector_set_zero(cv); gsl_vector_set_zero(atv);
	for(i=rowstart; i<=rowend; i++) {
		double tmp = 0;
		for(j=colstart; j<=colend; j++) {
			tmp = tmp + gsl_vector_get(v, j) * gsl_matrix_get(C, i, j);
		}
		gsl_vector_set(cv, i, tmp);
	}

	for(i=0; i<m; i++) {
		double tmp = 0;
		if(trans == CblasNoTrans) 
		{ for(j=colstart; j<=colend; j++) tmp=tmp+gsl_vector_get(v, j)*gsl_matrix_get(A, j, i); } // Fixed on 08-21-2012
		else
		{ for(j=colstart; j<=colend; j++) tmp=tmp+gsl_vector_get(v, j)*gsl_matrix_get(A, i, j); } // Fixed on 08-21-2012
		gsl_vector_set(atv, i, tmp);
	}
	for(i=rowstart; i<=rowend; i++) {
		double tmp = 0;
		if(trans == CblasNoTrans)
		{ for(j=0; j<m; j++) tmp=tmp+gsl_vector_get(atv, j)*gsl_matrix_get(A, i, j)*alpha; }
		else
		{ for(j=0; j<m; j++) tmp=tmp+gsl_vector_get(atv, j)*gsl_matrix_get(A, j, i)*alpha; }
		gsl_vector_set(aatv, i, tmp);
	}
	ret = Is_GSL_Vector_Equal(aatv, cv);
	free(cv); free(aatv); free(atv); free(v);
} FTV_REAL_CATCH(0) {} FTV_REAL_END(0);
	return ret;
}

noinline
int check_aat_submatrix(const gsl_matrix* C, int rowstart, int rowend, int colstart, int colend, double alpha, const gsl_matrix* A, CBLAS_TRANSPOSE_t trans) {
	MY_SET_SIGSEGV_HANDLER();
	return check_aat_submatrix_actual(C, rowstart, rowend, colstart, colend, alpha, A, trans);
}

int Print_Hello(void) {
	printf("haha! ver 03-15\n");
	return 0;
}

noinline
void GSL_BLAS_DGEMM_FT3(CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, 
			double alpha, const gsl_matrix* matA, const gsl_matrix* matB, 
			double beta, gsl_matrix* matC)

{
	/* Protect the pointers to matrix A, B, (the original) C */
	unsigned long matA_0, matA_1, matA_2, 
	     matB_0, matB_1, matB_2,
	     matC_0, matC_1, matC_2;
	TRIPLICATE(matA, matA_0, matA_1, matA_2);
	TRIPLICATE(matB, matB_0, matB_1, matB_2);

	/* Protect the tda, size1 and size2 of matC_bak */
	size_t mcb0, mcb1, mcb2;
	size_t mcs10, mcs11, mcs12;
	size_t mcs20, mcs21, mcs22;
	size_t mbs10, mbs11, mbs12;
	size_t mbs20, mbs21, mbs22;
	size_t mas10, mas11, mas12;
	size_t mas20, mas21, mas22;
	TRIPLICATE_SIZE_T(matB->size1, mbs10, mbs11, mbs12);
	TRIPLICATE_SIZE_T(matB->size2, mbs20, mbs21, mbs22);
	TRIPLICATE_SIZE_T(matA->size1, mas10, mas11, mas12);
	TRIPLICATE_SIZE_T(matA->size2, mas20, mas21, mas22);
	
	double sumA, sumB, sumC;
	sumA=my_sum_matrix(matA); sumB=my_sum_matrix(matB); sumC=my_sum_matrix(matC); 
	/* Protect the pointer to ECC codes. */
	void* ecMatA, *ecMatB, *ecMatC, *ecMatC_2;
	unsigned long ema_0, ema_1, ema_2, emb_0, emb_1, emb_2, emc_0, emc_1, emc_2;
	
	nonEqualCount=0;
	MY_SET_SIGSEGV_HANDLER();
	gsl_matrix* matC_bak = gsl_matrix_alloc(matC->size1, matC->size2);
	/* Protect matC_bak -> data */
	unsigned long mCbd_0, mCbd_1, mCbd_2;
	TRIPLICATE(matC_bak, matC_0, matC_1, matC_2);
	TRIPLICATE(matC_bak->data, mCbd_0, mCbd_1, mCbd_2);
	gsl_matrix_memcpy(matC_bak, matC);
	TRIPLICATE_SIZE_T(matC_bak->tda, mcb0, mcb1, mcb2);
	TRIPLICATE_SIZE_T(matC_bak->size1, mcs10, mcs11, mcs12);
	TRIPLICATE_SIZE_T(matC_bak->size2, mcs20, mcs21, mcs22);
	int jmpret, isEqual;
	encode((matB->data), (matB->size1*matB->size2), &ecMatB);
	TRIPLICATE(ecMatB, emb_0, emb_1, emb_2);
	int emc_size = encode((matC_bak->data), (matC->size1*matC->size2), &ecMatC);
	ecMatC_2 = malloc(sizeof(double) * emc_size);
	memcpy(ecMatC_2, ecMatC, sizeof(double) * emc_size);
	TRIPLICATE(ecMatC, emc_0, emc_1, emc_2);
	encode((matA->data), (matA->size1*matA->size2), &ecMatA); 
	TRIPLICATE(ecMatA, ema_0, ema_1, ema_2);
	DBG(printf("[DGEMM_FT3] ECC Code Addr: A=%lx, B=%lx, C=%lx\n", (unsigned long)ecMatA, (unsigned long)ecMatB, (unsigned long)ecMatC));
	ResetIsEqualKnob();
	SUPERSETJMP("After encoding");
	if(jmpret != 0) {
		kk:
		/* Chasing the wind, only */
		DBG(printf("[DGEMM_FT3]Recovering input...\n"));
		{
			/* Recover ptr to matA and matA's ecc */
			TRI_RECOVER(matA_0, matA_1, matA_2);
			size_t *mas1 = &(((gsl_matrix*)matA)->size1),
			       *mas2 = &(((gsl_matrix*)matA)->size2);
			TRI_RECOVER_SIZE_T((*mas1), mas10, mas11, mas12);
			TRI_RECOVER_SIZE_T((*mas2), mas20, mas21, mas22);
			TRI_RECOVER(ema_0, ema_1, ema_2);
			if((long)ecMatA != ema_0) ecMatA = (void*)ema_0;
			MY_MAT_CHK_RECOVER_POECC(sumA, ecMatA, (gsl_matrix*)matA_0);

			/* Recover ptr to matB and matB's ecc */
			TRI_RECOVER(matB_0, matB_1, matB_2);
			size_t *mbs1 = &(((gsl_matrix*)matB)->size1),
			       *mbs2 = &(((gsl_matrix*)matB)->size2);
			TRI_RECOVER_SIZE_T((*mbs1), mbs10, mbs11, mbs12);
			TRI_RECOVER_SIZE_T((*mbs2), mbs20, mbs21, mbs22);
			TRI_RECOVER(emb_0, emb_1, emb_2);
			if((long)ecMatB != emb_0) ecMatB = (void*)emb_0;
			MY_MAT_CHK_RECOVER_POECC(sumB, ecMatB, (gsl_matrix*)matB_0);

			/* Recover ptr to mat C and matC's ecc */
			/* Mat C is special: its size1 and size2 are also 
			 * vulnerable to soft errors therefore needs extra care
			 */
			TRI_RECOVER(matC_0, matC_1, matC_2);
			TRI_RECOVER(mCbd_0, mCbd_1, mCbd_2);
			if(matC_bak->data != (double*)mCbd_0)
				matC_bak->data = (double*)mCbd_0;
			TRI_RECOVER_SIZE_T(matC_bak->size1, mcs10, mcs11, mcs12);
			TRI_RECOVER_SIZE_T(matC_bak->size2, mcs20, mcs21, mcs22);
			TRI_RECOVER_SIZE_T(matC->size1, mcs10, mcs11, mcs12);
			TRI_RECOVER_SIZE_T(matC->size2, mcs20, mcs21, mcs22);
			TRI_RECOVER_SIZE_T(matC_bak->tda, mcb0, mcb1, mcb2);
			TRI_RECOVER(emc_0, emc_1, emc_2);
			if((long)ecMatC != emc_0) ecMatC = (void*)emc_0;
			int n_retry = 0;
chk_rec_c:
			;  // <---------- Used to fix the bug with llvm-gcc. Compiling with llvm-g++ won't have this bug.
			float fr = MY_MAT_CHK_RECOVER_POECC(sumC, ecMatC, (gsl_matrix*)matC_0);
			if(fr > 0.5 && n_retry < 5) {
				DBG(printf("[DGEMM_FT3] Oh! Is ecMatC damaged? (%d/%d)\n", n_retry, 5));
				memcpy(ecMatC, ecMatC_2, sizeof(double) * emc_size);
				n_retry++;
				goto chk_rec_c;
			}

			//And there we have completed recovering it....
			DBG(printf("[DGEMM_FT3]Copying back input...\n"));
			int ec;
			if((ec = gsl_matrix_memcpy(matC, matC_bak)) != GSL_SUCCESS) {
				printf("[DGEMM_FT3] ! memcpy returned %d\n", ec);
			}
		}

	}
	DBG(printf("[DGEMM_FT3]Normal call to dgemm.. nonEqualCount=%d\n", nonEqualCount));
	SW3START;
	SUPERSETJMP("Just before dgemm");
	if(jmpret != 0) {
		goto kk;
	}
	gsl_blas_dgemm(TransA, TransB, alpha, matA, matB, beta, matC);
	SW3STOP;
	DBG(printf("[DGEMM_FT3] Check Results..\n"));
	isEqual = Is_GSL_MM_Equal(TransA, TransB, alpha, matA, matB, beta, matC_bak, matC);
	if(isEqual==1) { DBG(printf("[DGEMM_FT3]Result: Equal\n")); }
	else {
		my_stopwatch_checkpoint(6);
		printf("[DGEMM_FT3]Result: NOT Equal\n");
		nonEqualCount = nonEqualCount + 1;
		num_total_retries += 1;
		if(nonEqualCount < NUM_OF_RERUN) {
			printf("[DGEMM_FT3]Restart calculation.\n");
			if((nonEqualCount & 4)==4) AdjustIsEqualKnob(1);
			goto kk;
		}
	}
	ResetIsEqualKnob();
	my_stopwatch_stop(6);
	printf("[DGEMM_FT3]Releasing memory for mat_bak's and deleting temp files\n");
	SUPERSETJMP("Just before free");
	if(jmpret == 0) {
		TRI_RECOVER(matC_0, matC_1, matC_2);
		gsl_matrix_free((gsl_matrix*)matC_0);
		printf("Released.\n");
		
		{
			TRI_RECOVER(ema_0, ema_1, ema_2);
			if((unsigned long)ecMatA != ema_0) ecMatA = (void*)ema_0;
			free(ecMatA);
		}
		{
			TRI_RECOVER(emb_0, emb_1, emb_2);
			if((unsigned long)ecMatB != emb_0) ecMatB = (void*)emb_0;
			free(ecMatB);
		}
		{
			TRI_RECOVER(emc_0, emc_1, emc_2);
			if((unsigned long)ecMatC != emc_0) ecMatC = (void*)emc_0;
			free(ecMatC);
		}
		
		// Last, restore the sizes of elements in case there are corruptions
		{
			size_t *mas1 = &(((gsl_matrix*)matA)->size1),
			       *mas2 = &(((gsl_matrix*)matA)->size2);
			TRI_RECOVER_SIZE_T((*mas1), mas10, mas11, mas12);
			TRI_RECOVER_SIZE_T((*mas2), mas20, mas21, mas22);
			size_t *mbs1 = &(((gsl_matrix*)matB)->size1),
			       *mbs2 = &(((gsl_matrix*)matB)->size2);
			TRI_RECOVER_SIZE_T((*mbs1), mbs10, mbs11, mbs12);
			TRI_RECOVER_SIZE_T((*mbs2), mbs20, mbs21, mbs22);
			TRI_RECOVER_SIZE_T(matC->size1, mcs10, mcs11, mcs12);
			TRI_RECOVER_SIZE_T(matC->size2, mcs20, mcs21, mcs22)
		}

	} else {
		printf("Oops. There is a segfault while trying to free. giving up freeing.");
	}
}

noinline
void GSL_BLAS_DGEMV_FT3(CBLAS_TRANSPOSE_t Trans, double alpha, 
			const gsl_matrix* matA, const gsl_vector* vecX, 
			double beta, gsl_vector* vecY)
{
	/* Protect pointers to matrix A, vectors X and (the original) Y */
	unsigned long matA0, matA1, matA2,
	              vecX0, vecX1, vecX2,
	              vecY0, vecY1, vecY2,
	              vecYbak0, vecYbak1, vecYbak2;
	TRIPLICATE(matA, matA0, matA1, matA2); // A
	TRIPLICATE(vecX, vecX0, vecX1, vecX2); // X
	TRIPLICATE(vecY, vecY0, vecY1, vecY2); // Y
	
	/* Protect size_t of A, X and Y */
	size_t mas10, mas11, mas12;  // A
	size_t mas20, mas21, mas22;  // A
	size_t matda0,matda1,matda2; // A
	size_t vxs0, vxs1, vxs2;  // X
	size_t vxst0,vxst1,vxst2; // X
	size_t vys0, vys1, vys2;  // Y
	size_t vyst0,vyst1,vyst2; // Y
	size_t vybkd0, vybkd1, vybkd2; // Y_bak->data
	TRIPLICATE_SIZE_T(matA->size1, mas10, mas11, mas12);
	TRIPLICATE_SIZE_T(matA->size2, mas20, mas21, mas22);
	TRIPLICATE_SIZE_T(matA->tda,   matda0,matda1,matda2);
	TRIPLICATE_SIZE_T(vecX->size,  vxs0,  vxs1,  vxs2);
	TRIPLICATE_SIZE_T(vecX->stride,vxst0, vxst1, vxst2);
	TRIPLICATE_SIZE_T(vecY->size,  vys0,  vys1,  vys2);
	TRIPLICATE_SIZE_T(vecY->stride,vyst0, vyst1, vyst2);


	MY_SET_SIGSEGV_HANDLER();
	double sumA, sumX, sumY; 
	DBG(printf("[DGEMV_FT3] 1. Get sum of Mat and Vec\n"));
	sumA=my_sum_matrix(matA); sumX=my_sum_vector(vecX); sumY=my_sum_vector(vecY); 

	/* Protect ptrs to ECC codes. */
	void* ecMatA, *ecVecX, *ecVecY;
	unsigned long ema0, ema1, ema2, evx0, evx1, evx2, evy0, evy1, evy2;
	nonEqualCount = 0;

	DBG(printf("[DGEMV_FT3] 2. Encoding input\n"));

	gsl_vector* vecY_bak = gsl_vector_alloc(vecY->size);
	gsl_vector_memcpy(vecY_bak, vecY);
	TRIPLICATE(vecY_bak, vecYbak0, vecYbak1, vecYbak2);
	TRIPLICATE(vecY_bak->data, vybkd0, vybkd1, vybkd2);

	encode((matA->data), (matA->size1*matA->size2), &ecMatA); 
	TRIPLICATE(ecMatA, ema0, ema1, ema2);
	encode((vecX->data), (vecX->size), &ecVecX); 
	TRIPLICATE(ecVecX, evx0, evx1, evx2);
	encode((vecY->data), (vecY->size), &ecVecY); 
	TRIPLICATE(ecVecY, evy0, evy1, evy2);
	
	nonEqualCount=0;
	int jmpret, isEqual;
	DBG(printf("[DGEMV_FT3] 3. Setjmp\n"));
	ResetIsEqualKnob();
	SUPERSETJMP("[DGEMV_FT3] After encoding\n");
	if(jmpret != 0) {
		kk: 
		DBG(printf("[DGEMV_FT3]Recovering input...\n"));
		{
			{ // Recover input data A
				DBG(printf("[DGEMV_FT3]Recovering A ... \n"));
				size_t *mas1 = &(((gsl_matrix*)matA)->size1);
				size_t *mas2 = &(((gsl_matrix*)matA)->size2);
				TRI_RECOVER_SIZE_T((*mas1), mas10, mas11, mas12);
				TRI_RECOVER_SIZE_T((*mas2), mas20, mas21, mas22);
				TRI_RECOVER(matA0, matA1, matA2);
				TRI_RECOVER(ema0, ema1, ema2);
				MY_MAT_CHK_RECOVER_POECC(sumA, (void*)ema0, (gsl_matrix*)matA0);
			}

			{ // Recover input data X
				DBG(printf("[DGEMV_FT3]Recovering X ... \n"));
				TRI_RECOVER(vecX0, vecX1, vecX2);
				size_t *vxs = &(((gsl_vector*)vecX)->size);
				TRI_RECOVER_SIZE_T((*vxs), vxs0, vxs1, vxs2);
				TRI_RECOVER(evx0, evx1, evx2);
				MY_VEC_CHK_RECOVER_POECC(sumX, (void*)evx0, (gsl_vector*)vecX0);
			}

			{ // Recover input/output data Y's backup
			  // and copying back to Y
			  	DBG(printf("[DGEMV_FT3]Recovering Y ... \n"));
				TRI_RECOVER(vecYbak0, vecYbak1, vecYbak2);
				if((unsigned long)vecY_bak!=vecYbak0) vecY_bak=(gsl_vector*)vecYbak0; 
				size_t *vybs = &(((gsl_vector*)vecY_bak)->size);
				TRI_RECOVER_SIZE_T((*vybs), vys0, vys1, vys2); // Y and Y_bak share
				size_t *vybst= &(((gsl_vector*)vecY_bak)->stride);
				TRI_RECOVER_SIZE_T((*vybst),vyst0,vyst1,vyst2);
				TRI_RECOVER(vybkd0,vybkd1,vybkd2);
				if(vecY_bak->data != (double*)vybkd0) vecY_bak->data=(double*)vybkd0;
				TRI_RECOVER(evy0, evy1, evy2);                 // the same dims
				MY_VEC_CHK_RECOVER_POECC(sumY, (void*)evy0, (gsl_vector*)vecYbak0);
				
				TRI_RECOVER(vecY0, vecY1, vecY2);
				if((unsigned long)vecY != vecY0) vecY = (gsl_vector*)vecY0;
				size_t *vys  = &(((gsl_vector*)vecY)->size);
				TRI_RECOVER_SIZE_T((*vys),  vys0, vys1, vys2);
				size_t* vyst = &(((gsl_vector*)vecY)->stride);
				TRI_RECOVER_SIZE_T((*vyst), vyst0,vyst1,vyst2);
				gsl_vector_memcpy(vecY, vecY_bak);
			}
		} 
	}
	DBG(printf("[DGEMV_FT3]Normal call to dgemv.. nonEqualCount=%d\n", nonEqualCount));
	SW3START; 
	SUPERSETJMP("Just before DGEMV");
	if(jmpret != 0) {
		goto kk;
	}
	gsl_blas_dgemv(Trans, alpha, matA, vecX, beta, vecY);
	SW3STOP; 
	my_stopwatch_checkpoint(11); 
	isEqual = Is_GSL_DGEMV_Equal(Trans, alpha, matA, vecX, beta, vecY_bak, vecY);
	my_stopwatch_stop(11); 
	if(isEqual==1) { DBG(printf("[DGEMV_FT3]Result: Equal\n")); }
	else {
		my_stopwatch_checkpoint(6); 
		DBG(printf("[DGEMV_FT3]Result: NOT Equal\n"));
		nonEqualCount = nonEqualCount + 1;
		if(nonEqualCount < NUM_OF_RERUN) {
			DBG(printf("[DGEMV_FT3]Restart calculation.\n"));
			if((nonEqualCount % 4)==4) AdjustIsEqualKnob(1);
			goto kk; 
		}
	}
	my_stopwatch_stop(6); 
	ResetIsEqualKnob();
	DBG(printf("[DGEMV_FT3]Releasing memory for vec_bak's\n"));
	SUPERSETJMP("Just before free");
	if(jmpret == 0) {
		gsl_vector_free(vecY_bak);
		{
			TRI_RECOVER(ema0, ema1, ema2);
			if((unsigned long)ecMatA != ema0) ecMatA = (void*)ema0;
			free(ecMatA);
		}
		{
			TRI_RECOVER(evx0, evx1, evx2);
			if((unsigned long)ecVecX != evx0) ecVecX = (void*)evx0;
			free(ecVecX);
		}
		{
			TRI_RECOVER(evy0, evy1, evy2);
			if((unsigned long)ecVecY != evy0) ecVecY = (void*)evy0;
			free(ecVecY);
		}
		DBG(printf("Released.\n"));
	} else {
		DBG(printf("Oops. Something wrong trying to free the memory. Giving up freeing."));
	}
	MY_REMOVE_SIGSEGV_HANDLER();
}

noinline
void GSL_BLAS_DSYRK_FT3(CBLAS_UPLO_t uplo, CBLAS_TRANSPOSE_t trans, 
	           double alpha, const gsl_matrix* A,
		   double beta, gsl_matrix* C)
{
	/* Protect pointers to matrix A and matrix C */
	unsigned long matA0, matA1, matA2,
	              matC0, matC1, matC2,
		      matCbk0, matCbk1, matCbk2;
	TRIPLICATE(A, matA0, matA1, matA2);
	TRIPLICATE(C, matC0, matC1, matC2);

	MY_SET_SIGSEGV_HANDLER();
	nonEqualCount=0;

	/* Protect size_t of A and C */
	size_t mas10, mas11, mas12, mas20, mas21, mas22;
	size_t mcs10, mcs11, mcs12, mcs20, mcs21, mcs22;
	TRIPLICATE_SIZE_T(A->size1, mas10, mas11, mas12);
	TRIPLICATE_SIZE_T(A->size2, mas20, mas21, mas22);
	TRIPLICATE_SIZE_T(C->size1, mcs10, mcs11, mcs12);
	TRIPLICATE_SIZE_T(C->size2, mcs20, mcs21, mcs22);

	double sumA, sumC;
	sumA=my_sum_matrix(A); sumC=my_sum_matrix(C); 

	/* Protect ptrs to ECC codes. */
	void *ecMatA, *ecMatC;
	unsigned long ema0, ema1, ema2, emc0, emc1, emc2;

	gsl_matrix* C_bak = gsl_matrix_alloc(C->size1, C->size2);
	gsl_matrix_memcpy(C_bak, C);
	TRIPLICATE(C_bak, matCbk0, matCbk1, matCbk2);

	encode((A->data), (A->size1*A->size2), &ecMatA); 
	TRIPLICATE(ecMatA, ema0, ema1, ema2);
	encode((C_bak->data), (C_bak->size1*C_bak->size2), &ecMatC); 
	TRIPLICATE(ecMatC, emc0, emc1, emc2);

	int jmpret, isEqual;
	SUPERSETJMP("[DSYRK_FT3] After encoding\n");
	ResetIsEqualKnob(); // RK's checkers requires larger tolerance, perhaps, perhaps ...
	AdjustIsEqualKnob(3);
	if(jmpret!=0) {
		kk: 
		DBG(printf("[DSYRK_FT3]Recovering matrices from error correction data (jmpret=%d)\n", jmpret));
		size_t *mas1 = &(((gsl_matrix*)A)->size1);
		size_t *mas2 = &(((gsl_matrix*)A)->size2);
		TRI_RECOVER_SIZE_T((*mas1), mas10, mas11, mas12);
		TRI_RECOVER_SIZE_T((*mas2), mas20, mas21, mas22);
		TRI_RECOVER(matA0, matA1, matA2);
		TRI_RECOVER(ema0, ema1, ema2);
		MY_MAT_CHK_RECOVER_POECC(sumA, (void*)ema0, (gsl_matrix*)matA0);

		size_t *mcs1 = &(((gsl_matrix*)C)->size1);
		size_t *mcs2 = &(((gsl_matrix*)C)->size2);
		TRI_RECOVER_SIZE_T((*mcs1), mcs10, mcs11, mcs12);
		TRI_RECOVER_SIZE_T((*mcs2), mcs20, mcs21, mcs22);
		TRI_RECOVER(matCbk0, matCbk1, matCbk2);
		TRI_RECOVER(emc0, emc1, emc2);
		MY_MAT_CHK_RECOVER_POECC(sumC, (void*)emc0, (gsl_matrix*)matCbk0);
		gsl_matrix_memcpy(C, C_bak); 
	}
	DBG(printf("[DSYRK_FT3]Normal call to dsyrk.. nonEqualCount=%d\n", nonEqualCount));
	my_stopwatch_checkpoint(3); 
	SUPERSETJMP("Just before dsyrk");
	if(jmpret != 0) { goto kk; }
	gsl_blas_dsyrk(uplo, trans, alpha, A, beta, C);
	my_stopwatch_stop(3);
	my_stopwatch_checkpoint(11); // Time spent in checks
	isEqual = Is_GSL_DSYRK_Equal2(uplo, trans, alpha, A, beta, C_bak, C); // May use Equal or Equal2
	my_stopwatch_stop(11);
	if(isEqual==1) { DBG(printf("[DSYRK_FT3]Result: Equal\n")); }
	else {
		my_stopwatch_checkpoint(6); 
		DBG(printf("[DSYRK_FT3]Result: NOT Equal\n"));
		nonEqualCount = nonEqualCount + 1;
		if(nonEqualCount < NUM_OF_RERUN) {
			DBG(printf("[DSYRK_FT3]Restart calculation.\n"));
			if((nonEqualCount % 4) == 0) AdjustIsEqualKnob(1);
			goto kk; 
		}
	}
	ResetIsEqualKnob();
	my_stopwatch_stop(6); 
	DBG(printf("[DSYRK_FT3]Releasing memory for C_bak\n"));
	SUPERSETJMP("Just before free");
	if(jmpret == 0) {
		{
			TRI_RECOVER(emc0, emc1, emc2);
			if((unsigned long)ecMatC != emc0) ecMatC = (void*)emc0;
			free(ecMatC);
		}
		{
			TRI_RECOVER(ema0, ema1, ema2);
			if((unsigned long)ecMatA != ema0) ecMatA = (void*)ema0;
			free(ecMatA);
		}
		gsl_matrix_free(C_bak);
		DBG(printf("Released.\n"));
	} else {
		
	}
	MY_REMOVE_SIGSEGV_HANDLER();
}

noinline
void GSL_BLAS_DTRSV_FT3(CBLAS_UPLO_t uplo, CBLAS_TRANSPOSE_t TransA,
			CBLAS_DIAG_t Diag, const gsl_matrix* A, 
			gsl_vector* X)
{
	/* Protect pointers to A and X */
	unsigned long matA_0, matA_1, matA_2,  // A
	              vecX_0, vecX_1, vecX_2,  // X
				  vecXb_0,vecXb_1,vecXb_2; // X_bak
	TRIPLICATE(A, matA_0, matA_1, matA_2);
	TRIPLICATE(X, vecX_0, vecX_1, vecX_2);

	/* Protect the tda, size1, size2, tda of A,
	   stride and size of X */
	size_t mas10, mas11, mas12, mas20, mas21, mas22, matda0, matda1, matda2; // A
	size_t vxst0, vxst1, vxst2, vxs0, vxs1, vxs2;                            // X
	size_t vxbst0, vxbst1, vxbst2, vxbs0, vxbs1, vxbs2;                      // X_bak
	TRIPLICATE_SIZE_T(A->size1, mas10, mas11, mas12);   // A
	TRIPLICATE_SIZE_T(A->size2, mas20, mas21, mas22);   // A
	TRIPLICATE_SIZE_T(A->tda,  matda0, matda1, matda2); // A
	TRIPLICATE_SIZE_T(X->size, vxs0, vxs1, vxs2);       // X
	TRIPLICATE_SIZE_T(X->stride, vxst0, vxst1, vxst2);  // X
	
	MY_SET_SIGSEGV_HANDLER();
	double sumA, sumX;
	DBG(printf("[DTRSV_FT3] 1. Get Sum of Mat and Vec\n"));
	sumA=my_sum_matrix(A); sumX=my_sum_vector(X); 
	void *ecMatA, *ecVecX;
	/* Protect pointers to ECC codes */
	unsigned long ema0, ema1, ema2, evx0, evx1, evx2;
	nonEqualCount=0;
	
	DBG(printf("[DTRSV_FT3] 2. Encoding Input\n"));
	
	encode(A->data, (A->size1*A->size2), &ecMatA); 
	TRIPLICATE(ecMatA, ema0, ema1, ema2); // A

	gsl_vector* X_bak = gsl_vector_alloc(X->size);
	TRIPLICATE(X_bak, vecXb_0, vecXb_1, vecXb_2); // X_bak
	TRIPLICATE(X_bak->size, vxbs0, vxbs1, vxbs2); // X_bak
	TRIPLICATE(X_bak->stride, vxbst0, vxbst1, vxbst2); // X_bak
	encode(X_bak->data, (X_bak->size), &ecVecX); 
	TRIPLICATE(ecVecX, evx0, evx1, evx2); // X
	gsl_vector_memcpy(X_bak, X);
	
	DBG(printf("[DTRSV_FT3] ECC Code Addr: A=%lx, X=%lx\n", (unsigned long)ecMatA,
		(unsigned long)ecVecX));
	int jmpret=0, isEqual;
	SUPERSETJMP("After encoding");
	if(jmpret != 0) {
		kk: 
		DBG(printf("[DTRSV_FT]Recovering matrix and vector from file\n"));
		{ // Recover input data A
			TRI_RECOVER(matA_0, matA_1, matA_2);
			size_t *mas1 = &(((gsl_matrix*)A)->size1),
			       *mas2 = &(((gsl_matrix*)A)->size2), 
				   *matda = &(((gsl_matrix*)A)->tda);
			TRI_RECOVER_SIZE_T((*mas1), mas10, mas11, mas12);
			TRI_RECOVER_SIZE_T((*mas2), mas20, mas21, mas22);
			TRI_RECOVER_SIZE_T((*matda),matda0,matda1,matda2);
			TRI_RECOVER(ema0, ema1, ema2);
			if((unsigned long)ecMatA != ema0) ecMatA=(void*)ema0;
			MY_MAT_CHK_RECOVER_POECC(sumA, ecMatA, (gsl_matrix*)matA_0);
		}
		{ // Recover input data X's backup
			TRI_RECOVER(vecXb_0, vecXb_1, vecXb_2);
			size_t *vxbs = &(X_bak->size), *vxbst = &(X_bak->stride);
			TRI_RECOVER_SIZE_T((*vxbs), vxbs0, vxbs1, vxbs2);
			TRI_RECOVER_SIZE_T((*vxbst),vxbst0,vxbst1,vxbst2);
			TRI_RECOVER(evx0, evx1, evx2);
			if((unsigned long)ecVecX != evx0) ecVecX=(void*)evx0;
			MY_VEC_CHK_RECOVER_POECC(sumX, ecVecX, (gsl_vector*)X_bak);
		}
		{ // Copying back from X_bak to X
			TRI_RECOVER(vecX_0, vecX_1, vecX_2);
			size_t *vxs = &(X->size), *vxst = &(X->stride);
			TRI_RECOVER_SIZE_T((*vxs), vxs0, vxs1, vxs2);
			TRI_RECOVER_SIZE_T((*vxst), vxst0, vxst1, vxst2);
			gsl_vector_memcpy(X, X_bak);
		}
	}
	DBG(printf("[DTRSV_FT]Normal call to dtrsv.. nonEqualCount=%d\n", nonEqualCount));
	SUPERSETJMP("Just before dtrsv");
	if(jmpret != 0) { goto kk; }
	SW3START; 
	int ret = gsl_blas_dtrsv(uplo, TransA, Diag, A, X);
	SW3STOP; 
	if(ret != GSL_SUCCESS) { DBG(printf("[DTRSV_FT3] GSL_ERROR occurred\n")); isEqual=1; /* To force the routine to quit */ } 
	else isEqual = Is_GSL_DTRSV_Equal(uplo, TransA, Diag, A, X_bak, X);
	if(isEqual==1) { DBG(printf("[DTRSV_FT]Result: Equal\n")); }
	else {
		my_stopwatch_checkpoint(6); 
		DBG(printf("[DTRSV_FT]Result: NOT Equal\n"));
		nonEqualCount = nonEqualCount + 1;
		if(nonEqualCount < NUM_OF_RERUN) {
			DBG(printf("[DTRSV_FT]Restart calculation.\n"));
			goto kk; 
		}
	}
	my_stopwatch_stop(6); 
	DBG(printf("[DTRSV_FT]Releasing memory for X_bak\n"));
	SUPERSETJMP("Just before freeing memory");
	if(jmpret == 0) {
		{
			TRI_RECOVER(ema0, ema1, ema2);
			if((unsigned long)ecMatA != ema0) ecMatA = (void*)ema0;
			free(ecMatA);
		}
		{
			TRI_RECOVER(evx0, evx1, evx2);
			if((unsigned long)ecVecX != evx0) ecVecX = (void*)evx0;
			free(ecVecX);
		}
		gsl_vector_free(X_bak);
	}
	DBG(printf("Released.\n"));
	MY_REMOVE_SIGSEGV_HANDLER();
}

noinline
void GSL_LINALG_CHOLESKY_DECOMP_FT3(gsl_matrix* A)
{
	// This one is easier, we only have one matrix to protect
	double sumA; 
	sumA=my_sum_matrix(A);
	
	/* Protect pointer to matrix A */
	unsigned long matA0, matA1, matA2,
	                matAb0, matAb1, matAb2;
	TRIPLICATE(A, matA0, matA1, matA2);
	
	/* Protect A->{size1, size2, tda, data} */
	size_t mas10, mas11, mas12; // A->size1
	size_t mas20, mas21, mas22; // A->size2
	size_t matda0,matda1,matda2;// A->tda
	size_t mad0,  mad1,  mad2;  // A->data
	TRIPLICATE_SIZE_T(A->size1, mas10, mas11, mas12);
	TRIPLICATE_SIZE_T(A->size2, mas20, mas21, mas22);
	TRIPLICATE_SIZE_T(A->tda,   matda0,matda1,matda2);
	TRIPLICATE(A->data,  mad0,  mad1,  mad2);
	
	/* Protect A_bak -> data
	 * A_bak and A have the same dimensions so we save some space */
	size_t mabd0,  mabd1,  mabd2;
	
	void* ecMatA;
	/* Protect pointer to ECC data of matrix A */
	unsigned long ema0, ema1, ema2;
	
	nonEqualCount=0; 
	MY_SET_SIGSEGV_HANDLER(); 
	encode(A->data, (A->size1*A->size2), &ecMatA);
	TRIPLICATE(ecMatA, ema0, ema1, ema2);
	
	/* Back up A */
	gsl_matrix* A_bak = gsl_matrix_alloc(A->size1, A->size2);
	TRIPLICATE(A_bak, matAb0, matAb1, matAb2);
	TRIPLICATE(A_bak->data, mabd0, mabd1, mabd2);
	gsl_matrix_memcpy(A_bak, A);
	
	int jmpret = 0, isEqual; 
	ResetIsEqualKnob();
	SUPERSETJMP("[GSL_LINALG_CHOLESKY_DECOMP_FT3] After encoding\n");
	if(jmpret != 0) {
		kk: 
		
		{
			// Recover A's fields and A_bak's fields
			size_t *mas1 = &(((gsl_matrix*)A)->size1);
			size_t *mabs1= &(((gsl_matrix*)A_bak)->size1);
			size_t *mas2 = &(((gsl_matrix*)A)->size2);
			size_t *mabs2= &(((gsl_matrix*)A_bak)->size2);
			size_t *matda= &(((gsl_matrix*)A)->tda);
			size_t *mabtda=&(((gsl_matrix*)A_bak)->tda);
			TRI_RECOVER_SIZE_T((*mas1), mas10, mas11, mas12);
			TRI_RECOVER_SIZE_T((*mabs1), mas10, mas11, mas12);
			TRI_RECOVER_SIZE_T((*mas2), mas20, mas21, mas22);
			TRI_RECOVER_SIZE_T((*mabs2), mas20, mas21, mas22);
			TRI_RECOVER_SIZE_T((*matda), matda0, matda1, matda2);
			TRI_RECOVER_SIZE_T((*mabtda), matda0, matda1, matda2);
			TRI_RECOVER(mad0, mad1, mad2);
			if(A->data != (double*)mad0) A->data = (double*)mad0;
			TRI_RECOVER(mabd0, mabd1, mabd2);
			if(A_bak->data != (double*)mabd0) A_bak->data=(double*)mabd0;
		}
		{
			// Correct A_bak
			TRI_RECOVER(ema0, ema1, ema2);
			if((unsigned long)ecMatA != ema0) ecMatA=(void*)ema0;
			TRI_RECOVER(matAb0, matAb1, matAb2);
			if((unsigned long)A_bak != matAb0) A_bak = (gsl_matrix*)matAb0;
			MY_MAT_CHK_RECOVER_POECC(sumA, ecMatA, A_bak);
		}
		{
			// Copy back
			TRI_RECOVER(matA0, matA1, matA2);
			if((unsigned long)A != matA0) A = (gsl_matrix*)matA0;
			gsl_matrix_memcpy(A, A_bak);
		}
	}
	DBG(printf("[GSL_LINALG_CHOLESKY_DECOMP_FT3] normal call to cholesky_decomp. nonEqualCount=%d\n", nonEqualCount)); 
	SW3START; 
	SUPERSETJMP("Just before cholesky_decomp");
	if(jmpret != 0) {
		goto kk;
	}
	int ret = gsl_linalg_cholesky_decomp(A); 
	SW3STOP; 
	if(ret != GSL_SUCCESS) { DBG(printf("[GSL_LINALG_CHOLESKY_DECOMP_FT3] GSL_ERROR occurrec\n")); isEqual=1; /* To force the routine to quit */} 
	else isEqual = Is_GSL_linalg_cholesky_decomp_Equal(A_bak, A); 
	if(isEqual==1) { DBG(printf("[GSL_LINALG_CHOLESKY_DECOMP_FT3] Result: Equal\n"));} 
	else {
		my_stopwatch_checkpoint(6); 
		DBG(printf("[GSL_LINALG_CHOLESKY_DECOMP_FT3] Result: NOT Equal\n")); 
		nonEqualCount = nonEqualCount+1; 
		if(nonEqualCount < NUM_OF_RERUN) {
			DBG(printf("[GSL_LINALG_CHOLESKY_DECOMP_FT3] Restart calculation.\n"));
			if((nonEqualCount & 4)==4) AdjustIsEqualKnob(1);
			goto kk; 
		} 
	} 
	ResetIsEqualKnob();
	my_stopwatch_stop(6);
	DBG(printf("[GSL_LINALG_CHOLESKY_DECOMP_FT3] Releasing memory\n"));
	SUPERSETJMP("Just before free");
	if(jmpret == 0)
	{
		TRI_RECOVER(ema0, ema1, ema2);
		if((unsigned long)ecMatA != ema0) ecMatA=(void*)ema0;
		free(ecMatA);
		
		TRI_RECOVER(matAb0, matAb1, matAb2);
		if((unsigned long)A_bak != matAb0) A_bak = (gsl_matrix*)matAb0;
		gsl_matrix_free(A_bak);
	}
	MY_REMOVE_SIGSEGV_HANDLER(); 
}

// Taken from tommy.h
#ifndef POECC_VERSION1
// Returns how much could not be recovered (0 to 1)
noinline
float MY_MAT_CHK_RECOVER_POECC(double sum, void* ecMat, gsl_matrix* mat) {
	MY_SET_SIGSEGV_HANDLER();
	// Protect mat and ecMat
	unsigned long mat0, mat1, mat2;
	unsigned long em0, em1, em2;
	TRIPLICATE(mat, mat0, mat1, mat2);
	TRIPLICATE(ecMat, em0, em1, em2);
	
	float ret = 0;
	my_stopwatch_checkpoint(9);
	double c = my_sum_matrix(mat);
	DBG(printf("[Matrix RECOVER_POECC] Sum=%g vs %g\n", sum, c));
	if(c!=sum) {
		my_stopwatch_checkpoint(8);
		int jmpret=0;
		SUPERSETJMP("MAT_CHK_RECOVER_POECC");
		/* Recover mat and ecMat */
		TRI_RECOVER(mat0, mat1, mat2);
		if((unsigned long)mat != mat0) mat = (gsl_matrix*)mat0;
		TRI_RECOVER(em0, em1, em2);
		unsigned int failed_corr = decode(mat->data, mat->size1*mat->size2, (const double*)em0);
		int n_blks = (int)((mat->size1 * mat->size2)/BLK_LEN/BLK_LEN+1);
		DBG(printf("[Matrix RECOVER_POECC] FAILED corrections = %d out of %d blocks\n", failed_corr, n_blks));
		ret = (float)failed_corr / n_blks;
		my_stopwatch_stop(8);
	} else {
		DBG(printf("[Matrix RECOVER_POECC] Sums equal for ecMat=%lx, mat=%lx\n", (long)ecMat, (long)mat));
	}
	my_stopwatch_stop(9);
	return ret;
}

noinline
float MY_VEC_CHK_RECOVER_POECC(double sum, void* ecVec, gsl_vector* vec) {
	my_stopwatch_checkpoint(9);
	double* ecv1, *ecv2; ecv1=(double*)ecVec; ecv2=(double*)ecVec;
	double c = my_sum_vector(vec);
	if(c!=sum) {
		DBG(printf("[Vector RECOVER_POECC] Sums not equal for ecVec=%lx, vec=%lx.\n", (long)ecVec, (long)vec));
		gsl_vector *vec1, *vec2; vec1=(gsl_vector*)vec; vec2=(gsl_vector*)vec;
		my_stopwatch_checkpoint(8);
		MY_SET_SIGSEGV_HANDLER();
		int jmpret = 0;
		SUPERSETJMP("VEC_CHK_RECOVER_POECC");
		while(vec!=TRIO_DOUBLE(vec,vec1,vec2)) vec=TRIO_DOUBLE(vec,vec1,vec2);
		while(ecVec!=TRIO_DOUBLE(ecVec,ecv1,ecv2)) ecVec=TRIO_DOUBLE(ecVec,ecv1,ecv2);
		decode(vec->data, vec->size, (const double*)ecVec);
		my_stopwatch_stop(8);
	} else {
		DBG(printf("[Vector RECOVER_POECC] Sums equal for ecVec=%lx, vec=%lx.\n", (long)ecVec, (long)vec));
	}
	my_stopwatch_stop(9);
}

noinline double GetIsEqualKnob() {
	return is_equal_knob;
}

noinline double ResetIsEqualKnob() {
	is_equal_knob = FT_TOLERANCE;
	DBG(printf("[ResetIsEqualKnob] Knob reset to %g\n", is_equal_knob));
}

noinline double AdjustIsEqualKnob(int pwr10) {
	int i;
	if(pwr10 > 0) {
		for(i=0; i<pwr10; i++) is_equal_knob *= 10;	
	} else {
		pwr10 = -pwr10;
		for(i=0; i<pwr10; i++) is_equal_knob /= 10;
	}
	DBG(printf("[AdjustIsEqualKnob] Knob adjusted to %g\n", is_equal_knob));
}

#undef PROTECT_IDX_I

#endif
