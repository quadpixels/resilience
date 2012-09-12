#ifndef _TRIPLICATE_H
#define _TRIPLICATE_H

#include <gsl/gsl_matrix.h>

#define noinline __attribute__((noinline))
noinline void trick_me_jr(int);
noinline unsigned long trick_me_ptr(unsigned long, int);
noinline size_t trick_me_size_t(size_t, size_t);
noinline void trick_blackhole(size_t);
noinline unsigned int GetMatrixChecksum(const gsl_matrix*);

// The same as above. We can't change the field of a const gsl_matrix.
#define TRIPLICATE_SIZE_T(in, s0, s1, s2) { \
	s0 = trick_me_size_t(in, 0); \
	s1 = trick_me_size_t(in, 123); \
	s2 = trick_me_size_t(in, 456); \
}

// This time we must correct the field of the const gsl_matrix.
#define TRMSG DBG(printf("[TRI_RECOVER_SIZE_T] Corrected 1 size_t"))
#define TRI_RECOVER_SIZE_T(in, s0, s1, s2) { \
	if((s0 != s1-123) && (s1-123 == s2-456)) { s0 = s1-123; TRMSG; } \
	if((s0 != s1-123) && (s0 == s2-456)) { s1 = s0 + 123; TRMSG; } \
	if((s1-123 != s2-456) && (s0 == s1-123)) { s2 = s0 + 456; TRMSG; } \
	if(in != s0) in = s0; \
}
	
// in = whatever, r1, r2 are unsigned long's
// We can't change the const gsl_matrix*, so we need make three copies,
// and change these copies.
#define TRIPLICATE(in, r0, r1, r2) \
	r0 = trick_me_ptr((unsigned long)in, 0); \
	r1 = trick_me_ptr((unsigned long)in, 123); \
	r2 = trick_me_ptr((unsigned long)in, 456);

#define TRI_RECOVER(r0, r1, r2) { \
	if((r0 != r1-123) && (r1-123 == r2-456)) { r0 = r1-123; DBG(printf("[TRI_RECOVER] Corrected 1 element\n"));} \
	else if((r0 != r1-123) && (r0 == r2-456)) { r1 = r0+123; DBG(printf("[TRI_RECOVER] Corrected 1 element\n"));} \
	else if((r0 != r2-456) && (r0 == r1-123)) { r2 = r0+456; DBG(printf("[TRI_RECOVER] Corrected 1 element\n"));} \
}

#endif
