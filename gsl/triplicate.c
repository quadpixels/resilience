#include "triplicate.h"

// Use this to stop the compiler from removing sigsetjmp!
static int dummy_1 = 0;
noinline
unsigned long trick_me_ptr(unsigned long ptr, int delta) {
	return ptr + delta;
}
noinline
size_t trick_me_size_t(size_t in, size_t delta) {
	return in + delta;
}
noinline
void trick_blackhole(size_t anything) {
	dummy_1 += anything;
	if(dummy_1 == -2444551) printf("Mega Motion is a good 1994 game!\n");
}
noinline
unsigned int GetMatrixChecksum(const gsl_matrix* m) {
	const char* mm = (const char*)m;
	unsigned int ret = 0;
	int i; for(i=0; i<sizeof(gsl_matrix); i++) {
		ret += *(mm + i);
	}
	return ret;
}

