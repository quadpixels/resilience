#include "tommy.h"
#include "tommy_config.h"

// Used for sigmentation fault handler and long jmp.

#ifdef __cplusplus
extern "C" {
#endif

jmp_buf buf, buf_1, buf_2;

int _count = 0;
int _countmax = 20;
int num_total_retries = 0;
int num_retries_mm = 0;
int num_retries_mv = 0;
int num_retries_rk = 0;
int num_retries_sv = 0;
int num_retries_cd = 0;

noinline void trick_me_jr(int jmpret) {
	if(jmpret == 999) printf("Jack is in the box!\n");
}

noinline
unsigned int GetJmpBufChecksum(jmp_buf* buf) {
	unsigned int ret = 0; int i=0;
	for(i=0; i<sizeof(jmp_buf); i++) {
		ret += *((char*)(buf) + i);
	}
	return ret;
}

noinline
void MY_SET_SIGSEGV_HANDLER() {
	struct sigaction sa;
	memset(&sa, 0, sizeof(struct sigaction));
	sa.sa_handler = my_action;
	sigemptyset(&sa.sa_mask);
	sigaction(SIGSEGV, &sa, NULL);
	sigaction(SIGBUS,  &sa, NULL);
	sigaction(SIGILL,  &sa, NULL);
}

noinline
void MY_REMOVE_SIGSEGV_HANDLER() {
	struct sigaction sa;
	memset(&sa, 0, sizeof(struct sigaction));
	sa.sa_handler = SIG_DFL;
	sigemptyset(&sa.sa_mask);
	sigaction(SIGSEGV, &sa, NULL);
	sigaction(SIGBUS,  &sa, NULL);
}

noinline
void my_action(int sig) {
	my_stopwatch_checkpoint(5);
	_count = _count + 1;
	printf(" >> Caught SIGSEGV signal (%d out of %d allowed)\n", _count, NUM_OF_SIGSEGV);
	time_t rawtime; struct tm *timeinfo; time(&rawtime); timeinfo = localtime(&rawtime);
	printf("    Time: %s", asctime(timeinfo)); // Structure tm is statically allocated.
	
	/* Backtrace stuff */
	void* array[10];
	char** strings;
	size_t size;
	size = backtrace(array, 10);
	printf(" >> Stack contents (level=%ld):\n", size);
	backtrace_symbols_fd(array, size, 2);
#ifdef UNWIND
	printf(" >> Stack contents by libunwind\n");
	do_backtrace();
	printf(" << End of stack contents <<<<<<\n");
#endif
	if(_count >= NUM_OF_SIGSEGV) { 
		printf("[SIGSEGV_HANDLER] Before aborting, re-calculation has taken place %d times.\n", num_total_retries);
		abort(); 
	}
	else { 
		my_stopwatch_stop(5);
		unsigned long sum1 = GetJmpBufChecksum(&buf), sum2 = GetJmpBufChecksum(&buf_1),
		sum3 = GetJmpBufChecksum(&buf_2);
		DBG(printf("Byzantine tolerance of buf. sum: %lu, %lu, %lu.\n", sum1, sum2, sum3));
		if(sum1 != sum2 && sum2 == sum3) memcpy(&buf, &buf_1, sizeof(jmp_buf));
		else if(sum2 != sum1 && sum1 == sum3) memcpy(&buf_1, &buf, sizeof(jmp_buf));
		else if(sum3 != sum2 && sum2 == sum1) memcpy(&buf_2, &buf, sizeof(jmp_buf));
		DBG(printf("Now, sum(buf) = %u\n", GetJmpBufChecksum(&buf)));
		printf("\n");
		siglongjmp(buf, 1); 
	}
	printf(" >> Execution not expected here.\n");
	abort();
}

#ifdef UNWIND
#include "libunwind.h"
void do_backtrace() {
	unw_cursor_t cursor;
	unw_context_t context;
	int ret;
	ret = unw_getcontext(&context);
	if(ret != 0) { printf("[do_backtrace] unw_getcontext returned %d.\n", ret); }
	ret = unw_init_local(&cursor, &context);
	if(ret != 0) { printf("[do_backtrace] unw_init_local returned %d.\n", ret); }

	do {
		unw_word_t offset, pc;
		char fname[64];
		unw_get_reg(&cursor, UNW_REG_IP, &pc);

		fname[0]='\0';
		(void) unw_get_proc_name(&cursor, fname, sizeof(fname), &offset);
		printf("%p : (%s+0x%x) [%p]\n", (void*)pc, fname, (unsigned int)offset, (void*)pc);
	}
	while(unw_step(&cursor) > 0);
}
#endif

void MY_SIGSEGV_HANDLER_SUMMARY() {
	printf("[MySigsegvHandler] # of faults handled: %d\n", _count);
	printf("[FT3] # of retries - MM: %d\n", num_retries_mm);
	printf("[FT3] # of retries - MV: %d\n", num_retries_mv);
	printf("[FT3] # of retries - RK: %d\n", num_retries_rk);
	printf("[FT3] # of retries - SV: %d\n", num_retries_sv);
	printf("[FT3] # of retries - Cholesky Decomp: %d\n", num_retries_cd);
}

#ifdef __cplusplus
}
#endif

