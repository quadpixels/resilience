#include "tommy.h"
#include <sys/time.h>

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

#define SWSZ 20 // SWSZ = StopWatch SiZe
		// Stopwatch [9] is reserved for counting time spent in checking + recovery
		// Stopwatch [8] is reserved for counting time spent in recovery
		// Stopwatch [7] is reserved for counting time spent in preparation (checking ECC codes or backing up to file)
		// Stopwatch [6] is reserved for counting re-runs
		// Stopwatch [5] is reserved for sigsegv handlers
		// Stopwatch [4] is reserved for result checkers
		// Stopwatch [3] is reserved for counting time spent in calculation


char _swMessage[][80] = {
	"Stopwatch 0",
	"Stopwatch 1",
	"Stopwatch 2",
	"Stopwatch 3 (calculation time)",
	"Stopwatch 4 (result chkr time)",
	"Stopwatch 5 (sigsegv handler time)",
	"Stopwatch 6 (re-run time)",
	"Stopwatch 7 (preparation time)",
	"Stopwatch 8 (input data recovery time)",
	"Stopwatch 9 (input data checking and recovery time) "
	"Stopwatch 10",
	"Stopwatch 11",
	"Stopwatch 12",
	"Stopwatch 13",
	"Stopwatch 14",
	"Stopwatch 15",
	"Stopwatch 16",
	"Stopwatch 17",
	"Stopwatch 18",
	"Stopwatch 19",
	};

struct timeval g_tvTick[SWSZ], g_tvTock[SWSZ];
struct timeval* g_ptv[SWSZ];
unsigned long g_time[SWSZ];

//#define PROFILE_DISABLED
#ifndef PROFILE_DISABLED
noinline
void my_stopwatch_checkpoint(int id) {
	STOPWATCH_INCREMENT(id);
	DBG(printf("[my_stopwatch_checkpoint #%d] Tick:[%lu,%lu], Tock:[%lu,%lu]\n",
		id,
		(unsigned long)g_tvTick[id].tv_sec, (unsigned long)g_tvTick[id].tv_usec,
		(unsigned long)g_tvTock[id].tv_sec, (unsigned long)g_tvTock[id].tv_usec));
}

noinline
void my_stopwatch_stop(int id) {
	STOPWATCH_INCREMENT(id);
	ptv_late->tv_usec = 0; ptv_late->tv_sec = 0;
	ptv_early->tv_usec = 0; ptv_early->tv_sec = 0;
	DBG(printf("[my_stopwatch_stop #%d] Tick:[%lu,%lu], Tock:[%lu,%lu]\n",
		id,
		(unsigned long)g_tvTick[id].tv_sec, (unsigned long)g_tvTick[id].tv_usec,
		(unsigned long)g_tvTock[id].tv_sec, (unsigned long)g_tvTock[id].tv_usec));

}

noinline unsigned long my_stopwatch_get(int id) {
	return g_time[id];
}

noinline void my_stopwatch_show(int id, char* extra_msg, char writeToFile) {
	char k[100];
	sprintf(k, "[my_stopwatch_show #%d] %lu microseconds (%s|%s).\n", id, my_stopwatch_get(id), _swMessage[id], extra_msg);
	printf("%s", k);
	if(writeToFile==1) {
		FILE* f = fopen("My_Stopwatch_Log", "a");
		fprintf(f, "%s", k);
		fclose(f);
	}
	return;
}
#else
noinline void my_stopwatch_checkpoint(int id) {}
noinline void my_stopwatch_stop(int id) {}
noinline unsigned long my_stopwatch_get(int id) {}
noinline void my_stopwatch_show(int id, char* extra_msg, char writeToFile) {}
#endif


#define SW3START my_stopwatch_checkpoint(3)
#define SW3STOP my_stopwatch_stop(3)

