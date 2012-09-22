// Configuration header file for soft error analysis //

// 1. Enable debug verbose output?
//#define DEBUG
#ifdef DEBUG
#define DBG(call) { call; }
#else   
#define DBG(call) {}
#endif

// 1.5 Log all calls to (FT) GSL routines ?
#define LOGCALLS
#ifdef LOGCALLS
#define LOGCALL(call) { printf(">>>---[PID=%d]", getpid()); call; }
#else
#define LOGCALL(call) {}
#endif

// 2. Number of maximum consecutive re-calculations
#define NUM_OF_RERUN 35	// Up to how many re-runs ? If the number is exceeded the program would return the last result

// 3. Number of maximum segmentation faults to handle
#define NUM_OF_SIGSEGV 1500 // Up to how many SIGSEGV's ? If the number is exceeded the program would be aborted.

// 4. Algorithmic fault checker error detection threshold
#ifdef DTRSV_EXPERIMENT
#define FT_TOLERANCE 0.01
#else
#define FT_TOLERANCE       0.0001
#endif

// 5. Disable profiling ?
//#define PROFILE_DISABLED

// 6. Enable input checksum for MV ?
//#define ENABLE_POECC_MV
