// Taken from RELAX framework.
// Used in conjunction with spoof_real.h.
#ifndef __REAL_H__
#define __REAL_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
void __POKE__ProfilingBegin(int);
void __POKE__ProfilingEnd(int);

#define ESCAPE_DECLARE(TYPE) \
extern TYPE llvm_real_escape_declare_##TYPE(unsigned id, TYPE val)
ESCAPE_DECLARE(float);
ESCAPE_DECLARE(double);
ESCAPE_DECLARE(uint32_t);
ESCAPE_DECLARE(uint16_t);
ESCAPE_DECLARE(uint8_t);
ESCAPE_DECLARE(int32_t);
ESCAPE_DECLARE(int16_t);
ESCAPE_DECLARE(int8_t);
#undef ESCAPE_DECLARE

#define ESCAPE_ASSIGN_ATOMIC_LAST(TYPE) \
extern TYPE llvm_real_escape_assign_atomic_last_##TYPE(unsigned id, TYPE val)
ESCAPE_ASSIGN_ATOMIC_LAST(float);
ESCAPE_ASSIGN_ATOMIC_LAST(double);
ESCAPE_ASSIGN_ATOMIC_LAST(uint32_t);
ESCAPE_ASSIGN_ATOMIC_LAST(uint16_t);
ESCAPE_ASSIGN_ATOMIC_LAST(uint8_t);
ESCAPE_ASSIGN_ATOMIC_LAST(int32_t);
ESCAPE_ASSIGN_ATOMIC_LAST(int16_t);
ESCAPE_ASSIGN_ATOMIC_LAST(int8_t);
#undef ESCAPE_ASSIGN_ATOMIC_LAST

#define ESCAPE_ASSIGN_ATOMIC(TYPE) \
extern TYPE llvm_real_escape_assign_atomic_##TYPE(unsigned id, TYPE val)
ESCAPE_ASSIGN_ATOMIC(float);
ESCAPE_ASSIGN_ATOMIC(double);
ESCAPE_ASSIGN_ATOMIC(uint32_t);
ESCAPE_ASSIGN_ATOMIC(uint16_t);
ESCAPE_ASSIGN_ATOMIC(uint8_t);
ESCAPE_ASSIGN_ATOMIC(int32_t);
ESCAPE_ASSIGN_ATOMIC(int16_t);
ESCAPE_ASSIGN_ATOMIC(int8_t);
#undef ESCAPE_ASSIGN_ATOMIC

#define ESCAPE_ASSIGN(TYPE) \
extern TYPE llvm_real_escape_assign_##TYPE(unsigned id, TYPE val)
ESCAPE_ASSIGN(float);
ESCAPE_ASSIGN(double);
ESCAPE_ASSIGN(uint32_t);
ESCAPE_ASSIGN(uint16_t);
ESCAPE_ASSIGN(uint8_t);
ESCAPE_ASSIGN(int32_t);
ESCAPE_ASSIGN(int16_t);
ESCAPE_ASSIGN(int8_t);
#undef ESCAPE_ASSIGN

extern void llvm_real_try(uint32_t);
extern void llvm_real_split(uint32_t);
extern void llvm_real_catch(uint32_t);
extern void llvm_real_end(uint32_t);

extern uint32_t dummy;
#define REAL_TRY(LABEL) \
llvm_real_split(LABEL); \
if (dummy) goto jump##LABEL; \
goto real_catch##LABEL; \
jump##LABEL: \
llvm_real_split(LABEL); \
if (dummy) goto real_try##LABEL; \
goto real_end##LABEL; \
real_try##LABEL: llvm_real_try(LABEL);

#define REAL_CATCH(LABEL) \
goto real_end##LABEL; \
real_catch##LABEL: llvm_real_catch(LABEL);

#define REAL_END(LABEL) \
real_end##LABEL: llvm_real_end(LABEL);

#define REAL_ESCAPE_DECLARE(LABEL, TYPE, ID, VAL) \
TYPE ID = llvm_real_escape_declare_##TYPE(LABEL, VAL); 

#define REAL_ESCAPE_ASSIGN(LABEL, TYPE, ID, VAL) \
ID = llvm_real_escape_assign_##TYPE(LABEL, VAL);

#define REAL_ESCAPE_ASSIGN_ATOMIC(LABEL, TYPE, ID, VAL) \
ID = llvm_real_escape_assign_atomic_##TYPE(LABEL, VAL);

#define REAL_ESCAPE_ASSIGN_ATOMIC_LAST(LABEL, TYPE, ID, VAL) \
ID = llvm_real_escape_assign_atomic_last_##TYPE(LABEL, VAL);

#ifdef __cplusplus
}
#endif

#endif


