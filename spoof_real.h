// Used for full-application fault injection.
// Including this file will cause REAL_TRY and REAL_CATCH to become ineffective.
#ifndef __REAL_H__
#define __REAL_H__

#ifdef __cplusplus
extern "C" {
#endif

#define REAL_TRY(LABEL)
#define REAL_CATCH(LABEL)
#define REAL_END(LABEL)

#ifdef __cplusplus
}
#endif

#endif


