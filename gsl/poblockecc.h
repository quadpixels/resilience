// This library is deprecated !
// 08-09-2012
#ifdef GSL
#include <gsl/gsl_matrix.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif
	void* make_PoEC_correction_data(double* array, int len);
	void PoEC_correct_restore(void* ecdata, double* array);
	#ifdef GSL
	void* make_PoEC_data_matrix(const gsl_matrix* mat);
	void* make_PoEC_data_vector(const gsl_vector* vec);
	void PoEC_restore_matrix(void* ecdata, gsl_matrix* mat); // This guy is different --- in-place recovery!
	void PoEC_restore_vector(void* ecdata, gsl_vector* vec); // This guy is different --- in-place recovery!
	#endif
#ifdef __cplusplus
}
#endif
