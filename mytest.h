// Usage of schifra in here is deprecated!
// 2012-08-09

#ifdef GSL
#include <gsl/gsl_matrix.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif
	void* make_error_correction_data(double* array, int len);
	void error_correct_restore(void* ecdata, double* array);
	int hello();
	#ifdef GSL
	void* make_EC_data_matrix(gsl_matrix* mat);
	void* make_EC_data_vector(gsl_vector* vec);
	gsl_matrix* EC_restore_matrix(void* ecdata);
	gsl_vector* EC_restore_vector(void* ecdata);
	#endif
#ifdef __cplusplus
}
#endif
