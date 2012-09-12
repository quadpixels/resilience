#include "mytest.h"
#include "tommy.h"
#include "mytest_deps.h"
#include <gsl/gsl_matrix.h>
#define forceinline __inline__ __attribute__((always_inline))

int LEN;
double* data;
#define MYMASK 0x00
/* 04-03-2012: Tweak these two parameters to balance between protection and performance. */
#define CODE_LENGTH 255
#define FEC_LENGTH 8
typedef std::vector<schifra::reed_solomon::block<255, 8 > > MYCONTAINER;

int AtoI(char* input) {
   int ret; std::stringstream ssout(input);
   ssout>>ret;
   return ret;
}

//template <std::size_t code_length, std::size_t fec_length>
//MYCONTAINER double_array_encode(double* d, int length) {
MYCONTAINER double_array_decode(double* d, int length) {
	std::cout << "[encode]";
	printf("Length: %d\n", length);
	std::size_t field_descriptor = 8;
	std::size_t generator_poly_index = 120;
	std::size_t generator_poly_root_count = FEC_LENGTH;
	
	schifra::galois::field field(field_descriptor,
				     schifra::galois::primitive_polynomial_size06,
				     schifra::galois::primitive_polynomial06);
	schifra::galois::field_polynomial generator_polynomial(field);
	schifra::sequential_root_generator_polynomial_creator(field,
							      generator_poly_index,
							      generator_poly_root_count,
							      generator_polynomial);
	schifra::reed_solomon::encoder<CODE_LENGTH, FEC_LENGTH> encoder(field, generator_polynomial);
	MYCONTAINER ret;
	printf("1");
	int i, j, k=0;
	int nElem = (int)((CODE_LENGTH - FEC_LENGTH)/sizeof(double));
	std::stringstream sstream;
	int elemsDone = 0;

	double sz2 = (double)length;
	for(j=0; j<(int)sizeof(double); j++) {
		unsigned char c = *((unsigned char*)&sz2 + j);
		c = c ^ MYMASK;
		sstream << c;
	}
	elemsDone = elemsDone + 1;
	printf("2");
	for(i=0; i<length; i++) {
		elemsDone = elemsDone + 1;
		double dd = d[i];
//		printf(" %g=", dd);
		for(j=0; j<(int)sizeof(double); j++) {
			unsigned char c = *(unsigned char*)((unsigned char*)&dd + j);
			c = c ^ MYMASK;
//			printf("%d ", c);
			sstream << c;
		}
		if((elemsDone>0 && elemsDone%nElem==0) || (i==length-1)) {
//			printf("[encode]elemsDone=%d\n", (unsigned int)elemsDone);
			sstream.flush();
			std::string s = sstream.str();
//			std::cout<<"Encoding block #" << k <<" (length=" << s.length() << ")";
			k = k + 1;
//			schifra::reed_solomon::block<code_length, fec_length> block;
			schifra::reed_solomon::block<CODE_LENGTH, FEC_LENGTH> block;
			if(!encoder.encode(s, block)) {
				std::cout << "ERROR: Critical encoding failure!!!" << std::endl;
			} else {
//				std::cout << "ok." << std::endl;
//				schifra::reed_solomon::block<code_length, fec_length> *block2 = new schifra::reed_solomon::block<code_length, fec_length>(block);

		//04-08-2012
				schifra::reed_solomon::block<CODE_LENGTH, FEC_LENGTH> *block2 = new schifra::reed_solomon::block<CODE_LENGTH, FEC_LENGTH>(block);
//				MYCONTAINER *block2 = new schifra::reed_solomon::block<CODE_LENGTH, FEC_LENGTH>(block);
				ret.push_back(*block2);
				sstream.str("");
			}
		}
	}
	printf("[encode]len=%d, elemsDone=%d\n", (unsigned int)length, (unsigned int)elemsDone);
	return ret;
}

//template <std::size_t code_length, std::size_t fec_length>
//void double_array_decode(std::list<schifra::reed_solomon::block<code_length, fec_length> >list, double* d) {
double* double_array_decode(MYCONTAINER list, double* d) {
	std::cout << "[decode]# of blocks: " << list.size() << std::endl;
	std::size_t field_descriptor = 8;
	std::size_t generator_poly_index = 120;
	std::size_t generator_poly_root_count = FEC_LENGTH;
	
	schifra::galois::field field(field_descriptor,
				     schifra::galois::primitive_polynomial_size06,
				     schifra::galois::primitive_polynomial06);
	schifra::galois::field_polynomial generator_polynomial(field);
	schifra::sequential_root_generator_polynomial_creator(field,
							      generator_poly_index,
							      generator_poly_root_count,
							      generator_polynomial);
	schifra::reed_solomon::decoder<CODE_LENGTH, FEC_LENGTH> decoder(field, generator_poly_index);
	/* Iterate through the list and decode */

	MYCONTAINER::iterator itr_list; // Can't use templates "code_length" and "fec_length" here.....
	int blksDecoded = 0, message_len=0, pArray=0, arrayLen=0;
	bool isNewArray = false;
	for(itr_list = list.begin(); itr_list != list.end(); itr_list ++) {
//		schifra::reed_solomon::block<code_length, fec_length> blk = *itr_list;
		schifra::reed_solomon::block<CODE_LENGTH, FEC_LENGTH> blk = *itr_list;
		if(!decoder.decode(blk)) {
//			printf("Decode Error.\n");
		} else {
			int ptr=0;
//			printf("Decode OK!\n");
			std::string message = std::string((CODE_LENGTH - FEC_LENGTH), static_cast<unsigned char>(0xC0));
			blk.data_to_string(message);
			const char* c=message.c_str();
//			char* c = strdup(message.c_str());
			message_len = message.length();
//			printf("[decode] String Length: %d\n", message_len);
			if(blksDecoded==0) {
				double k;
				int i; for(i=0; i<(int)(sizeof(double)); i++) {
					unsigned char p = *((unsigned char*)c + ptr);
//					printf("%d ", p);
					p = p ^ MYMASK;
					*(unsigned char*)((unsigned char*)&k+i)=p;
					ptr++;
				}
				printf("[decode]Array Length: %d\n", (int)k);
				arrayLen = (int)k;
				if(d==NULL) { printf("[decode] allocating result. %d elems.\n", arrayLen);
						d=(double*)malloc(sizeof(double)*arrayLen); isNewArray = true;}
			}
			int nElemsRounded = sizeof(double) * (int)((CODE_LENGTH - FEC_LENGTH)/sizeof(double));
			while(ptr < nElemsRounded) {
				double k;
				int i; for(i=0; i<(int)(sizeof(double)); i++) {
					unsigned char p = *((unsigned char*)c + ptr);
//					printf("%d ", p);
					p = p ^ MYMASK;
					*(unsigned char*)((unsigned char*)&k+i)=p;
					ptr++;
				}
//				printf("=%g ", k);
				if(pArray <= arrayLen-1) {
					if(isNewArray==false && d[pArray] != k) { printf("[decode] Detected difference: d[%d]=%g (should be %g)\n}", pArray, d[pArray], k); }
					d[pArray]=k;
					pArray++;
				}
			}
			blksDecoded++;
		}
	}
	return d;
}

template <typename T>
std::string To_String(T* array, int len) {
   printf("elem size: %d, elem count: %d\n", (int)sizeof(T), len);
   std::stringstream sstr;
   int i, j;
   for(i=0; i<len; i++) {
      for(j=0; j<sizeof(T); j++)
         sstr<<static_cast<unsigned char>(*(char*)(&array[i]+j));
   }
   std::string ret = sstr.str();
   std::cout << "Contents: " << ret << std::endl;
   return ret;
}

int hello() 
//int main(int argc, char** argv)
{
   double k[200];
   int x=0; for(x=0; x<200; x++) k[x]=x;
//   MYCONTAINER list1 = double_array_encode(k, 200);
   MYCONTAINER list1 = double_array_decode(k, 200);
   double_array_decode(list1, k);
   return 0;
}

void* make_error_correction_data(double* array, int len) {
	MYCONTAINER list1 = double_array_decode(array, len);
	MYCONTAINER* ret = new MYCONTAINER(list1);
	printf("Encoded..\n");
	printf("Size of block: %ld\n", sizeof(schifra::reed_solomon::block<CODE_LENGTH, FEC_LENGTH>));
	return (void*)ret;
}

void error_correct_restore(void* ecdata, double* array) {
	MYCONTAINER list = *(MYCONTAINER*)ecdata;
	std::cout << "[decode]# of blocks: " << list.size() << std::endl;
//	double_array_decode<CODE_LENGTH, FEC_LENGTH>(list, array);
	double_array_decode(list, array);
}

void* make_EC_data_matrix(gsl_matrix* mat) {
	my_stopwatch_checkpoint(7);
	int len = mat->size1 * mat->size2;
	printf("[make_EC_data_matrix]GSL_Matrix encoded.\n");
	double* array = (double*)malloc(sizeof(double) * (len + 2));
	array[0] = (double)(mat->size1);
	array[1] = (double)(mat->size2);
	int i,j,k; for(i=0; i<(int)(mat->size1); i++) {
		for(j=0; j<(int)(mat->size2); j++) {
			k = i * mat->size2 + j + 2;
			array[k] = gsl_matrix_get(mat, i, j);
		}
	}
	void* ecdata = make_error_correction_data(array, (len+2));
	free(array);
	my_stopwatch_stop(7);
	return ecdata;
}

void* make_EC_data_vector(gsl_vector* vec) {
	my_stopwatch_checkpoint(7);
	int len = vec->size;
	printf("[make_EC_data_vector]GSL_Vector encoded.\n");
	double* array = (double*)malloc(sizeof(double) * (len + 1));
	array[0] = (double)len;
	int i, j; for(i=0; i<len; i++) {
		j = i+1;
		array[j] = gsl_vector_get(vec, i);
	}
	void* ecdata = make_error_correction_data(array, (len+1));
	free(array);
	my_stopwatch_stop(7);
	return ecdata;
}

gsl_matrix* EC_restore_matrix(void* ecdata) {
	MYCONTAINER list = *(MYCONTAINER*)ecdata;
	printf("[EC_Restore_Matrix] ecdata=%x\n", *(int*)ecdata);
	double* k = double_array_decode(list, NULL);
	int size1 = k[0], size2 = k[1], i, j, l;
	printf("[EC_Restore_Matrix] size=%d x %d\n", size1, size2);
	gsl_matrix* ret = gsl_matrix_alloc(size1, size2);
	for(i=0; i<size1; i++) {
		for(j=0; j<size2; j++) {
			l = i*size2 + j + 2;
			gsl_matrix_set(ret, i, j, k[l]);
		}}
	printf("[Matrix EC Restore] Matrix size: %d x %d\n", size1, size2);
	free(k);
	return ret;
}

gsl_vector* EC_restore_vector(void* ecdata) {
	MYCONTAINER list = *(MYCONTAINER*)ecdata;
	printf("[EC_Restore_Vector] ecdata=%x\n", *(int*)ecdata);
	double* k = double_array_decode(list, NULL);
	int size = k[0], i, j;
	printf("[EC_Restore_Vector] size=%d\n", size);
	gsl_vector* ret = gsl_vector_alloc(size);
	for(i=0; i<size; i++) {
		j = i + 1;
		gsl_vector_set(ret, i, k[j]);
	}
	printf("[Vector EC Restore] Vector size: %d\n", size);
	free(k);
	return ret;
}
