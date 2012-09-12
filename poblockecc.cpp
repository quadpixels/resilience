// Usage of this library is deprecated !!
// 2012-08-09
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include "poblockecc.h"
#include "tommy.h"

#define BLKSIZE 3

//#define DEBUG

#ifdef DEBUG
#define DBG(call) { call; }
#else
#define DBG(call) {}
#endif

template<typename DATA, unsigned int blocklen>
class PoBlocksum {
public:
	PoBlocksum() {
		rowsums=(DATA*)malloc(sizeof(DATA) * blocklen);
		colsums=(DATA*)malloc(sizeof(DATA) * blocklen);
	}
	~PoBlocksum() {
		free(rowsums);
		free(colsums);
	}

	/* Push blocklen*blocklen data entries */
	void PushBlock(DATA* block) {
		unsigned int i, j;
		/* Compute row sums and column sums */
		DATA rsum, csum, allsum=0;
		for(i=0; i<blocklen; i++) {
			rsum=0; csum=0;
			for(j=0; j<blocklen; j++) {
				rsum = rsum + block[i*blocklen + j];
				csum = csum + block[j*blocklen + i];
			}
			allsum = allsum + rsum;
			rowsums[i] = rsum;
			colsums[i] = csum;
//			DBG(printf("Row sum #%d = %g\n", i, rsum));
//			DBG(printf("Col sum #%d = %g\n", i, csum));
		}
	}

	/* Check against stored rowsums and colsums and see if everything is correct */
	/* The positions of the corrected elements are pushed into pos. 
	 * The corrected numbers are pushed into naosu. 
	 * Returns: Number of elements corrected
	 * */
	unsigned int CheckAndCorrect(DATA* newblk, std::vector<int>* pos, std::vector<DATA>* naosu) {
		unsigned int i, j, ret=0;
		DATA rsum, csum, allsum=0;
		std::vector<unsigned int> errRows;
		std::vector<unsigned int> errCols;
		for(i=0; i<blocklen; i++) {
			rsum=0; csum=0;
			for(j=0; j<blocklen; j++) {
				rsum = rsum + newblk[i*blocklen + j];
				csum = csum + newblk[j*blocklen + i];
			}
			allsum = allsum + rsum;
			if(rsum != rowsums[i]) {
				DBG(printf("[PoBlocksum::CheckAndCorrect] Row %d is erroneous (sum=%g, should be %g)\n", i, rsum, rowsums[i]));
				errRows.push_back(i);
			}
			if(csum != colsums[i]) {
				DBG(printf("[PoBlocksum::CheckAndCorrect] Col %d is erroneous (sum=%g, should be %g)\n", i, csum, colsums[i]));
				errCols.push_back(i);
			}
		}
		if(errRows.size() > 0 || errCols.size() > 0) DBG(printf("[PoBlocksum::CheckAndCorrect] Different rolsums=%d, colsums=%d\n", (int)(errRows.size()), (int)(errCols.size())));
		if(errRows.size() > 1 && errCols.size() > 1) {
			DBG(printf(" There are more than 1 wrong rowsums and more than 1 colsums -- cannot recover.\n"));
		} else {
			if(errRows.size()==1) {
				int rowid=errRows[0];
				for(i=0; i<errCols.size(); i++) {
					int colid=errCols[i];
					DATA errElem = newblk[rowid*blocklen + colid];
					DATA elem = colsums[colid]; // subtract the colsum of the erroneous column to all elements except the one on the erroneous row
					for(j=0; j<blocklen; j++) {
						if(j!=rowid) elem = elem - newblk[j*blocklen + colid];
					}
					DBG(printf(" element [%d, %d] has been corrected from %g to %g\n", rowid, colid, errElem, elem));
					pos->push_back(rowid*blocklen + colid);
					naosu->push_back(elem);
					ret++;
				}
			}
			else if(errCols.size()==1) {
				int colid=errCols[0];
				for(i=0; i<errRows.size(); i++) {
					int rowid=errRows[i];
					DATA errElem = newblk[rowid*blocklen + colid];
					DATA elem = rowsums[rowid];
					for(j=0; j<blocklen; j++) {
						if(j!=colid) elem = elem - newblk[rowid*blocklen + j];
					}
					DBG(printf(" element [%d, %d] has been corrected from %g to %g\n", rowid, colid, errElem, elem));
					pos->push_back(rowid*blocklen + colid);
					naosu->push_back(elem);
					ret++;
				}
			}
		}
		return ret;
	}
private:
	DATA* rowsums;
	DATA* colsums;
	DATA allsum;
};

template<typename DATA, unsigned int blocklen>
class PoBlock {
public:
	enum State {
		POBLOCK_MAKING_ECCCODE,
		POBLOCK_READY_TO_CORRECT,
		POBLOCK_CORRECTING
	};
	PoBlock(unsigned int sz) {
		arraysize = sz;
		nPushed = 0; pC=0;
		cache = (DATA*)malloc(sizeof(DATA) * blocklen * blocklen);
		DBG(printf("PoBlock for an array of size %d\n", sz));
		state = POBLOCK_MAKING_ECCCODE;
	}
	void PushDataForBackup(DATA d) {
		if(nPushed == arraysize) {
			DBG(printf("[PoBlock::PushData] This ECC word is already filled.\n"));
			return;
		}
		cache[pC] = d;
		pC++; nPushed++;
		if((pC == blocklen * blocklen) ||
		   (nPushed == arraysize)) {
		   	if(pC < blocklen*blocklen) {
				DBG(printf("."));
				unsigned int i; for(i=pC; i<blocklen*blocklen; i++) cache[i]=0;
			}
		   	PoBlocksum<DATA, blocklen>* newblocksum = new PoBlocksum<DATA, blocklen>();
			newblocksum->PushBlock(cache);
			blocksums.push_back(newblocksum);
			pC=0;
		   }
		if(nPushed==arraysize) { 
			DBG(printf("Ma-ru!\n")); 
			state = POBLOCK_READY_TO_CORRECT;
			pC=0; nPushed=0; pB=0;
			DBG(printf("[PoBlock::PushData] State changed to \"READY_TO_CORRECT\".\n"));
		}
	}
	void SetDataArrayPtrForCorrection(DATA* ptr) {
		if(state==POBLOCK_READY_TO_CORRECT) {
			sourceptr = ptr;
			state = POBLOCK_CORRECTING; 
			DBG(printf("[PoBlock::SetDataArrayPtrForCorrection] Starting to correct\n"));
			pC=0; nPushed=0; pB=0;
		} else {
			DBG(printf("[PoBlock::SetDataArrayPtrForCorrection] State not correct!\n"));
			DBG(printf(" must be POBLOCK_READY_TO_CORRECT to proceed\n"));
		}
	}
	void PushDataForChecking(DATA d) {
		if(state!=POBLOCK_CORRECTING) { 
			DBG(printf("[PoBlock::PushDataForChecking] State not correct!\n"));
			DBG(printf("        Must be POBLOCK_CORRECTING to proceed.\n"));
			return;
		}
		cache[pC]=d;
		pC++; nPushed++;
		std::vector<int> pos;
		std::vector<DATA> naosu;
		if((pC == blocklen*blocklen) || (nPushed==arraysize)) {
			if(pC < blocklen*blocklen) {
				unsigned int i; for(i=pC; i<blocklen*blocklen; i++) cache[i]=0;
			}
			PoBlocksum<DATA, blocklen>* blk = blocksums[pB];
			pos.clear();
			naosu.clear();
			unsigned int n = blk->CheckAndCorrect(cache, &pos, &naosu);
			if(n>0) {
				int i; for(i=0; i<n; i++) {
					unsigned int ptr = pB * blocklen * blocklen + pos.at(i);
					DATA orig = sourceptr[ptr];
					DATA elem = naosu.at(i);
					sourceptr[ptr] = elem;
					DBG(printf("[PoBlock::PushDataForChecking] Corrected array[%d] from %g to %g\n", ptr, orig, elem));
				}
			}
			pB++;
			pC=0;
		}
		if(nPushed==arraysize) {
			DBG(printf("[PoBlock::PushDataForChecking] Tried my best to correct the data.\n"));
			state = POBLOCK_READY_TO_CORRECT;
		}
	}
	~PoBlock() {
		free(cache);
	}
private:
	unsigned int nPushed;
	unsigned int pC;      // Pointer to cache
	unsigned int pB;      // Pointer to PoBlockSum
	unsigned int arraysize;
	DATA* cache;
	DATA* sourceptr; // The array to be corrected
	std::vector<PoBlocksum<DATA, blocklen>* > blocksums;
	State state;
};

int main123(int argc, char** argv) {
	printf("Very simple error correcting code\n");
	printf("Usage: ./a.out [TestArraySize] [NumOfDataToDestroy]\n");
	int sz, blksize, pct;
	if(argc >= 2) sz=atoi(argv[1]);
	else sz = 1032;
	if(argc >= 3) pct=atoi(argv[2]);
	else pct = 10;

	srand(time(0));
	double* dummyarray = (double*)malloc(sizeof(double) * sz);
	double* dummyarray_gold = (double*)malloc(sizeof(double)*sz);
	int i; for(i=0; i<sz; i++) dummyarray[i] = dummyarray_gold[i] = ((int)rand() % 99999)/100.0f;
	PoBlock<double, BLKSIZE> *block0 = new PoBlock<double, BLKSIZE>(sz);
	for(i=0; i<sz; i++) block0->PushDataForBackup(dummyarray[i]);

	/* Tamper with input data */
	int nDataToDestroy = pct;
	for(i=0; i<nDataToDestroy; i++) {
		int p = (int)(rand() % sz);
		double orig = dummyarray[p];
		double ouch = ((int)rand() % 99999) / 100.0f;
		dummyarray[p] = ouch;
		printf("array[%d] has been changed from %g to %g\n", p, orig, ouch);

	}
	block0->SetDataArrayPtrForCorrection(dummyarray);
	for(i=0; i<sz; i++) block0->PushDataForChecking(dummyarray[i]);

	int same = 0;
	for(i=0; i<sz; i++) {
		if(dummyarray[i] == dummyarray_gold[i]) same++;
	}
	printf("After correction, %d/%d are now identical.\n", same, sz);

	delete(block0);
	return 0;
};

#ifdef GSL
void* make_PoEC_data_matrix(const gsl_matrix* mat) {
	my_stopwatch_checkpoint(7);
	DBG(printf("[make_PoEC_data_matrix] Generating PoBlock ECC. Block size: %dx%d\n", BLKSIZE, BLKSIZE));
	int sz = mat->size1 * mat->size2;
	PoBlock<double, BLKSIZE>* ret = new PoBlock<double, BLKSIZE>(sz);
	int i; for(i=0; i<sz; i++) { ret->PushDataForBackup(mat->data[i]); }
	my_stopwatch_stop(7);
	return ret;
}

void* make_PoEC_data_vector(const gsl_vector* vec) {
	my_stopwatch_checkpoint(7);
	DBG(printf("[make_PoEC_data_vector] Generating PoBlock ECC. Block size: %dx%d\n", BLKSIZE, BLKSIZE));
	int sz = vec->size;
	PoBlock<double, BLKSIZE>* ret = new PoBlock<double, BLKSIZE>(sz);
	int i; for(i=0; i<sz; i++) { ret->PushDataForBackup(vec->data[i]); }
	my_stopwatch_stop(7);
	return ret;
}

void PoEC_restore_matrix(void* ecdata, gsl_matrix* mat) {
	/* Checksum is defined in the macro */
	int sz = mat->size1 * mat->size2;
	PoBlock<double, BLKSIZE>* blk = (PoBlock<double, BLKSIZE>*)ecdata;
	blk->SetDataArrayPtrForCorrection(mat->data);
	int i; for(i=0; i<sz; i++) { blk->PushDataForChecking(mat->data[i]); }
	return;
}

void PoEC_restore_vector(void* ecdata, gsl_vector* vec) {
	/* Checksum is defined in the macro */
	int sz = vec->size;
	PoBlock<double, BLKSIZE>* blk = (PoBlock<double, BLKSIZE>*)ecdata;
	blk->SetDataArrayPtrForCorrection(vec->data);
	int i; for(i=0; i<sz; i++) { blk->PushDataForChecking(vec->data[i]); }
	return;
}

#endif
