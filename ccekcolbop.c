#include "ccekcolbop.h"
#include "real.h"
#include "tommy.h"
#include "triplicate.h"
#include <sys/signal.h>
#include <signal.h>

#define FT_ROUTINES_VULN
#ifdef FT_ROUTINES_VULN
	#ifndef FT3TEST
		#define FTV_REAL_TRY(label) REAL_TRY(label)
		#define FTV_REAL_CATCH(label) REAL_CATCH(label)
		#define FTV_REAL_END(label) REAL_END(label)
	#else
		#define FTV_REAL_TRY(label) ;
		#define FTV_REAL_CATCH(label) ;
		#define FTV_REAL_END(label) ;
	#endif
#else
	#define FTV_REAL_TRY(label) ;
	#define FTV_REAL_CATCH(label) ;
	#define FTV_REAL_END(label) ;
#endif

// When fault rate is above 1e-05, we must enforce fault tolerance mechanisms on the encoding library
#define FT_ENCODE
#ifdef FT_ENCODE
	#define PROTECT_IDX_I  \
	((i!=i1-1) && (i1-1==i2-2) && (i=i1-1)), \
	((i1-1!=i2-2) && (i2-2==i) && (i1=i+1)), \
	((i2-2!=i) && (i==i1-1) && (i2=i+2))
	#define PROTECT_IDX_COLID_COLSTART \
	((colId!=colId1-1) && (colId1-1==colId2-2) && (colId=colId1-1)), \
	((colId1-1!=colId2-2) && (colId2-2==colId) && (colId1=colId+1)), \
	((colId2-2!=colId) && (colId==colId1-1) && (colId2=colId+2))
	#define PROTECT_IDX_P_COL_START \
	((pColStart!=pColStart1-1) && (pColStart1-1==pColStart2-2) && (pColStart=pColStart1-1)), \
	((pColStart1-1!=pColStart2-2) && (pColStart2-2==pColStart) && (pColStart1=pColStart+1)), \
	((pColStart2-2!=pColStart) && (pColStart==pColStart1-1) && (pColStart2=pColStart+2))
#else
#endif

// Debug outputs
#ifndef DBG
#ifdef DEBUG
#define DBG(call) {call;}
#else
#define DBG(call) {}
#endif
#endif

long poecc_num_encoded = 0;
long poecc_num_corrected = 0;

/* The doctor is examining the patient. If patient's diseases get discovered, cure them. */
/* It's possible to speculate the range of the size of the patient from the size of the 
 * doctor and vice versa, the exact sizes still need be fed into the call. */
noinline
unsigned int do_decode(double* patient, const int offsetPatient, const int lenPatient, 
	const double* doc, const int offsetDoc, const int lenDoc) /* lenDoc in num of elements */
{
	
unsigned int failed_corr = 0;
//FTV_REAL_TRY(0) {
	double colSums[BLK_LEN], rowSums[BLK_LEN], tileSum, grandTotal, rowSum=0, colSum=0;
	unsigned int isColDiff[BLK_LEN], isRowDiff[BLK_LEN]; /* If isColDiff[2] is 1, then column sum 2 is different */
	int nTiles = (lenDoc - 1) / (2*BLK_LEN + 1), pRSum, pCSum;
	int i, i1, i2, j, k, p, rowId, colId, pColStart, pRowStart;


	/* Iterating through the patient */
	#ifndef FT_ENCODE
		for(i=0; i<nTiles; i++) {
	#else
		unsigned long patient0, patient1, patient2;
		TRIPLICATE(patient, patient0, patient1, patient2);
		for(i=0, i1=1, i2=2; i<nTiles; i++, i1++, i2++, PROTECT_IDX_I) {
		TRI_RECOVER(patient0, patient1, patient2);
		if(patient != (double*)patient0) patient=(double*)patient0;
	#endif
//		printf("i=%d ", i);
		for(j=0; j<BLK_LEN; j++) { colSums[j] = rowSums[j] = isRowDiff[j] = isColDiff[j] = 0; }
		tileSum = 0;
		int pPTStart = offsetPatient + i*BLKSIZE; // Patient Tile Start
		int pPTEnd   = pPTStart + BLKSIZE; if(pPTEnd > offsetPatient + lenPatient) pPTEnd = offsetPatient + lenPatient;
		for(j=0, p=pPTStart; (j<BLKSIZE && p<pPTEnd); j++, p++) {
			double elem = patient[p];
			rowSum = rowSum + elem;
			colSums[j%BLK_LEN] += elem;
			tileSum += elem;
			grandTotal += elem;
			if( (((j+1)%BLK_LEN)==0) || p==pPTEnd-1 ) {
				rowId = j/BLK_LEN;
				rowSums[rowId] = rowSum;
				rowSum = 0;
			}
		}
		/* Calculate the colsums and rowsums of the patient */
		int lastDiffCol = 0, lastDiffRow = 0;
		int nDiffCols = 0, nDiffRows = 0;
		for(j=0; j<BLK_LEN; j++) {
			int pCSum = j + BLK_LEN + offsetDoc + i*(BLK_LEN*2+1);
			int pRSum = j + offsetDoc + i*(BLK_LEN*2+1);
			if(doc[pRSum] != rowSums[j]) {
#ifdef DEBUG_ECC
				printf("[decode] tile %d, sum of row %d different: should be %f but is %f\n",
					i, j, doc[pRSum], rowSums[j]);
#endif
				isRowDiff[j]=1;
				lastDiffRow = j; nDiffRows++;
			}
			if(doc[pCSum] != colSums[j]) {
#ifdef DEBUG_ECC
				printf("[decode] tile %d, column %d difference: should be %f but is %f\n",
					i, j, doc[pCSum], colSums[j]);
#endif
				isColDiff[j]=1;
				lastDiffCol = j; nDiffCols++;
			}
		}
		if((nDiffRows > 1) && (nDiffCols > 1)) {
#ifdef DEBUG_ECC
			printf("[decode] There are more than 1 row and more than 1 column differences, can't recover..\n");
#endif
			failed_corr += 1;
		} else {
			if(nDiffRows==1) {
				for(colId=0, pColStart=pPTStart; colId<BLK_LEN; colId++, pColStart++) {
					if(isColDiff[colId]==0) continue;

					pCSum = colId + offsetDoc + BLK_LEN + i*(2*BLK_LEN + 1);
					double sum = doc[pCSum];
					for(rowId=0, p=pColStart; (rowId<BLK_LEN && p<pPTEnd); rowId++, p+=BLK_LEN)
						{ if(rowId != lastDiffRow) {sum = sum - patient[p]; }
					p = pPTStart + lastDiffRow * BLK_LEN + colId;
#ifdef DEBUG_ECC
					printf("[debug] Correcting patient[%d] from %f to %f\n",
					p, patient[p], sum);
#endif
					patient[p] = sum;
				}
			}
			if(nDiffCols==1) {
				for(rowId=0, pRowStart=pPTStart; (rowId<BLK_LEN && pRowStart<pPTEnd); rowId++, pRowStart+=BLK_LEN) {
					if(isRowDiff[rowId]==0) continue;
					pRSum = rowId + offsetDoc + i*(2*BLK_LEN + 1);
					double sum = doc[pRSum];
					for(colId=0, p=pRowStart; (colId<BLK_LEN && p<pPTEnd); colId++, p+=1)
						{ if(colId != lastDiffCol) { sum = sum - patient[p]; } }
					p = pPTStart + lastDiffCol + rowId * BLK_LEN;
#ifdef DEBUG_ECC
					printf("[debug] Correcting patient[%d] from %f to %f\n",
					p, patient[p], sum);
#endif
					patient[p] = sum;
				}
			}
			}
		}
	}
//	}FTV_REAL_CATCH(0) {} FTV_REAL_END(0);
	return failed_corr;
}

noinline
unsigned int decode(double* patient, const int lenPatient, const double* doc) {
unsigned int failed_corr = 0;
	DBG(printf("[decode] patient=%lx lenPatient=%d doc=%lx\n", (unsigned long)patient,
		lenPatient, (unsigned long)doc));
	int jmpret = 0;
	SUPERSETJMP("Block ECC decode()");
	
	#ifdef TOMMY_H
	my_stopwatch_checkpoint(8);
	#endif

FTV_REAL_TRY(0) {
	/* Perhaps protecting patient and doctor also? */
	unsigned long p0, p1, p2;
	unsigned long d0, d1, d2;
	TRIPLICATE(patient, p0, p1, p2);
	TRIPLICATE(doc, d0, d1, d2);

	const int nTiles = ((lenPatient-1) / BLKSIZE) + 1;
	const int eccSize = ((nTiles*BLK_LEN*2) + nTiles + 1);
	const int nTilesEccEcc = ((eccSize-1) / BLKSIZE) + 1;
	const int eccEccSize = ((nTilesEccEcc*BLK_LEN*2) + nTilesEccEcc) + 1;

	/* If we need to correct the ECC itself first.... */
	int i; for(i=0; i<ECCECC; i++) {
		TRI_RECOVER(d0, d1, d2);
		do_decode((double*)d0, 0, eccSize, (const double*)d0, eccSize+eccEccSize*i, eccEccSize);
	}
	TRI_RECOVER(p0, p1, p2);
	TRI_RECOVER(d0, d1, d2);
	failed_corr = do_decode((double*)p0, 0, lenPatient, (const double*)d0, 0, eccSize); /* lenDoc in num of elements */
} FTV_REAL_CATCH(0) {} FTV_REAL_END(0);
#ifdef DEBUG
	if(failed_corr > 0) {
		print_array(patient, "Input Data", lenPatient);
		print_array(doc, "Doctor", lenPatient);
	}
#endif
	#ifdef TOMMY_H
	my_stopwatch_stop(8);
	#endif

	return failed_corr;
}

noinline
unsigned int encode(const double* array, const int len, void** backup) {
	poecc_num_encoded += len;
	int jmpret = 0;
	SUPERSETJMP("Block ECC encode()");
	printf("[encode] array=%lx len=%d backup=%lx\n", (unsigned long)array, len, (unsigned long)backup);

	/* Profiling */
	#ifdef TOMMY_H
	my_stopwatch_checkpoint(8);
	#endif
	
	/* 1. Some preparations. */
	const int nTiles = ((len-1) / BLKSIZE) + 1;
	const int eccSize = ((nTiles*BLK_LEN*2) + nTiles + 1);
	const int nTilesEccEcc = ((eccSize-1) / BLKSIZE) + 1;
	const int eccEccSize = ((nTilesEccEcc*BLK_LEN*2) + nTilesEccEcc) + 1;
	const int retSize = eccSize + eccEccSize*ECCECC;

	/* 2. Allocating space */
	double* ret = (double*)malloc(sizeof(double) * retSize);
	#ifdef FT_ENCODE
		unsigned long ret0, ret1, ret2;
		TRIPLICATE(ret, ret0, ret1, ret2);
		printf("[encode] Triplicated ret\n");
		if(jmpret != 0) {
			TRI_RECOVER(ret0, ret1, ret2);
			if(ret != (double*)ret0) ret = (double*)ret0;
		}
	#endif
	
FTV_REAL_TRY(0) {
	/* 3. Calculate the ECC of the input array and the ECC of the ECC. */
	do_encode(array, 0, len, ret, 0);
#ifdef FT_ENCODE
		TRI_RECOVER(ret0, ret1, ret2);
		if(ret != (double*)ret0) ret = (double*)ret0;
#endif
	if(ECCECC>0) do_encode(ret, 0, eccSize, ret, eccSize); // Hierarchical ECC
#ifdef DEBUG_REF
	double* ret2 = (double*)malloc(sizeof(double) * retSize);
	do_encode_2(array,0,len,ret2, 0);
	if(ECCECC>0) do_encode_2(ret2,0,eccSize, ret2, eccSize);
#endif

	if(ECCECC > 1) {
		int i, j; for(i=1; i<ECCECC; i++) {
			for(j=0; j<eccEccSize; j++) {
			ret[eccSize+eccEccSize*i+j] = ret[eccSize+j];
#ifdef DEBUG_REF
			ret2[eccSize+eccEccSize*i+j]= ret2[eccSize+j];
#endif
			}
		}
	}
#ifdef DEBUG_REF
	{ int i; int diff=0; for(i=0; i<retSize; i++) {
		double diff = ret[i] - ret2[i]; if(diff < 0) diff = -diff;
		if(diff > 1e-5) {DBG(printf(" >> [encode] ret and ret2 differ at position %d (%g vs %g.\n", i, 
			ret[i], ret2[i])); diff++;}
		}
		if(diff==0) DBG(printf(" >> [encode] ret and ret2 agree with each other!!!!\n"));
	}
#endif

	/* 4. Return */
	DBG(printf("[encode] input data size=%d, encoded data size=%d \
	(%d tiles * %d + %d + 1 + %d x %d ECCs of ECC)\n",
		len, retSize, nTiles, (2*BLK_LEN), nTiles, eccEccSize, ECCECC));
	DBG(print_array(array, "Input Data", len));
	DBG(print_array(ret, "Encoded Data", retSize));
#ifdef DEBUG_REF
	DBG(print_array(ret2, "Encoded Data 2", retSize));
#endif

	*backup = ret;
	#ifdef TOMMY_H
	my_stopwatch_stop(8);
	#endif
	return retSize;

} FTV_REAL_CATCH(0) {} FTV_REAL_END(0);
}

noinline
void do_encode_2(const double* in, const int offsetIn, const int lenIn, double* out, const int offsetOut) {
	/* 1. Some preparations. */
	int i, rowId, colId, j, k, pRowStart, pColStart, p; /* Pointer to array */
	int pCSum, pRSum, pTileSum, pGrandTotal; /* Pointer to ECC code */
	double colSum, rowSum, tileSum, grandTotal=0; /* I think we should avoid segfaults by using as few mallocs as possible. */
	const int nTiles = ((lenIn-1) / BLKSIZE) + 1;

	for(i=0; i<nTiles; i++) {
		int pStart = offsetIn + i*BLKSIZE, pEnd = pStart + BLKSIZE;
		if(pEnd > offsetIn + lenIn) pEnd=offsetIn + lenIn;

		double colSums[BLK_LEN]; for(j=0; j<BLK_LEN; j++) colSums[j]=0;
		rowId = 0;
		for(j=0, p=pStart; (j<BLKSIZE && p<pEnd); j++, p++) {
			double elem = in[p];
			rowSum += elem; colSums[j%BLK_LEN] += elem; tileSum += elem;
				if( (((j+1)%BLK_LEN)==0) || (p==pEnd-1)) {
				rowId = j/BLK_LEN;
				pRSum = offsetOut + i*(BLK_LEN*2 + 1) + rowId;
				out[pRSum] = rowSum;
				rowSum=0;
				rowId = rowId+1; 
			}
		}

		/* 2. Fill column sums */
		for(colId=0; colId<BLK_LEN; colId++) {
			pCSum = offsetOut + i*(BLK_LEN*2 + 1) + colId + BLK_LEN;
			out[pCSum] = colSums[colId];
		}
		/* 3. Fill remaining row sums */
		for(rowId; rowId < BLK_LEN; rowId++) {
			pRSum = offsetOut + i*(BLK_LEN*2 + 1) + rowId;
			out[pRSum] = 0;
		}

		pTileSum = offsetOut + i*(BLK_LEN*2 + 1) + 2*BLK_LEN;
		out[pTileSum] = tileSum;
		grandTotal += tileSum;
		tileSum = 0;
	}
	/* 3. The grand total */
	pGrandTotal = offsetOut + nTiles*(BLK_LEN*2+1);
	out[pGrandTotal] = grandTotal;
	DBG(printf(" >> [do_encode_2] Done.\n"));
}

/* Do encoding of array in[offset:len], store it into out */
noinline
void do_encode(const double* in, const int offsetIn, const int lenIn, double* out, const int offsetOut) {

FTV_REAL_TRY(0) {
	/* 1. Some preparations. */
	int i, rowId, colId, j, k, pRowStart, pColStart, p; /* Pointer to array */
	#ifdef FT_ENCODE
		int i1, i2, j1, j2, k1, k2;
		int colId1, colId2, pColStart1, pColStart2;
		unsigned long out0, out1, out2;
		TRIPLICATE(out, out0, out1, out2);
	#endif
	int pCSum, pRSum, pTileSum, pGrandTotal; /* Pointer to ECC code */
	double colSum, rowSum, tileSum, grandTotal=0; /* I think we should avoid segfaults by using as few mallocs as possible. */
	const int nTiles = ((lenIn-1) / BLKSIZE) + 1;

	#ifndef FT_ENCODE
		for(i=0; i<nTiles; i++) {
	#else
		for(i=0, i1=1, i2=2; i<nTiles; i++, i1++, i2++, PROTECT_IDX_I) {
			TRI_RECOVER(out0, out1, out2);
			if(out!=(double*)out0) out=(double*)out0;
	#endif
		int pStart = offsetIn + i*BLKSIZE, pEnd = pStart + BLKSIZE;
		if(pEnd > offsetIn + lenIn) pEnd=offsetIn + lenIn;
		/* 2.1 Column sums 
		 * We use row major here, so a[1,2,3,4,5] is a row,
		 * a[1,6,11,16,21] is a column. */
		#ifndef FT_ENCODExx
		// Replicating colId is no good
		for(colId=0, pColStart = pStart+colId; (colId<BLK_LEN && pColStart<pEnd); colId++, pColStart++) {
		#else
		for(colId=0, colId1=1, colId2=2,
			pColStart=pStart+colId, pColStart=pColStart+1,	pColStart2=pColStart+2;
			(colId<BLK_LEN && pColStart<pEnd);
			colId++, colId1++, colId2++, pColStart++, pColStart1++, pColStart2++,
			PROTECT_IDX_COLID_COLSTART) {
		#endif
			colSum = 0;
			pCSum = offsetOut + i*(BLK_LEN*2 + 1) + colId + BLK_LEN;
			/* Each element in the column. If there is no element due to being out of bounds then
			 * return zero */
			for(rowId=0, p=pColStart; (rowId<BLK_LEN && p<pEnd); rowId++, p+=BLK_LEN) { colSum += in[p]; }
			out[pCSum] = colSum;
			tileSum += colSum;
		}

		/* 2.2 Row sums */
		for(rowId=0, pRowStart = pStart+rowId*BLK_LEN; (rowId<BLK_LEN && pRowStart<pEnd); rowId++, pRowStart+=BLK_LEN) {
			rowSum = 0;
			pRSum = offsetOut + i*(BLK_LEN*2 + 1) + rowId;
			/* Each element in the row. */
			for(colId=0, p=pRowStart; (colId<BLK_LEN && p<pEnd && p<pRowStart+BLK_LEN); colId++, p+=1) { rowSum += in[p]; }
			out[pRSum] = rowSum;
		}

		/* 2.3 Tile sum */
		pTileSum = offsetOut + i*(BLK_LEN*2 + 1) + 2*BLK_LEN;
		out[pTileSum] = tileSum;
		grandTotal += tileSum;
		tileSum = 0;
	}
	
	/* 3. The grand total */
	pGrandTotal = offsetOut + nTiles*(BLK_LEN*2+1);
	out[pGrandTotal] = grandTotal;
	DBG(printf(" >> [do_encode] Done.\n"));
} FTV_REAL_CATCH(0) {} FTV_REAL_END(0);
}

/* main() is for testing purposes. */
#ifdef DEBUG2
int main(int argc, char** argv) {
	printf("Very simple error correcting code recoded in C!\n");
	printf("Usage: ./a.out [TestArraySize] [NumOfDataToDestroy]\n");
	int sz, pct;
	if(argc >= 2) sz = atoi(argv[1]); else sz=71;
	if(argc >= 3) pct = atoi(argv[2]); else pct=sz/200;

	srand(time(0));
	/* 1. Allocation and initing */
	double* array = (double*)malloc(sizeof(double) * sz);
	double* arrayGold = (double*)malloc(sizeof(double) * sz);
	int i; 
	for(i=0; i<sz; i++) array[i] = arrayGold[i] = ((int)rand() % 99999)/100.0f;
//	for(i=0; i<sz; i++) array[i] = arrayGold[i] = (i%4+2);
	void* arrayPoEcc;
	encode(array, sz, &arrayPoEcc);

	/* 2. Tamper with the input */
	for(i=0; i<pct; i++) {
		int p = (int)(rand() % sz);
		double ouch = ((int)rand() % 99999) / 100.0f;
		array[p] = ouch;
	}
	printf(">> Tempered with %d/%d inputs.\n", pct, sz);
#ifdef TAMPER_WITH_ECC_ITSELF
	int nTiles = ((sz-1)/BLKSIZE) + 1;
	int eccSize = (nTiles*2*BLK_LEN + nTiles + 1)*sizeof(double);
	for(i=0; i<pct; i++) {
		*((char*)arrayPoEcc + rand()%eccSize)=0x00;
	}
#endif


//	print_array(array, "Tampered Input", sz);

	/* 3. Attempt to recover */
	decode(array, sz, arrayPoEcc);

	/* 4. Compare recovery results */
	int same=0;
	for(i=0; i<sz; i++) if(arrayGold[i]==array[i])same++;

//	print_array(array, "Recovered Input", sz);
	DBG(printf(">> After correction, %d/%d elements were identical (recovered %d).\n", same, sz, pct-sz+same));

	return 0;
}
#endif

void POECC_SUMMARY() {
	printf("[[ PoECC Ver 2 SUMMARY]]\n");
	printf(">> Elems encoded: %d\n", poecc_num_encoded);
	printf(">> Elems corrected: %d\n", poecc_num_corrected);
}

#undef DBG
