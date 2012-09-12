#include "ccekcolbop.h"
#define FT_ROUTINES_VULN
#ifdef FT_ROUTINES_VULN
#include "real.h"
#else
#include "spoof_real.h"
#endif

long poecc_num_encoded_f = 0;
long poecc_num_corrected_f = 0;

/* The doctor is examining the patient. If patient's diseases get discovered, cure them. */
/* It's possible to speculate the range of the size of the patient from the size of the 
 * doctor and vice versa, the exact sizes still need be fed into the call. */
noinline
void do_decode_float(float* patient, const int offsetPatient, const int lenPatient, 
	const float* doc, const int offsetDoc, const int lenDoc) /* lenDoc in num of elements */
{
REAL_TRY(0) {
	float colSums[BLK_LEN], rowSums[BLK_LEN], tileSum, grandTotal, rowSum=0, colSum=0;
	unsigned int isColDiff[BLK_LEN], isRowDiff[BLK_LEN]; /* If isColDiff[2] is 1, then column sum 2 is different */
	int nTiles = (lenDoc - 1) / (2*BLK_LEN + 1), pRSum, pCSum;
	int i, j, k, p, rowId, colId, pColStart, pRowStart;


	/* Iterating through the patient */
	for(i=0; i<nTiles; i++) {
		for(j=0; j<BLK_LEN; j++) { colSums[j] = rowSums[j] = isRowDiff[j] = isColDiff[j] = 0; }
		tileSum = 0;
		int pPTStart = offsetPatient + i*BLKSIZE; // Patient Tile Start
		int pPTEnd   = pPTStart + BLKSIZE; if(pPTEnd > offsetPatient + lenPatient) pPTEnd = offsetPatient + lenPatient;
		for(j=0, p=pPTStart; (j<BLKSIZE && p<pPTEnd); j++, p++) {
			float elem = patient[p];
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
#ifdef DEBUG
				printf("[decode_float] tile %d, sum of row %d different: should be %f but is %f\n",
					i, j, doc[pRSum], rowSums[j]);
#endif
				isRowDiff[j]=1;
				lastDiffRow = j; nDiffRows++;
			}
			if(doc[pCSum] != colSums[j]) {
#ifdef DEBUG
				printf("[decode_float] tile %d, column %d difference: should be %f but is %f\n",
					i, j, doc[pCSum], colSums[j]);
#endif
				isColDiff[j]=1;
				lastDiffCol = j; nDiffCols++;
			}
		}
		if((nDiffRows > 1) && (nDiffCols > 1)) {
#ifdef DEBUG
			printf("[decode_float] There are more than 1 row and more than 1 column differences, can't recover..\n");
#endif
		} else {
			if(nDiffRows==1) {
				for(colId=0, pColStart=pPTStart; colId<BLK_LEN; colId++, pColStart++) {
					if(isColDiff[colId]==0) continue;

					pCSum = colId + offsetDoc + BLK_LEN + i*(2*BLK_LEN + 1);
					float sum = doc[pCSum];
					for(rowId=0, p=pColStart; (rowId<BLK_LEN && p<pPTEnd); rowId++, p+=BLK_LEN)
						{ if(rowId != lastDiffRow) {sum = sum - patient[p]; }
					p = pPTStart + lastDiffRow * BLK_LEN + colId;
#ifdef DEBUG
					printf("[debug] Correcting patient[%d] from %f to %f\n",
					p, patient[p], sum);
#endif
					patient[p] = sum; poecc_num_corrected_f++;
				}
			}
			if(nDiffCols==1) {
				for(rowId=0, pRowStart=pPTStart; (rowId<BLK_LEN && pRowStart<pPTEnd); rowId++, pRowStart+=BLK_LEN) {
					if(isRowDiff[rowId]==0) continue;
					pRSum = rowId + offsetDoc + i*(2*BLK_LEN + 1);
					float sum = doc[pRSum];
					for(colId=0, p=pRowStart; (colId<BLK_LEN && p<pPTEnd); colId++, p+=1)
						{ if(colId != lastDiffCol) { sum = sum - patient[p]; } }
					p = pPTStart + lastDiffCol + rowId * BLK_LEN;
#ifdef DEBUG
					printf("[debug] Correcting patient[%d] from %f to %f\n",
					p, patient[p], sum);
#endif
					patient[p] = sum; poecc_num_corrected_f++;
				}
			}
			}
		}
	}
	}REAL_CATCH(0) {} REAL_END(0);
}

noinline
void decode_float(float* patient, const int lenPatient, const float* doc) {
REAL_TRY(0) {
	const int nTiles = ((lenPatient-1) / BLKSIZE) + 1;
	const int eccSize = ((nTiles*BLK_LEN*2) + nTiles + 1);
	const int nTilesEccEcc = ((eccSize-1) / BLKSIZE) + 1;
	const int eccEccSize = ((nTilesEccEcc*BLK_LEN*2) + nTilesEccEcc) + 1;

	/* If we need to correct the ECC itself first.... */
	int i; for(i=0; i<ECCECC; i++) {
		do_decode_float((void*)doc, 0, eccSize, doc, eccSize+eccEccSize*i, eccEccSize);
	}

	do_decode_float(patient, 0, lenPatient, doc, 0, eccSize); /* lenDoc in num of elements */
} REAL_CATCH(0) {} REAL_END(0);
}

noinline
void encode_float(const float* array, const int len, void** backup) {
	poecc_num_encoded_f += len;
REAL_TRY(0) {
	/* 1. Some preparations. */
	const int nTiles = ((len-1) / BLKSIZE) + 1;
	const int eccSize = ((nTiles*BLK_LEN*2) + nTiles + 1);
	const int nTilesEccEcc = ((eccSize-1) / BLKSIZE) + 1;
	const int eccEccSize = ((nTilesEccEcc*BLK_LEN*2) + nTilesEccEcc) + 1;
	const int retSize = eccSize + eccEccSize*ECCECC;

	/* 2. Allocating space */
	float* ret = (float*)malloc(sizeof(float) * retSize);
	
	/* 3. Calculate the ECC of the input array and the ECC of the ECC. */
	do_encode_float(array, 0, len, ret, 0);
	if(ECCECC>0) do_encode_float(ret, 0, eccSize, ret, eccSize); // Hierarchical ECC
#ifdef DEBUG_REFERENCE_ENCODER
	float* ret2 = (float*)malloc(sizeof(float) * retSize);
	do_encode_2_float(array,0,len,ret2, 0);
	if(ECCECC>0) do_encode_2_float(ret2,0,eccSize, ret2, eccSize);
#endif

	if(ECCECC > 1) {
		int i, j; for(i=1; i<ECCECC; i++) {
			for(j=0; j<eccEccSize; j++) {
			ret[eccSize+eccEccSize*i+j] = ret[eccSize+j];
#ifdef DEBUG_REFERENCE_ENCODER
			ret2[eccSize+eccEccSize*i+j]= ret2[eccSize+j];
#endif
			}
		}
	}
#ifdef DEBUG_REFERENCE_ENCODER
	{ int i; int diff=0; for(i=0; i<retSize; i++) {
		if(ret[i]!=ret2[i]) {
			printf(" >> [encode_float] ret and ret2 differ at position %d (%g vs %g).\n", i,
				ret[i], ret2[i]);
			diff++;
			
			}
		}
		if(diff==0) printf(" >> [encode_float] ret and ret2 agree with each other!!!!\n");
	}
#endif

	/* 4. Return */
#ifdef DEBUG_REFERENCE_ENCODER
	printf("[encode_float] input data size=%d, encoded data size=%d \
	(%d tiles * %d + %d + 1 + %d x %d ECCs of ECC)\n",
		len, retSize, nTiles, (2*BLK_LEN), nTiles, eccEccSize, ECCECC);
	print_array(array, "Input Data", len);
	print_array(ret, "Encoded Data", retSize);
	print_array(ret2, "Encoded Data 2", retSize);
#endif

	*backup = ret;
	#ifdef TOMMY_H
	my_stopwatch_stop(8);
	#endif

} REAL_CATCH(0) {} REAL_END(0);

}

noinline
void do_encode_2_float(const float* in, const int offsetIn, const int lenIn, float* out, const int offsetOut) {
	/* 1. Some preparations. */
	int i, rowId, colId, j, k, pRowStart, pColStart, p; /* Pointer to array */
	int pCSum, pRSum, pTileSum, pGrandTotal; /* Pointer to ECC code */
	float colSum, rowSum, tileSum, grandTotal=0; /* I think we should avoid segfaults by using as few mallocs as possible. */
	const int nTiles = ((lenIn-1) / BLKSIZE) + 1;

	for(i=0; i<nTiles; i++) {
		int pStart = offsetIn + i*BLKSIZE, pEnd = pStart + BLKSIZE;
		if(pEnd > offsetIn + lenIn) pEnd=offsetIn + lenIn;

		float colSums[BLK_LEN]; for(j=0; j<BLK_LEN; j++) colSums[j]=0;
		rowId = 0;
		for(j=0, p=pStart; (j<BLKSIZE && p<pEnd); j++, p++) {
			float elem = in[p];
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
void do_encode_float(const float* in, const int offsetIn, const int lenIn, float* out, const int offsetOut) {
REAL_TRY(0) {
	/* 1. Some preparations. */
	int i, rowId, colId, j, k, pRowStart, pColStart, p; /* Pointer to array */
	int pCSum, pRSum, pTileSum, pGrandTotal; /* Pointer to ECC code */
	float colSum, rowSum, tileSum, grandTotal=0; /* I think we should avoid segfaults by using as few mallocs as possible. */
	const int nTiles = ((lenIn-1) / BLKSIZE) + 1;

	for(i=0; i<nTiles; i++) {
		int pStart = offsetIn + i*BLKSIZE, pEnd = pStart + BLKSIZE;
		if(pEnd > offsetIn + lenIn) pEnd=offsetIn + lenIn;
		/* 2.1 Column sums 
		 * We use row major here, so a[1,2,3,4,5] is a row,
		 * a[1,6,11,16,21] is a column. */
		for(colId=0, pColStart = pStart+colId; (colId<BLK_LEN && pColStart<pEnd); colId++, pColStart++) {
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
	printf(" >> [do_encode_float] Done. Grand Total=%g\n", grandTotal);
} REAL_CATCH(0) {} REAL_END(0);
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
	float* array = (float*)malloc(sizeof(float) * sz);
	float* arrayGold = (float*)malloc(sizeof(float) * sz);
	int i; 
	for(i=0; i<sz; i++) array[i] = arrayGold[i] = ((int)rand() % 99999)/100.0f;
//	for(i=0; i<sz; i++) array[i] = arrayGold[i] = (i%4+2);
	void* arrayPoEcc;
	encode_float(array, sz, &arrayPoEcc);

	/* 2. Tamper with the input */
	for(i=0; i<pct; i++) {
		int p = (int)(rand() % sz);
		float ouch = ((int)rand() % 99999) / 100.0f;
		array[p] = ouch;
	}
	printf(">> Tempered with %d/%d inputs.\n", pct, sz);
#ifdef TAMPER_WITH_ECC_ITSELF
	int nTiles = ((sz-1)/BLKSIZE) + 1;
	int eccSize = (nTiles*2*BLK_LEN + nTiles + 1)*sizeof(float);
	for(i=0; i<pct; i++) {
		*((char*)arrayPoEcc + rand()%eccSize)=0x00;
	}
#endif


//	print_array(array, "Tampered Input", sz);

	/* 3. Attempt to recover */
	decode_float(array, sz, arrayPoEcc);

	/* 4. Compare recovery results */
	int same=0;
	for(i=0; i<sz; i++) if(arrayGold[i]==array[i])same++;

//	print_array(array, "Recovered Input", sz);
	printf(">> After correction, %d/%d elements were identical (recovered %d).\n", same, sz, pct-sz+same);

	return 0;
}
#endif

__attribute__((noinline))
void POECC_FLOAT_SUMMARY() {
	printf("[[ PoECC Ver 2 SUMMARY (single precision ver) ]]\n");
	printf(">> Elems encoded: %d\n", poecc_num_encoded_f);
	printf(">> Elems corrected: %d\n", poecc_num_corrected_f);
}
