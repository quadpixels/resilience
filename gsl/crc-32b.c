#include "crc-32b.h"

unsigned long crcTable[256];

/* global calls */
static unsigned long g_crc;

void mycrc_reset() { crcgen(); g_crc=0xFFFFFFFF; }

void mycrc_pushd(double d) {
	char* k = (char*)&d;
	int i; for(i=0; i<sizeof(double)/sizeof(char); i++) {
		char c=((*(k+i))&0xFF);
		g_crc=((g_crc >> 8) & 0x00FFFFFF) ^ crcTable[(g_crc^c) & 0xFF];
	}
}

unsigned long mycrc_get() { return g_crc ^ 0xFFFFFFFF; }

/**
 * Use this code to get a CRC code.
 */
unsigned long mygetcrc(char* data, unsigned long len) {
	register unsigned long crc;
	printf("[mygetcrc] File size: %ld\n", len);
	crc=0xFFFFFFFF;
	unsigned long i; for(i=0; i<len; i++) {
		crc = ((crc>>8) & 0x00FFFFFF) ^ crcTable[ (crc^data[i]) & 0xFF ];
//		printf("[mycrc] CRC: %08lx\n", crc);
	}
	return( crc ^ 0xFFFFFFFF);
}

/****************************************************************************/
unsigned long getcrc( FILE* fp )
{
	register unsigned long crc;
	int c;

	crc = 0xFFFFFFFF;
	while( (c=getc(fp)) != EOF ) {
		crc = ((crc>>8) & 0x00FFFFFF) ^ crcTable[ (crc^c) & 0xFF ];
//		printf("[getcrc] CRC: %08lx\n", crc);
	}
	return( crc^0xFFFFFFFF );
}

/****************************************************************************/
void crcgen( )
{
//	if(g_crcGened==1) return;
	if(crcTable[0]!=0) return;
	printf("[crcgen]\n");
	unsigned long	crc, poly;
	int	i, j;

	poly = 0xEDB88320L;
	for (i=0; i<256; i++) {
		crc = i;
		for (j=8; j>0; j--) {
			if (crc&1) {
				crc = (crc >> 1) ^ poly;
			} else {
				crc >>= 1;
			}
		}
		crcTable[i] = crc;
	}
}

