/*  CRC-32b version 1.03 by Craig Bruce, 27-Jan-94
**
**  Based on "File Verification Using CRC" by Mark R. Nelson in Dr. Dobb's
**  Journal, May 1992, pp. 64-67.  This program DOES generate the same CRC
**  values as ZMODEM and PKZIP
**
**  v1.00: original release.
**  v1.01: fixed printf formats.
**  v1.02: fixed something else.
**  v1.03: replaced CRC constant table by generator function.
*/
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
unsigned long getcrc(FILE*);
unsigned long mygetcrc(char*, unsigned long);
void crcgen();

unsigned long mycrc_get();
void mycrc_reset();
void mycrc_pushd(double);
void mycrc_pushc(char);
