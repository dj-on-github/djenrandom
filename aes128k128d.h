/****************************************************************/
/* 802.11 AES Block Cipher                                      */
/* Implements AES128. I.E. 128 bit key with 128 bit data.       */
/*                                                              */
/* Originally (c) 2002, David Johnston                          */
/* This code has been placed in the public domain by the author */
/* Author: David Johnston                                       */
/* Email (home) : dj@deadhat.com                                */
/* Email (general) : david.johnston@ieee.org                    */
/* Email (work) : dj.johnston@intel.com                         */
/****************************************************************/
 
#include <stdlib.h>
void xor_128(unsigned char *a, unsigned char *b, unsigned char *out);
void aes128k128d(unsigned char *key, unsigned char *data, unsigned char *ciphertext);


