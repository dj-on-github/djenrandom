
#include <stdlib.h>
#include <stdio.h>
#include "aes128k128d.h"
#include "cmac.h"

/****************************************/
/* leftshifts by 1 bit a 128 bit vector */
/****************************************/

/* Take note that in the NIST CMAC specification
   the MSB is byte index 0, so the carries go the
   other way to the way you might expect */

void leftshift128(unsigned char *a, unsigned char *out) {
	unsigned char x;
	int carry;
	int i;
	carry=0;

	for (i=15; i >= 0 ; i--)
	{	
		x = a[i] << 1;
		if (carry == 1)
		{
			x = x | 0x01;
		}
		else
		{
			x = x & 0xfe;
		}

		out[i] = x;

		carry = 0;
		if ((a[i] & 0x80)==0x80) carry = 1;
	}

}

/********************************************************/
/* int cmac()                                           */
/* Computes the 128 bit CMAC of the plaintext           */
/********************************************************/

int cmac(   unsigned char *key,
		unsigned char *k1,
		unsigned char *k2,
		unsigned char *plaintext,
                unsigned int length,
                unsigned char *t)
{

	unsigned char l[16];
	unsigned char lshl[16];

	unsigned char r[16] = {0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x87 };
	unsigned char zeroes[16] = {0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00 };

	unsigned char c[16] = {0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00 };
	unsigned char temp[16] = {0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00 };
	unsigned char mn[16] = {0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00 };
	unsigned char *m;

	unsigned int n;
	unsigned int i;
	unsigned int remainder;

	/* Compute subkeys K1 and K2 */

	aes128k128d(key, zeroes, l);
	if ((l[0] & 0x80) == 0x80)  /* if MSB of L is set */
	{
		leftshift128(l,lshl);
		xor_128(lshl,r,k1);
	}
	else
	{
		leftshift128(l,k1);
	};

	
	if ((k1[0] & 0x80) == 0x80)  /* if MSB of K1 is set */
	{
		leftshift128(k1,lshl);
		xor_128(lshl,r,k2);
	}
	else
	{
		leftshift128(k1,k2);
	};

	/* Compute the MAC */

	n = length/16;
	remainder = length % 16;

	if (remainder != 0) { n++; };

	for (i=0; i<16; i++) {c[i] = 0x00;};	/* Set c0 to all zeroes */

	for (i = 0; i< n-1; i++)
	{
		m = (unsigned char *)(plaintext+(16*i));
		xor_128(c,m,temp);
		aes128k128d(key, temp, c);
	};

					/* Final block */
	if (remainder == 0)
	{
		m = plaintext+(16*(n-1));
		xor_128(m, k1, mn);
		xor_128(c, mn, temp);
		aes128k128d(key,temp, t);
	}
	else
	{
		for (i=0; i<16; i++) temp[i]=0x00; /* clear temp */

						/* copy last fragment into temp */
		m = plaintext+(16*(n-1));
		for (i=0; i<remainder; i++)
		{
			temp[i] = m[i];	
		}

						/* set the leftmost bit of the padding */
		temp[remainder] = 0x80;

		xor_128(temp, k2, mn);
		xor_128(mn, c, temp);

		aes128k128d(key, temp, t);

	}

	

};

