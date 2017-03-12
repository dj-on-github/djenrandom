
/*
    genrandom - A utility to generate random numbers.
    
    Copyright (C) 2017  David Johnston

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    -----
    
    Contact. David Johnston dj@deadhat.com
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>

#include "aes128k128d.h"
#include "genrandommodel.h"
#include "rdrand.h"

#define RIGHT_VARIANCE 1
#define LEFT_VARIANCE 1

void nondeterministic_bytes(const size_t byte_len, void* byte_buf, t_rngstate *rngstate) {
    if (rngstate->rdrand_available) {
        rdrand_get_bytes_step(byte_len, byte_buf);
    }
    else if (rngstate->devurandom_available) {
        fread(rngstate->rngbits,16,1,rngstate->devrandom);
    }
    else {
        perror("Error: -s : Nondeterministic seed not available. Neither /dev/urandom nor RdRand instruction is available.");
        exit(1);
    }
}

/* AES-CTR based PRNG, reseeded with /dev/random if chosen. Outputs 16 bits at a time */ 
/* So invoked AES once every 8 outputs, since AES makes 128 bits at a time.           */

void xor16(unsigned char *a, unsigned char *b, unsigned char *c)
{
	 int i;
	for (i=0;i<16;i++)
	{
		c[i] = a[i] ^ b[i];
	}
}

/* return 16 bytes of uniform random numbers to rngstate */
int getrand16(t_rngstate* rngstate)
{
	int index;
	unsigned char realrand[16];
	unsigned char realrand2[16];
	unsigned char realrand3[16];
	unsigned char realrand4[16];
	unsigned char temprand[16];
	unsigned char temprand2[16];
	
	unsigned char out2[16];
	int i;
	int j;
	unsigned long int theint;

	/* CTR variables for random number gen */
	unsigned char out[16];
	
	/* Make a uniform Random number.              */
	/* put the random bits into a long int        */

	index = rngstate->temp;
	theint = 0;
	theint = (unsigned long int)(rngstate->rngbits[(index*2)])+ 256*((unsigned long int)(rngstate->rngbits[(index*2)+1]));

	if (rngstate->temp == 7)
	{
		rngstate->c = (rngstate->c)+1;
		if (rngstate->c == rngstate->c_max)
		{
			if (rngstate->randseed==1)
			{
                
                nondeterministic_bytes(16, realrand, rngstate);
                nondeterministic_bytes(16, realrand2, rngstate);
                nondeterministic_bytes(16, realrand3, rngstate);
                nondeterministic_bytes(16, realrand4, rngstate);
				//fread(realrand, 16, 1 , rngstate->devrandom);
				//fread(realrand2, 16, 1 , rngstate->devrandom);
				//fread(realrand3, 16, 1 , rngstate->devrandom);
				//fread(realrand4, 16, 1 , rngstate->devrandom);
                //#
				aes128k128d(rngstate->kprime, rngstate->pool0, temprand);
				xor16(realrand, temprand, temprand);
				aes128k128d(rngstate->kprime, temprand, temprand2);
				xor16(realrand2, temprand2, temprand);
				aes128k128d(rngstate->kprime, temprand, rngstate->pool0);

				aes128k128d(rngstate->kprime, rngstate->pool1, temprand);
				xor16(realrand3, temprand, temprand);
				aes128k128d(rngstate->kprime, temprand, temprand2);
				xor16(realrand4, temprand2, temprand);
				aes128k128d(rngstate->kprime, temprand, rngstate->pool1);

				for (i=0;i<16;i++) rngstate->k[i] ^= rngstate->pool0[i];
				for (i=0;i<16;i++) rngstate->v[i] ^= rngstate->pool1[i];

				rngstate->c = 0;
			}
			else
			{
				aes128k128d(rngstate->k,rngstate->v,out);
				aes128k128d(rngstate->k,out,rngstate->v);
				aes128k128d(out,rngstate->k,out2);
				aes128k128d(out2,rngstate->v,rngstate->k);
				rngstate->c = 0;
			}

		}

		aes128k128d(rngstate->k,rngstate->v,rngstate->rngbits);
		/*increment v*/
		for (j=0;j<15;j++)
		{
			if (rngstate->v[j] != 0xff)
			{
				rngstate->v[j]++;
				break;
			}
			else
			{
				rngstate->v[j]=0x00;
			}
		}
		rngstate->temp = 0;
	}
	else
	{
		rngstate->temp++;
	}

	return(theint);
}

double getNormal(t_modelstate *modelstate, t_rngstate* rngstate)
{
	unsigned long int theint;

	double u1;
	double u2;
	double y;
	double v1;
	double v2;
	double s;

	theint = getrand16(rngstate);
	u1 = ((double)theint)/65536.0;

	theint = getrand16(rngstate);
	u2 = ((double)theint)/65536.0;

	do
	{
		theint = getrand16(rngstate);
		u1 = ((double)theint)/65536.0;

		theint = getrand16(rngstate);
		u2 = ((double)theint)/65536.0;

		v1=2.0 * u1 -1.0;            /* V1=[-1,1] */
		v2=2.0 * u2 -1.0;            /* V2=[-1,1] */
		s=v1 * v1 + v2 * v2;
	} while(s >=1);
	
	/*x=sqrt(-2 * log(s) / s) * v1;*/
	y=sqrt(-2 * log(s) / s) * v2;

	/* x = modelstate->mean + (sqrt(modelstate->variance) * x);*/
	y = modelstate->mean + (sqrt(modelstate->variance) * y);
	return(y);
		
}

double getNormal_mv(double mean, double variance, t_rngstate* rngstate)
{
        unsigned long int theint;

        double u1;
        double u2;
        double y;
        double v1;
        double v2;
        double s;

        theint = getrand16(rngstate);
        u1 = ((double)theint)/65536.0;

        theint = getrand16(rngstate);
        u2 = ((double)theint)/65536.0;

        do
        {
                theint = getrand16(rngstate);
                u1 = ((double)theint)/65536.0;

                theint = getrand16(rngstate);
                u2 = ((double)theint)/65536.0;

                v1=2.0 * u1 -1.0;            /* V1=[-1,1] */
                v2=2.0 * u2 -1.0;            /* V2=[-1,1] */
                s=v1 * v1 + v2 * v2;
        } while(s >=1);

        /*x=sqrt(-2 * log(s) / s) * v1;*/
        y=sqrt(-2 * log(s) / s) * v2;

        /* x = mean + (sqrt(variance) * x); */
        y = mean + (sqrt(variance) * y);
        return(y);

}


/* Compute probablitity of a shift in state away from the centre*/
double smooth_prob_move_from_center(double t)
{
	double prob_shiftout;

	prob_shiftout = 0.5L * exp(-0.5L * t *t);
	return(prob_shiftout);
}

/* int entropysource(*j, stepsize, *k, *v);                */
/* Computes the next bit in the output bit sequence.       */
/* t is the current position                               */
/* right_stepsize is the stepsize when moving right        */
/* left_stepsize is the stepsize when moving left          */
/* k is a pointer to a 16 byte key for the simulation RNG  */
/* v is the current CTR vector for the simulation RNG      */
/* It returns the next state for j                         */
/* out is a pointer to a 16 byte array to hold rng state   */

int smoothsource(t_modelstate* modelstate, t_rngstate* rngstate)
{
	double tee;
	double pmfc;
	double randomnumber;
	int result;

	/* vars for converting from random bit to a float */
	unsigned long int theint;

	/* Get a uniform Random number.              */
	theint = getrand16(rngstate);

	tee = modelstate->t;
	
	/* Normalize it to 0-1 as a double precision floating point */
	randomnumber = theint/65536.0;

	/* evaluate the probability of a left or right */
	/* and compare with the random number to make  */
	/* a decision on a 1 or a 0.  */

	/* modelstate->sums_bias = randomnumber; */

	pmfc = smooth_prob_move_from_center(tee);

	if (pmfc > randomnumber)
	{
		/* If we're on the right, the keep going right */
		if (tee>0)
		{
			result = 1;
			tee = tee+modelstate->right_stepsize;
			modelstate->sums_bias = pmfc;
		}
		/* Else we're on the left, keep going left */
		else
		{
			result = 0;
			tee = tee-modelstate->left_stepsize;
			modelstate->sums_bias = 1.0-pmfc;
		}
	}
	/* Else we're moving towards the center */
	else
	{
		/* If we're on the right, move left */
		if (tee > 0)
		{
			result = 0;
			tee = tee-modelstate->left_stepsize;
			modelstate->sums_bias = pmfc;
		}
		/* Else we're on the left, move right */
		else
		{
			result = 1;
			tee = tee+modelstate->right_stepsize;
			modelstate->sums_bias = 1.0-pmfc;
		}
	}

	/* add noise to the step */
	if (modelstate->using_stepnoise == 1)
	{
		tee = tee + getNormal_mv(0.0,modelstate->stepnoise,rngstate);
	}

	modelstate->t = tee;

	if (modelstate->using_jfile ==1)
	{
		fprintf(modelstate->jfile,"%0.6f\n",tee);
	}

	return(result);
}

int puresource(t_modelstate* modelstate, t_rngstate* rngstate)
{
	int result;
	unsigned long int theint;

	/* get a uniform Random number.              */
	/* put the random bits into a long int        */

	theint = getrand16(rngstate);

	if (theint > 32767) result = 1;
	else result = 0;

	return(result);
}

int biasedsource(t_modelstate* modelstate, t_rngstate* rngstate)
{
	int result;
	double dthreshold;
	int threshold;
	unsigned long int theint;

	theint = getrand16(rngstate);

	dthreshold = 65536.0*(modelstate->bias);
	threshold = (int)dthreshold;
	if (theint < threshold) result = 1;
	else result = 0;

	if (modelstate->using_jfile ==1)
	{
		fprintf(modelstate->jfile,"%0.6f\n",modelstate->bias);
	}
	return(result);

}

int correlatedsource(t_modelstate *modelstate, t_rngstate* rngstate)
{
	int result;
	double dthreshold;
	int threshold;
	unsigned long int theint;
	double bias;

	/* get a uniform Random number.              */

	theint = getrand16(rngstate);

	bias = (modelstate->correlation + 1.0)/2.0;
	dthreshold = 65536.0*bias;
	threshold = (int)dthreshold;

	if (modelstate->lastbit==1)
	{
		if (theint < threshold) result = 1;
		else result = 0;
	}
	else
	{
		bias = 1.0 - bias;
		if (theint >= threshold) result = 1;
		else result = 0;
	}

	modelstate->sums_bias = bias;

	if (modelstate->using_jfile ==1)
	{
		fprintf(modelstate->jfile,"%0.6f\n",bias);
	}
	return(result);
}

int lcgsource(t_modelstate* modelstate, t_rngstate* rngstate)
{
    unsigned long long a;
    unsigned long long c;
    unsigned long long m;
    unsigned long long x;
    unsigned int truncate;
    unsigned int outbits;
    
    a = modelstate->lcg_a;
    c = modelstate->lcg_c;
    m = modelstate->lcg_m;
    x = modelstate->lcg_x;
    
    /* printf("Start: X = %llx\n",x); */
    truncate = modelstate->lcg_truncate;
    outbits = modelstate->lcg_outbits;
    
    x = (a*x + c) % m;  /* Compute next state */
    
    modelstate->lcg_x = x;
    
    /* Select the subset of bits to output */
    x = (x >> truncate);
    x = x & ((0x01ULL << outbits)-1);
    
    /* printf("End:  x = %llx, a = %llx, c = %llx, m = %llx\n",x,a,c,m); */
    return (int)x;
}

int pcg_lcg(t_modelstate* modelstate, t_rngstate* rngstate) {
    switch (modelstate->pcg_state_size) {
    case 16:
        modelstate->pcg16_state = modelstate->pcg16_state * modelstate->pcg16_multiplier;
        modelstate->pcg16_state += modelstate->pcg16_adder;
        break;
    case 32:
        modelstate->pcg32_state  = modelstate->pcg32_state * modelstate->pcg32_multiplier;
        modelstate->pcg32_state += modelstate->pcg32_adder;
        break;
    case 64:
        modelstate->pcg64_state  = modelstate->pcg64_state * modelstate->pcg64_multiplier;
        modelstate->pcg64_state += modelstate->pcg64_adder;
        break;
    case 128:
        /*modelstate->pcg16_state = modelstate->pcg16_state * modelstate->pcg16_multiplier;
        modelstate->pcg16_state += modelstate->pcg16_adder;*/
        break; 
    }
    return 0;
}

int pcg_mcg(t_modelstate* modelstate, t_rngstate* rngstate) {
    switch (modelstate->pcg_state_size) {
    case 16:
        modelstate->pcg16_state = modelstate->pcg16_state * modelstate->pcg16_multiplier;
        break;
    case 32:
        modelstate->pcg32_state  = modelstate->pcg32_state * modelstate->pcg32_multiplier;
        break;
    case 64:
        modelstate->pcg64_state  = modelstate->pcg64_state * modelstate->pcg64_multiplier;
        break;
    case 128:
        /*modelstate->pcg16_state = modelstate->pcg16_state * modelstate->pcg16_multiplier;*/
        break; 
    }
    return 0;
}

/* XSH_RS : Xor Shift, random xor shift. Stage in output function. Doesn't change state
            Outputs half the state.
*/
int pcg_xsh_rs(t_modelstate* modelstate, t_rngstate* rngstate) {
    uint64_t current;
    switch (modelstate->pcg_state_size) {
    case 16:
        current = ((modelstate->pcg16_state >> 7) ^ modelstate->pcg16_state);
        current = ((current >> ((modelstate->pcg16_state >> 14) + 3))) % 256;
        return current % 256;
    case 32:
        current = ((modelstate->pcg32_state >> 11u) ^ modelstate->pcg32_state);
        current = (current >> ((modelstate->pcg32_state >> 30u) + 11u)) % 65536;
        return current % 65536;
    case 64:
        current = ((modelstate->pcg64_state >> 22u) ^ modelstate->pcg64_state);
        current = (current >> ((modelstate->pcg64_state >> 61u) + 22u)) % 4294967296u;
        return current;
    }
    return 0;
}


/* Rotate 8 */
uint8_t rotate_uint8(uint8_t x, int n) {
    n = n % 8;
    if (n==0)
        return x;
    return ((x >> n) | (x << (8- n))); 
}

/* Rotate 16 */
uint16_t rotate_uint16(uint16_t x, int n) {
    n = n % 16;
    if (n==0)
        return x;
    return ((x >> n) | (x << (16 - n))); 
}

/* Rotate 32 */
uint32_t rotate_uint32(uint32_t x, int n) {
    n = n % 32;
    if (n==0)
        return x;
    return ((x >> n) | (x << (32 - n))); 
}

/* Rotate 64 */
uint64_t rotate_uint64(uint64_t x, int n) {
    n = n % 64;
    if (n==0)
        return x;
    return ((x >> n) | (x << (64 - n))); 
}

/* Rotate 128 bit number made of 2 64 bit uints                 */
/* The input array has two uint64s. Index 0 holds the low bits. */
void rotate_uint128(uint64_t* high_low_in, uint64_t* high_low_out, int n)
{
    uint64_t temp_high;
    uint64_t temp_low;
    
    /* limit to 128 */
    n = n % 128;

    /* If the rotate is >= 64, then it is more efficient to swap the values
       which is like shifting 64. So there n-64 shifts remaining */
    if (n  >= 64) {
        temp_low = high_low_in[1];
        temp_high = high_low_in[0];
        n = n - 64;
    }
    else {
        temp_low = high_low_in[0];
        temp_high = high_low_in[1];
    }

    if (n==0) return;

    high_low_out[0] = ((temp_low << n) | temp_high >> (64 - n));
    high_low_out[1] = ((temp_high << n) | temp_low >> (64 - n));
}

/* Shift right a 128 bit number made of 2 64 bit uints          */
/* The input array has two uint64s. Index 0 holds the low bits. */
void shift_right_uint128(uint64_t* high_low_in, uint64_t* high_low_out, int n)
{   
    /* limit to 128 */
    n = n % 128;

    /* If the rotate is > 127, returns zero */
    if (n  > 127) {
        high_low_out[0]=0;
        high_low_in[1]=0;
        return;
    }

    if (n==0) {
        high_low_out[0]=high_low_in[0];
        high_low_out[1]=high_low_in[1];
        return;
    }

    high_low_out[0] = ((high_low_in[0] >> n) | high_low_in[1] << (64 - n));
    high_low_out[1] = (high_low_in[1] >> n);
}

/* Shift left a 128 bit number made of 2 64 bit uints          */
/* The input array has two uint64s. Index 0 holds the low bits. */
void shift_left_uint128(uint64_t* high_low_in, uint64_t* high_low_out, int n)
{   
    /* limit to 128 */
    n = n % 128;

    /* If the rotate is > 127, returns zero */
    if (n  > 127) {
        high_low_out[0]=0;
        high_low_in[1]=0;
        return;
    }

    if (n==0) {
        high_low_out[0]=high_low_in[0];
        high_low_out[1]=high_low_in[1];
        return;
    }

    high_low_out[0] = (high_low_in[0] << n);
    high_low_out[1] = (high_low_in[1] << n) | (high_low_in[0] >> (64 - n));
}

/* XSH_RR : Xor Shift, random rotation. Stage in output function. Doesn't change state
            Outputs half the state.
    pcg_rotr_8(((state >> 5u) ^ state) >> 5u, state >> 13u)
    pcg_rotr_16(((state >> 10u) ^ state) >> 12u, state >> 28u);
    pcg_rotr_32(((state >> 18u) ^ state) >> 27u, state >> 59u);
    pcg_rotr_64(((state >> 29u) ^ state) >> 58u, state >> 122u);

*/
uint64_t pcg_xsh_rr(t_modelstate* modelstate, t_rngstate* rngstate) {
    uint16_t current16;
    uint32_t current32;
    uint64_t current64;
    int rotate_amount;
    switch (modelstate->pcg_state_size) {
    case 16:
        current16 = ((modelstate->pcg16_state >> 5) ^ modelstate->pcg16_state) >> 5;
        rotate_amount =  modelstate->pcg16_state >> 13u;
        /*printf("  16: in-%04x",(unsigned int)(current16)); */
        current16 = rotate_uint8(current16,rotate_amount); 
        /*printf("  %04x",(unsigned int)(current16)); */
        return current16 % 256;
    case 32:
        current32 = ((modelstate->pcg32_state >> 10) ^ modelstate->pcg32_state) >> 12;
        rotate_amount =  modelstate->pcg32_state >> 28u;
        /*printf("  32: in-%08x",(unsigned int)(current32));*/
        current32 = rotate_uint16(current32,rotate_amount); 
        /*printf("  %08x",(unsigned int)(current32));*/
        return current32 & 0xffffffff;
    case 64:
        current64 = ((modelstate->pcg64_state >> 18) ^ modelstate->pcg64_state) >> 27;
        rotate_amount =  modelstate->pcg64_state >> 59u;
        /*printf("  64: in-%016llx",(current64));*/
        current64 = current64 & 0xffffffff;
        current64 = rotate_uint32(current64,rotate_amount); 
        /*printf(" rot(%d) ",rotate_amount);
        printf("  %016llx",(current64));*/
        return current64;
    }
    return 0;
}

/* The PCG main routine which calls the internal state update function
   followed by the output function. The internal state update can be
   MCG or LCG. The output function can be XSH_RS or XSH_RS. In this
   implementation, the output size is always half the state size.*/
int pcgsource(t_modelstate* modelstate, t_rngstate* rngstate)
{   
    /* If we have no more bits to return, get a new value
     * from the PCG algorithm. Else Shift bits out.
     */
    if (modelstate->pcg_index == 0) {

        /* Update the internal state */
        if (modelstate->pcg_alg == PCG_MCG) {
            pcg_mcg(modelstate,rngstate);
            /*printf("New MCG 0x%016llx",modelstate->pcg64_state);*/
        }
        else {
            pcg_lcg(modelstate,rngstate);  
            /*printf("New LCG 0x%016llx",modelstate->pcg64_state); */  
        } 

        /* Call the output function */
        if (modelstate->pcg_of == XSH_RS) {
            modelstate->pcg_output = pcg_xsh_rs(modelstate,rngstate);
            /*printf("  RS 0x%08x\n",(unsigned int)(modelstate->pcg_output));*/
        }
        else if (modelstate->pcg_of == XSH_RR) {
            modelstate->pcg_output = pcg_xsh_rr(modelstate,rngstate);
            /*printf("  RR 0x%08x\n",(unsigned int)(modelstate->pcg_output));*/
        }
        
           
        /* Return one bit from the result and reset the index */
        modelstate->pcg_index = (modelstate->pcg_state_size >> 1) -1;
        return (int)(modelstate->pcg_output & 0x1);
    }
    else {
        modelstate->pcg_index -= 1;
        modelstate->pcg_output = modelstate->pcg_output >> 1;
        return (int)(modelstate->pcg_output & 0x1);
    }
}

int xorshiftsource(t_modelstate* modelstate, t_rngstate* rngstate)
{
    unsigned int x;
    
    if (modelstate->xorshift_size == 32)
    {
	    x = modelstate->xorshift_state_a;
	    x = x ^ (x << 13);
	    x = x ^ (x >> 17);
	    x = x ^ (x << 5);
	    modelstate->xorshift_state_a = x;
	    return (int)x;
    }
    else /* xorshift128 */
    {
    	x = modelstate->xorshift_state_d;
    	
	    x = x ^ (x << 11);
	    x = x ^ (x >> 8);
	    
	    modelstate->xorshift_state_d = modelstate->xorshift_state_c;
	    modelstate->xorshift_state_c = modelstate->xorshift_state_b;
	    modelstate->xorshift_state_b = modelstate->xorshift_state_a;

	    x = x ^ modelstate->xorshift_state_a;
	    x = x ^ (modelstate->xorshift_state_a >> 19);	
	    modelstate->xorshift_state_a = x;
	    return (int)x;
    }

}


int filesource(t_modelstate* modelstate, t_rngstate* rngstate)
{
	int result;
	int int_c;
	int doneit;

	doneit = 0;

	/* Fetch characters. '0' gives a 0, '1' gives a 1,
 	 * EOF returns a 0, others are skipped
 	 */

	do {
		int_c = fgetc(modelstate->infile);
		if (int_c == EOF) {
			result = 0;
			doneit = 1;
			rngstate->reached_eof = 1;
		}
		else if (int_c == '0') {
			result=0;
			doneit = 1;	
		}
		else if (int_c == '1') {
			result=1;
			doneit = 1;
		}
	} while (doneit==0);

	return(result);
}

int filesourcehex(t_modelstate* modelstate, t_rngstate* rngstate)
{
	int result;
	int int_c;
	int doneit;
	unsigned char myfilechar;
	doneit = 0;

	/* Fetch characters until we get a hex one.
 	 * Others are skipped. If we already have
 	 * a character and the bits are being shifted out
 	 * then the else clause is executed and another
 	 * bit is shifted out.
 	 */
	if (rngstate -> fileindex == 0)
	{
		do
		{
			int_c = fgetc(modelstate->infile);
			
			if (int_c == EOF)
			{
				rngstate->filechar = 0;
				doneit = 1;
				rngstate->reached_eof = 1;
			}
			else if (int_c == '0') {
				rngstate->filechar=0;
				doneit = 1;	
			}
			else if (int_c == '1') {
				rngstate->filechar=1;
				doneit = 1;
			}
			else if (int_c == '2') {
				rngstate->filechar=2;
				doneit = 1;
			}
			else if (int_c == '3') {
				rngstate->filechar=3;
				doneit = 1;
			}
			else if (int_c == '4') {
				rngstate->filechar=4;
				doneit = 1;
			}
			else if (int_c == '5') {
				rngstate->filechar=5;
				doneit = 1;
			}
			else if (int_c == '6') {
				rngstate->filechar=6;
				doneit = 1;
			}
			else if (int_c == '7') {
				rngstate->filechar=7;
				doneit = 1;
			}
			else if (int_c == '8') {
				rngstate->filechar=8;
				doneit = 1;
			}
			else if (int_c == '9') {
				rngstate->filechar=9;
				doneit = 1;
			}
			else if ((int_c == 'a')||(int_c == 'A')) {
				rngstate->filechar=10;
				doneit = 1;
			}
			else if ((int_c == 'b')||(int_c == 'B')) {
				rngstate->filechar=11;
				doneit = 1;
			}
			else if ((int_c == 'c')||(int_c == 'C')) {
				rngstate->filechar=12;
				doneit = 1;
			}
			else if ((int_c == 'd')||(int_c == 'D')) {
				rngstate->filechar=13;
				doneit = 1;
			}
			else if ((int_c == 'e')||(int_c == 'E')) {
				rngstate->filechar=14;
				doneit = 1;
			}
			else if ((int_c == 'f')||(int_c == 'F')) {
				rngstate->filechar=15;
				doneit = 1;
			}
		} while (doneit==0);
		result = ((rngstate->filechar & 0x08)>>3) & 0x01;
		rngstate->fileindex++;
	}
	else
	{
		myfilechar = rngstate->filechar;
		myfilechar = myfilechar << 1;
		rngstate->filechar = myfilechar;
		result = ((myfilechar & 0x08)>>3) & 0x01;
		rngstate->fileindex++;
		if ((rngstate->fileindex)==4)
		{
			rngstate->fileindex = 0;
		}
	}
	return(result);
}

double normalsource(t_modelstate *modelstate, t_rngstate* rngstate)
{
        double result;
        
	result = getNormal(modelstate, rngstate);

        return(result);
}

/*****************************************/
/* Init routines. Initializee the models */
/*****************************************/

void smoothinit(t_modelstate* modelstate, t_rngstate* rngstate)
{
	int i;
	unsigned char realrand[16];
	unsigned char out[16];

	/* Set k and v to some arbitary values */
	for (i=0;i<16;i++)
	{
		rngstate->k[i] = (unsigned char)55;
		rngstate->v[i] = (unsigned char)33;
	}

	if (rngstate->randseed==1)
	{
        nondeterministic_bytes(16, realrand, rngstate);
        for (i=0;i<16;i++) rngstate->k[i] = rngstate->k[i] ^ realrand[i];
        
        nondeterministic_bytes(16, realrand, rngstate);
        for (i=0;i<16;i++) rngstate->v[i] = rngstate->v[i] ^ realrand[i];
        
        //else /* Use /dev/random to read random data on Linux */
        //{
        //    fread(realrand, 16, 1 , rngstate->devrandom);
        //    for (i=0;i<16;i++) rngstate->k[i] = rngstate->k[i] ^ realrand[i];
        //
        //    fread(realrand, 16, 1 , rngstate->devrandom);
        //    for (i=0;i<16;i++) rngstate->v[i] = rngstate->v[i] ^ realrand[i];
        //}
	}

	/* Make k and v head off into the weeds */
	aes128k128d(rngstate->k,rngstate->v,out);
	aes128k128d(rngstate->k,out,rngstate->v);
	aes128k128d(out,rngstate->v,rngstate->k);
	aes128k128d(rngstate->k,rngstate->v,rngstate->kprime);
	
	rngstate->temp = 0;

	modelstate->t = 0;
	rngstate->reached_eof = 0;

}

void pureinit(t_modelstate* modelstate, t_rngstate* rngstate)
{
	unsigned char out[16];

	smoothinit(modelstate, rngstate);
	if (rngstate->randseed==1)
	{
	    nondeterministic_bytes(16, rngstate->rngbits, rngstate);
		/*fread(rngstate->rngbits,16,1,rngstate->devrandom);*/
	}
	else
	{
		aes128k128d(rngstate->k,rngstate->v,out);
		aes128k128d(out,rngstate->k,rngstate->rngbits);
	}

}
void biasedinit(t_modelstate *modelstate, t_rngstate *rngstate)
{
	unsigned char out[16];

	smoothinit(modelstate, rngstate);
	if (rngstate->randseed==1)
	{
	    nondeterministic_bytes(16, rngstate->rngbits, rngstate);
		/*fread(rngstate->rngbits,16,1,rngstate->devrandom);*/
	}
	else
	{
		aes128k128d(rngstate->k,rngstate->v,out);
		aes128k128d(out,rngstate->k,rngstate->rngbits);
	}

}
void correlatedinit(t_modelstate *modelstate, t_rngstate *rngstate)
{
	unsigned char out[16];

	smoothinit(modelstate, rngstate);
	if (rngstate->randseed==1)
	{
	    nondeterministic_bytes(16, rngstate->rngbits, rngstate);
		/*fread(rngstate->rngbits,16,1,rngstate->devrandom);*/
	}
	else
	{
		aes128k128d(rngstate->k,rngstate->v,out);
		aes128k128d(out,rngstate->k,rngstate->rngbits);
	}

}

void lcginit(t_modelstate *modelstate, t_rngstate *rngstate)
{
    unsigned long long state;
    
    modelstate->lcg_mask = (1 << (modelstate->lcg_m)) -1;
    
	if (rngstate->randseed==1)
	{
	    nondeterministic_bytes(sizeof(unsigned long long), &state, rngstate);
		modelstate->lcg_x = (modelstate->lcg_x ^ state) & (modelstate->lcg_mask);
	}
}

void pcginit(t_modelstate *modelstate, t_rngstate *rngstate)
{

    /* Multiplier constants for the LCG part for the different state sizes */
    /*modelstate->pcg8_multiplier    = 0x8D;*/
    modelstate->pcg16_multiplier   = 0x321D;
    modelstate->pcg32_multiplier   = 0x2C9277B5;
    modelstate->pcg64_multiplier   = 0x5851F42D4C957F2D;
    modelstate->pcg128_multiplier[0]  = 0x4385DF649FCCF645;
    modelstate->pcg128_multiplier[1]  = 0x2360ED051FC65DA4;
    
    /* Adder constants for the LCG part for the different state sizes */
    /*modelstate->pcg8_adder    = 0x4D;*/
    modelstate->pcg16_adder   = 0xBB75;
    modelstate->pcg32_adder   = 0xAC564B05;
    modelstate->pcg64_adder   = 0x14057B7EF767814F;
    modelstate->pcg128_adder[0]  = 0x14057B7EF767814F;
    modelstate->pcg128_adder[1]  = 0x5851F42D4C957F2D;
        
    /* Arbitary random start state gathered by running genrandom -s */
    /*modelstate->pcg8_state    = 0x6058;*/
    modelstate->pcg16_state   = 0x7109;
    modelstate->pcg32_state   = 0x4B55CD7E;
    modelstate->pcg64_state   = 0x71A56ADF832343F1;
    modelstate->pcg128_state[0] = 0x3C14A2F859655FEB;
    modelstate->pcg128_state[1] = 0x49BC03F5388BAD9D;
    
	if (rngstate->randseed==1)
	{
	    /*nondeterministic_bytes(sizeof(uint8_t), &(modelstate->pcg8_state), rngstate);*/
	    nondeterministic_bytes(sizeof(uint16_t), &(modelstate->pcg16_state), rngstate);
	    nondeterministic_bytes(sizeof(uint32_t), &(modelstate->pcg32_state), rngstate);
	    nondeterministic_bytes(sizeof(uint64_t), &(modelstate->pcg64_state), rngstate);
	    nondeterministic_bytes(sizeof(uint64_t), &(modelstate->pcg128_state[0]), rngstate);
	    nondeterministic_bytes(sizeof(uint64_t), &(modelstate->pcg128_state[1]), rngstate);
	}
}


void xorshiftinit(t_modelstate *modelstate, t_rngstate *rngstate)
{   
    modelstate->xorshift_state_a = 0xA634716A; 
    modelstate->xorshift_state_b = 0x998FCD1F;
    modelstate->xorshift_state_c = 0x6A9B90FE;
    modelstate->xorshift_state_d = 0x7344E998;
    
	if (rngstate->randseed==1)
	{
	    nondeterministic_bytes(sizeof(unsigned int), &(modelstate->xorshift_state_a), rngstate);
	    nondeterministic_bytes(sizeof(unsigned int), &(modelstate->xorshift_state_b), rngstate);
	    nondeterministic_bytes(sizeof(unsigned int), &(modelstate->xorshift_state_c), rngstate);
	    nondeterministic_bytes(sizeof(unsigned int), &(modelstate->xorshift_state_d), rngstate);
	}
}

void fileinit(t_modelstate* modelstate, t_rngstate* rngstate)
{
	rngstate->reached_eof = 0;
	rngstate->filechar = 0x00;
	rngstate->fileindex = 0;
	smoothinit(modelstate, rngstate);
}

void normalinit(t_modelstate *modelstate, t_rngstate *rngstate)
{
	unsigned char out[16];

	smoothinit(modelstate, rngstate);
	if (rngstate->randseed==1)
	{
	    nondeterministic_bytes(16, rngstate->rngbits, rngstate);
		/*fread(rngstate->rngbits,16,1,rngstate->devrandom);*/
	}
	else
	{
		aes128k128d(rngstate->k,rngstate->v,out);
		aes128k128d(out,rngstate->k,rngstate->rngbits);
	}
}

