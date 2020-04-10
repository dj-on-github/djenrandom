
/*
    djrandom - A utility to generate random numbers.
    
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
#include "djenrandommodel.h"
#include "rdrand.h"
#include "markov2p.h"

#define RIGHT_VARIANCE 1
#define LEFT_VARIANCE 1

extern int verbose_mode;

int bit_reverse_table[256] = {
    0x00,0x80,0x40,0xc0,0x20,0xa0,0x60,0xe0,0x10,0x90,0x50,0xd0,0x30,0xb0,0x70,0xf0,
    0x08,0x88,0x48,0xc8,0x28,0xa8,0x68,0xe8,0x18,0x98,0x58,0xd8,0x38,0xb8,0x78,0xf8,
    0x04,0x84,0x44,0xc4,0x24,0xa4,0x64,0xe4,0x14,0x94,0x54,0xd4,0x34,0xb4,0x74,0xf4,
    0x0c,0x8c,0x4c,0xcc,0x2c,0xac,0x6c,0xec,0x1c,0x9c,0x5c,0xdc,0x3c,0xbc,0x7c,0xfc,
    0x02,0x82,0x42,0xc2,0x22,0xa2,0x62,0xe2,0x12,0x92,0x52,0xd2,0x32,0xb2,0x72,0xf2,
    0x0a,0x8a,0x4a,0xca,0x2a,0xaa,0x6a,0xea,0x1a,0x9a,0x5a,0xda,0x3a,0xba,0x7a,0xfa,
    0x06,0x86,0x46,0xc6,0x26,0xa6,0x66,0xe6,0x16,0x96,0x56,0xd6,0x36,0xb6,0x76,0xf6,
    0x0e,0x8e,0x4e,0xce,0x2e,0xae,0x6e,0xee,0x1e,0x9e,0x5e,0xde,0x3e,0xbe,0x7e,0xfe,
    0x01,0x81,0x41,0xc1,0x21,0xa1,0x61,0xe1,0x11,0x91,0x51,0xd1,0x31,0xb1,0x71,0xf1,
    0x09,0x89,0x49,0xc9,0x29,0xa9,0x69,0xe9,0x19,0x99,0x59,0xd9,0x39,0xb9,0x79,0xf9,
    0x05,0x85,0x45,0xc5,0x25,0xa5,0x65,0xe5,0x15,0x95,0x55,0xd5,0x35,0xb5,0x75,0xf5,
    0x0d,0x8d,0x4d,0xcd,0x2d,0xad,0x6d,0xed,0x1d,0x9d,0x5d,0xdd,0x3d,0xbd,0x7d,0xfd,
    0x03,0x83,0x43,0xc3,0x23,0xa3,0x63,0xe3,0x13,0x93,0x53,0xd3,0x33,0xb3,0x73,0xf3,
    0x0b,0x8b,0x4b,0xcb,0x2b,0xab,0x6b,0xeb,0x1b,0x9b,0x5b,0xdb,0x3b,0xbb,0x7b,0xfb,
    0x07,0x87,0x47,0xc7,0x27,0xa7,0x67,0xe7,0x17,0x97,0x57,0xd7,0x37,0xb7,0x77,0xf7,
    0x0f,0x8f,0x4f,0xcf,0x2f,0xaf,0x6f,0xef,0x1f,0x9f,0x5f,0xdf,0x3f,0xbf,0x7f,0xff};

uint64_t choose_exponent(uint64_t start, t_rngstate* rngstate) {
    uint64_t e;

    e = start;
    do {
        if ((getrand64(rngstate) & 0x01) == 1) return e;
        e = e-1;
    } while (e > 0);
    return ((uint64_t)0);
}


double get_rand_double(t_rngstate* rngstate) {
    uint64_t start;
    uint64_t mantissa;
    uint64_t exponent;
    //uint64_t sign;
    uint64_t x;
    double *f;
    double result;

    start = 1022;
    //int i;
    //for (i=0;i<1000;i++) {
        mantissa = (getrand64(rngstate) & 0x07ffffffffffff) | 0x08000000000000;
        exponent = choose_exponent(start,rngstate);
        //sign = getrand64(rngstate) & 0x01;
        //sign = 0;
        //x = (sign << 63) | ((exponent & 0x7ff) << 52) | mantissa;
        x = ((exponent & 0x7ff) << 52) | mantissa;
        f = (double *)&x;
        //fprintf(stderr,"%f  exponent=%llu\n",*f,exponent);
    //}
    result = *f;
    //fprintf(stderr,"  GET_RAND_DOUBLT = %f\n",result);
    fflush(stdout);
    return result;
}


void nondeterministic_bytes(const size_t byte_len, void* byte_buf, t_rngstate *rngstate) {
    if (rngstate->rdrand_available) {
        rdrand_get_bytes_step(byte_len, byte_buf);
        rdrand_get_bytes_step(byte_len, byte_buf);
    }
    else if (rngstate->devurandom_available) {
        fread(rngstate->rngbits,16,1,rngstate->devrandom);
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

/* return 16 bits of uniform random numbers to rngstate */
int getrand16(t_rngstate* rngstate)
{
    //int i;
    int j;
    int index;
    //unsigned char realrand[16];
    //unsigned char realrand2[16];
    //unsigned char realrand3[16];
    //unsigned char realrand4[16];
    //unsigned char temprand[16];
    //unsigned char temprand2[16];
    
    /* CTR variables for random number gen */
    unsigned char out[16];
    unsigned char out2[16];

    //int i;
    //int j;
    unsigned long int theint;

    
    /* Make a uniform Random number.              */
    /* put the random bits into a long int        */

    index = rngstate->temp;
    theint = 0;
    theint = (unsigned long int)(rngstate->rngbits[(index*2)])+ (((unsigned long int)(rngstate->rngbits[(index*2)+1]))<<8);
    //if (verbose_mode) {
    //    fprintf(stderr, "  theint      = %lu\n",theint);
    //    fprintf(stderr, "  theint_pre1 = %02lx\n", (unsigned long int)(rngstate->rngbits[(index*2)]));
    //    fprintf(stderr, "  theint_pre2 = %02lx\n", (unsigned long int)(rngstate->rngbits[(index*2)+1]));
    // 
    //}
   
    if (rngstate->temp < 7) {
        rngstate-> temp++;
        return(theint);
    } 
    else
    {
        if (rngstate->randseed==1) {
            nondeterministic_bytes(16, rngstate->rngbits, rngstate);
            //nondeterministic_bytes(16, realrand2, rngstate);
            //nondeterministic_bytes(16, realrand3, rngstate);
            //nondeterministic_bytes(16, realrand4, rngstate);
            rngstate->temp = 0;
        }
        else {
        
            rngstate->c = (rngstate->c)+1;
            if (rngstate->c == rngstate->c_max)
            {
                aes128k128d(rngstate->k,rngstate->v,out);
                aes128k128d(rngstate->k,out,rngstate->v);
                aes128k128d(out,rngstate->k,out2);
                aes128k128d(out2,rngstate->v,rngstate->k);
                rngstate->c = 0;
            }

            aes128k128d(rngstate->k,rngstate->v,rngstate->rngbits);
            /*increment v*/
            for (j=0;j<15;j++) {
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
        rngstate->temp++;
    }
    return(theint);
}

uint64_t getrand64(t_rngstate* rngstate) {
    uint64_t therand;
    int quarterrand;
    int i;
    
    for(i=0;i<4;i++) {
        quarterrand = getrand16(rngstate);
        therand = (therand << 16) | (quarterrand & 0xffff);
    }
    return therand;
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
    double maxp;
    double entropy;

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

    /* modelstate->bias = randomnumber; */

    pmfc = smooth_prob_move_from_center(tee);

    if (pmfc > randomnumber)
    {
        /* If we're on the right, the keep going right */
        if (tee>0)
        {
            result = 1;
            tee = tee+modelstate->right_stepsize;
            modelstate->bias = pmfc;
        }
        /* Else we're on the left, keep going left */
        else
        {
            result = 0;
            tee = tee-modelstate->left_stepsize;
            modelstate->bias = 1.0-pmfc;
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
            modelstate->bias = pmfc;
        }
        /* Else we're on the left, move right */
        else
        {
            result = 1;
            tee = tee+modelstate->right_stepsize;
            modelstate->bias = 1.0-pmfc;
        }
    }

    /* add noise to the step */
    if (modelstate->using_stepnoise == 1)
    {
        tee = tee + getNormal_mv(0.0,modelstate->stepnoise,rngstate);
    }

    modelstate->t = tee;
    
    if (pmfc > 0.5)
        maxp = pmfc;
    else
        maxp = 1.0-pmfc;

    entropy = -log(maxp)/log(2);
    modelstate->averageentropy = entropy;
    modelstate->n = (modelstate->n)+1;
 
    if (modelstate->using_jfile ==1) {
        fprintf(modelstate->jfile,"%0.6f\n",entropy);
    }
    return(result);
}

int puresource(t_modelstate* modelstate, t_rngstate* rngstate)
{
    int result;
    unsigned long int theint;
    double entropy;

    /* get a uniform Random number.              */
    /* put the random bits into a long int        */

    theint = getrand16(rngstate);

    if (theint > 32767) result = 1;
    else result = 0;

    modelstate->averageentropy = 1.0;
    modelstate->n = (modelstate->n)+1;

    if (modelstate->using_jfile ==1)
    {
        entropy = 1.0;
        fprintf(modelstate->jfile,"%0.6f\n",entropy);
    }

    return(result);
}

int biasedsource(t_modelstate* modelstate, t_rngstate* rngstate)
{
    int result;
    double dthreshold;
    int threshold;
    unsigned long int theint;
    double maxp;
    double entropy;

    theint = getrand16(rngstate);

    dthreshold = 65536.0*(modelstate->bias);
    threshold = (int)dthreshold;
    if (theint < threshold) result = 1;
    else result = 0;

    if (modelstate->bias > 0.5)
        maxp = modelstate->bias;
    else
        maxp = 1.0 - modelstate->bias;
    
    entropy = -log(maxp)/log(2);
    modelstate->averageentropy = entropy;
    modelstate->n = (modelstate->n)+1;
 
    if (modelstate->using_jfile ==1) {
        fprintf(modelstate->jfile,"%0.6f\n",entropy);
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
    double maxp;
    double entropy;
    int delta;

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

    modelstate->bias = bias;

    if (bias > 0.5)
        maxp = bias;
    else
        maxp = 1.0 - bias;

    entropy = -log(maxp)/log(2);

    modelstate->n = (modelstate->n)+1;
    delta = entropy - (modelstate->averageentropy);
    modelstate->averageentropy = (modelstate->averageentropy) + (delta/(modelstate->n));

    if (modelstate->using_jfile == 1) {
        fprintf(modelstate->jfile,"%0.6f\n",entropy);
    }
    
    return(result);
}

int markov2psource(t_modelstate *modelstate, t_rngstate* rngstate)
{
    int result;
    double p01_dthreshold;
    double p10_dthreshold;
    int p10_threshold;
    int p01_threshold;
    unsigned long int theint;
    double p01;
    double p10;
    double maxp;
    double entropy;
    double bias;
    double delta;

    /* get a uniform Random number.              */

    theint = getrand16(rngstate);

    p01 = modelstate->p01;
    p10 = modelstate->p10;
    
    p01_dthreshold = 65536.0*p01;
    p10_dthreshold = 65536.0*p10;
    p01_threshold = (int)p01_dthreshold;
    p10_threshold = (int)p10_dthreshold;
    
            //fprintf(stderr,"bias = %f\n",modelstate.bias);
            //fprintf(stderr,"correlation = %f\n",modelstate.correlation);
            //fprintf(stderr,"p01 = %f\n",p01);
            //fprintf(stderr,"p10 = %f\n",p10);
            
    if (modelstate->lastbit==1)
    {
        if (theint < p10_threshold) result = 0;
        else result = 1;
        bias = p10;
    }
    else
    {
        if (theint < p01_threshold) result = 1;
        else result = 0;
        bias = p01;
    }

    modelstate->bias = bias;

    if (bias > 0.5)
        maxp = bias;
    else
        maxp = 1.0 - bias;

    entropy = -log(maxp)/log(2);

    modelstate->n = (modelstate->n)+1;
    delta = entropy - (modelstate->averageentropy);
    modelstate->averageentropy = (modelstate->averageentropy) + (delta/(modelstate->n));

    if (modelstate->using_jfile == 1) {
        fprintf(modelstate->jfile,"%0.6f\n",entropy);
    }

    return(result);
}

int markov2pfastsource(t_modelstate *modelstate, t_rngstate* rngstate)
{
    int symbol;
    int rsymbol;
    //int phase;
    unsigned long int theint;
    //double maxp;
    //double entropy;
    //double bias;
    //double delta;

    /* get a uniform Random number.              */

    theint = ((getrand16(rngstate) << 16) | (getrand16(rngstate))) & 0x000fffff;


            //fprintf(stderr,"bias = %f\n",modelstate.bias);
            //fprintf(stderr,"correlation = %f\n",modelstate.correlation);
            //fprintf(stderr,"p01 = %f\n",p01);
            //fprintf(stderr,"p10 = %f\n",p10);

    //phase = modelstate->markov2p_phase;
    //symbol = modelstate->markov2p_symbol;

    if (modelstate->lastbit==1)
    {
        symbol = modelstate->sampletable1[theint];
    } else {
        symbol = modelstate->sampletable0[theint];
    }
   
    //modelstate->bias = bias;
    modelstate->lastbit=(symbol >> 7) & 1;
    rsymbol=bit_reverse_table[symbol];
    return(rsymbol);
}

int markovsigmoidsource(t_modelstate *modelstate, t_rngstate* rngstate)
{
    int state;
    int chain_len;
    double p_left;
    double therand;
    int     result;
    double  maxp;
    double  entropy;
    double delta;
    
    /* get a uniform Random floating point number.              */
    therand = get_rand_double(rngstate);
    
    state = modelstate->sigmoid_state;
    p_left = modelstate->chain[state];
    chain_len = modelstate->states;
    
    if ((state > 0) && (p_left > therand)) { // move left
        state = state-1;
        result = 0;
    } else if ((state == 0) && (p_left > therand)) { // Stay at left
        state = 0;
        result = 0;
    } else if ((state < (chain_len-1)) && (p_left < therand)) { // move right
        state = state+1;
        result = 1;
    } else if ((state == (chain_len-1)) && (p_left < therand)) { // stay at right
        state = chain_len-1;
        result = 1;
    }
    else {     // do not be here
        if (verbose_mode==4) {
            fprintf(stderr,"Impossibru State state=%d  pleft=%0.3f therand=%0.4f\n",state,p_left,therand);
        }
        state = 0;
        result = 1;
    }

    modelstate->sigmoid_state = state;
    modelstate->sigmoid_bias = p_left;
    
    if (p_left > 0.5)
        maxp = p_left;
    else
        maxp = 1.0 - p_left;

    entropy = -log(maxp)/log(2);

    modelstate->n = (modelstate->n)+1;
    delta = entropy - (modelstate->averageentropy);
    modelstate->averageentropy = (modelstate->averageentropy) + (delta/(modelstate->n));

    if (modelstate->using_jfile == 1) {
        fprintf(modelstate->jfile,"%0.6f\n",entropy);
    }
    return(result);
}


int sinbiassource(t_modelstate *modelstate, t_rngstate* rngstate)
{
    int result;
    double dthreshold;
    int threshold;
    unsigned long int theint;
    double bias;
    double period;
    double amplitude;
    double offset;
    int t;
    double maxp;
    double entropy;
    double delta;

    /* get a uniform Random number.              */

    theint = getrand16(rngstate);

    period = modelstate->sinbias_period;
    amplitude = modelstate->sinbias_amplitude;
    offset = modelstate->sinbias_offset;
    t = modelstate->time;
    
    bias = offset + (amplitude*(sin(2.0*M_PI*t/period)));
    
    if (bias < 0.0) bias =0.0;
    if (bias > 1.0) bias = 1.0;
    
    dthreshold = 65536.0*bias;
    threshold = (int)dthreshold;


    if (theint < threshold) result = 1;
    else result = 0;

    modelstate->sinbias_bias = bias;

    if (bias > 0.5)
        maxp = bias;
    else
        maxp = 1.0 - bias;

    entropy = -log(maxp)/log(2);

    modelstate->n = (modelstate->n)+1;
    delta = entropy - (modelstate->averageentropy);
    modelstate->averageentropy = (modelstate->averageentropy) + (delta/(modelstate->n));

    if (modelstate->using_jfile == 1) {
        fprintf(modelstate->jfile,"%0.6f\n",entropy);
    }

    modelstate->time = t+1;
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
    /*printf ("  X = %llx, index=%d\n",x,modelstate->lcg_index);*/
    if ((modelstate->lcg_index)==0) {
        truncate = modelstate->lcg_truncate;
        outbits = modelstate->lcg_outbits;
    
        x = (a*x + c) % m;  /* Compute next state */
    
        modelstate->lcg_x = x;
    
        /* Select the subset of bits to output */
        x = (x >> truncate);
        modelstate->lcg_output = x & ((0x01ULL << outbits)-1);
        
        modelstate->lcg_index = outbits-1;
        return (modelstate->lcg_output & 1);
    }
    else {
        modelstate->lcg_index -= 1;
        modelstate->lcg_output = modelstate->lcg_output >> 1;
        return (modelstate->lcg_output & 1);
    }
    /* fprintf(stderr,"End:  x = %llx, a = %llx, c = %llx, m = %llx\n",x,a,c,m); */
    modelstate->averageentropy = 0.0;
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
        /*fprintf(stderr,"  16: in-%04x",(unsigned int)(current16)); */
        current16 = rotate_uint8(current16,rotate_amount); 
        /*fprintf(stderr,"  %04x",(unsigned int)(current16)); */
        return current16 % 256;
    case 32:
        current32 = ((modelstate->pcg32_state >> 10) ^ modelstate->pcg32_state) >> 12;
        rotate_amount =  modelstate->pcg32_state >> 28u;
        /*fprintf(stderr,"  32: in-%08x",(unsigned int)(current32));*/
        current32 = rotate_uint16(current32,rotate_amount); 
        /*fprintf(stderr,"  %08x",(unsigned int)(current32));*/
        return current32 & 0xffffffff;
    case 64:
        current64 = ((modelstate->pcg64_state >> 18) ^ modelstate->pcg64_state) >> 27;
        rotate_amount =  modelstate->pcg64_state >> 59u;
        /*fprintf(stderr,"  64: in-%016llx",(current64));*/
        current64 = current64 & 0xffffffff;
        current64 = rotate_uint32(current64,rotate_amount); 
        /*fprintf(stderr," rot(%d) ",rotate_amount);
        fprintf(stderr,"  %016llx",(current64));*/
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
    modelstate->averageentropy = 0.0;
    
    if (modelstate->pcg_index == 0) {

        /* Update the internal state */
        if (modelstate->pcg_alg == PCG_MCG) {
            pcg_mcg(modelstate,rngstate);
            /*fprintf(stderr,"New MCG 0x%016llx",modelstate->pcg64_state);*/
        }
        else {
            pcg_lcg(modelstate,rngstate);  
            /*fprintf(stderr,"New LCG 0x%016llx",modelstate->pcg64_state); */  
        } 

        /* Call the output function */
        if (modelstate->pcg_of == XSH_RS) {
            modelstate->pcg_output = pcg_xsh_rs(modelstate,rngstate);
            /*fprintf(stderr,"  RS 0x%08x\n",(unsigned int)(modelstate->pcg_output));*/
        }
        else if (modelstate->pcg_of == XSH_RR) {
            modelstate->pcg_output = pcg_xsh_rr(modelstate,rngstate);
            /*fprintf(stderr,"  RR 0x%08x\n",(unsigned int)(modelstate->pcg_output));*/
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
    
    modelstate->averageentropy = 0.0;
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
    modelstate->averageentropy = 0.0;

    return(result);
}

int filesourcehex(t_modelstate* modelstate, t_rngstate* rngstate)
{
    int result;
    int int_c;
    int doneit;
    unsigned char myfilechar;
    doneit = 0;
    
    modelstate->averageentropy = 0.0;

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

int filesourcebinary(t_modelstate* modelstate, t_rngstate* rngstate)
{
    int result;
    int int_c;
    /*int doneit;*/
    unsigned char myfilechar;
    /*doneit = 0;*/

    modelstate->averageentropy = 0.0;

    rngstate->reached_eof = 0;
    
    /* Fetch a byte then shift out the bits */ 
    if (rngstate -> fileindex == 0)
    {
        int_c = fgetc(modelstate->infile);
        if (int_c == EOF)
            {
                rngstate->filechar = 0;
                /*doneit = 1;*/
                rngstate->reached_eof = 1;
                return 0;
            }
            
        rngstate->filechar=int_c;
        rngstate->fileindex=7;
        result = ((rngstate->filechar & 0x80)>>7) & 0x01;
    }
    else
    {
        myfilechar = rngstate->filechar;
        myfilechar = myfilechar << 1;
        rngstate->filechar = myfilechar;
        result = ((myfilechar & 0x80)>>7) & 0x01;
        rngstate->fileindex--;
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

void init_rng(t_rngstate* rngstate) {
    int i;
    unsigned char realrand[16];
    unsigned char out[16];

    unsigned char cmac_key_a[16];
    unsigned char cmac_key_b[16];
    unsigned char cmac_k1[16];
    unsigned char cmac_k2[16];
    /* Set k and v to some arbitary values */
    for (i=0;i<16;i++)
    {
        rngstate->k[i] = (unsigned char)55;
        rngstate->v[i] = (unsigned char)33;
    }

    if (rngstate->got_detseed==1)
    {
        for(i=0;i<16;i++) {
            cmac_key_a[i]=(unsigned char)i;
            cmac_key_b[i]=(unsigned char)(i+16);
        }
        cmac(cmac_key_a, cmac_k1, cmac_k2, rngstate->detseed, strlen(rngstate->detseed), rngstate->k);
        cmac(cmac_key_b, cmac_k1, cmac_k2, rngstate->detseed, strlen(rngstate->detseed), rngstate->v);
    }

    if (rngstate->randseed==1)
    {
        nondeterministic_bytes(16, realrand, rngstate);
        for (i=0;i<16;i++) rngstate->k[i] = rngstate->k[i] ^ realrand[i];
        
        nondeterministic_bytes(16, realrand, rngstate);
        for (i=0;i<16;i++) rngstate->v[i] = rngstate->v[i] ^ realrand[i];
    }

    /* Make k and v head off into the weeds */
    aes128k128d(rngstate->k,rngstate->v,out);
    aes128k128d(rngstate->k,out,rngstate->v);
    aes128k128d(out,rngstate->v,rngstate->k);
    aes128k128d(rngstate->k,rngstate->v,rngstate->kprime);
    
    rngstate->temp = 0;
    rngstate->reached_eof = 0;

}

void smoothinit(t_modelstate* modelstate, t_rngstate* rngstate)
{
    unsigned char out[16];

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
    
    modelstate->t = 0;        
}

void pureinit(t_modelstate* modelstate, t_rngstate* rngstate)
{
    unsigned char out[16];

    //smoothinit(modelstate, rngstate);
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

    //smoothinit(modelstate, rngstate);
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

    //smoothinit(modelstate, rngstate);
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

void markov2pinit(t_modelstate *modelstate, t_rngstate *rngstate)
{
    unsigned char out[16];
    double epsilon;
    double mcv_prob;
    uint64_t mcv;

    //smoothinit(modelstate, rngstate);
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

    // Deal with the 3 parameter types
    if (modelstate->gotentropy==1) {
        epsilon = pow(2.0,-50);
        pick_point(&(modelstate->p01),&(modelstate->p10),modelstate->entropy,epsilon,modelstate->bitwidth,rngstate);
        modelstate->bias = modelstate->p01/(modelstate->p10+modelstate->p01);
        modelstate->correlation = 1.0 - modelstate->p01 - modelstate->p10;   
    }
    else if ((modelstate->gotbias == 1) || (modelstate->gotcorrelation==1)) {
        if (modelstate->gotbias==0){
             modelstate->bias = 0.5;
        }
        
        if (modelstate->gotcorrelation==0) {
            modelstate->correlation = 0.0;
        }
        modelstate->p01 = modelstate->bias * (1.0 - modelstate->correlation);
        modelstate->p10 = (1.0-modelstate->bias)*(1.0-modelstate->correlation);
        //fprintf(stderr,"bias = %f\n",modelstate.bias);
        //fprintf(stderr,"correlation = %f\n",modelstate.correlation);
        //fprintf(stderr,"p01 = %f\n",modelstate.p01);
        //fprintf(stderr,"p10 = %f\n",modelstate.p10);
        modelstate->entropy=p_to_entropy(modelstate->p01, modelstate->p10,8, &mcv_prob, &mcv);
    } else if ((modelstate->gotp01 == 1) || (modelstate->gotp10==1))  {
        if (modelstate->gotp01==0) modelstate->p01 = 0.5;
        if (modelstate->gotp10==0) modelstate->p10 = 0.5;
        modelstate->correlation = 1.0 - modelstate->p10 - modelstate->p01;
        modelstate->bias  = modelstate->p10/(modelstate->p10+modelstate->p01);
        modelstate->entropy=p_to_entropy(modelstate->p01, modelstate->p10,8, &mcv_prob, &mcv);
    } else if (modelstate->gotentropy==0) {
        modelstate->entropy = 1.0;
        modelstate->p01 = 0.5;
        modelstate->p10 = 0.5;
        modelstate->bias = 0.5;
        modelstate->correlation = 0.0;
    }

    if (verbose_mode > 0)
    {
        fprintf(stderr,"model=markov_2_param\n");
        fprintf(stderr,"  bias            = %f\n",modelstate->bias);
        fprintf(stderr,"  correlation     = %f\n",modelstate->correlation);
        fprintf(stderr,"  p01             = %f\n",modelstate->p01);
        fprintf(stderr,"  p10             = %f\n",modelstate->p10);
        fprintf(stderr,"  entropy         = %f\n",modelstate->entropy);
        //fprintf(stderr,"  MCV Prob        = %f\n",lmcv_prob);
        fprintf(stderr,"  Bits per symbol = %d\n",modelstate->bitwidth);
    }

    if (modelstate->using_yaml > 0)
    {
        fprintf(modelstate->yaml_file,"model:markov_2_param\n");
        if (modelstate->using_ofile==1)
            fprintf(modelstate->yaml_file,"filename:%s\n",modelstate->filename);
        fprintf(modelstate->yaml_file,"bias:%f\n",modelstate->bias);
        fprintf(modelstate->yaml_file,"scc:%f\n",modelstate->correlation);
        fprintf(modelstate->yaml_file,"p01:%f\n",modelstate->p01);
        fprintf(modelstate->yaml_file,"p10:%f\n",modelstate->p10);
        fprintf(modelstate->yaml_file,"minentropy:%f\n",modelstate->entropy);
        fprintf(modelstate->yaml_file,"bits_per_symbol:1\n");
    }

    if (modelstate->using_json > 0)
    {
        fprintf(modelstate->json_file,"{\n");
        fprintf(modelstate->json_file,"  \"model\"           : \"markov_2_param\"\n");
        if (modelstate->using_ofile==1)
            fprintf(modelstate->json_file,"  \"filename\"        : \"%s\"\n",modelstate->filename);
        fprintf(modelstate->json_file,"  \"bias\"            : \"%f\"\n",modelstate->bias);
        fprintf(modelstate->json_file,"  \"scc\"             : \"%f\"\n",modelstate->correlation);
        fprintf(modelstate->json_file,"  \"p01\"             : \"%f\"\n",modelstate->p01);
        fprintf(modelstate->json_file,"  \"p10\"             : \"%f\"\n",modelstate->p10);
        fprintf(modelstate->json_file,"  \"minentropy\"      : \"%f\"\n",modelstate->entropy);
        fprintf(modelstate->json_file,"  \"bits_per_symbol\" : \"1\"\n");
        fprintf(modelstate->json_file,"}\n");
    }

}

void markov2pfastinit(t_modelstate *modelstate, t_rngstate *rngstate)
{
    unsigned char out[16];
    double epsilon;
    double p01_dthreshold;
    double p10_dthreshold;
    int *sampletable0;
    int *sampletable1;
    double mcv_prob;
    uint64_t mcv;

    //smoothinit(modelstate, rngstate);
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

    // Deal with the 3 parameter types
    if (modelstate->gotentropy==1) {
        epsilon = pow(2.0,-50);
        pick_point(&(modelstate->p01),&(modelstate->p10),modelstate->entropy,epsilon,modelstate->bitwidth,rngstate);
        modelstate->bias = modelstate->p01/(modelstate->p10+modelstate->p01);
        modelstate->correlation = 1.0 - modelstate->p01 - modelstate->p10;   
    }
    else if ((modelstate->gotbias == 1) || (modelstate->gotcorrelation==1)) {
        if (modelstate->gotbias==0){
             modelstate->bias = 0.5;
        }
        
        if (modelstate->gotcorrelation==0) {
            modelstate->correlation = 0.0;
        }
        modelstate->p01 = modelstate->bias * (1.0 - modelstate->correlation);
        modelstate->p10 = (1.0-modelstate->bias)*(1.0-modelstate->correlation);
        //fprintf(stderr,"bias = %f\n",modelstate.bias);
        //fprintf(stderr,"correlation = %f\n",modelstate.correlation);
        //fprintf(stderr,"p01 = %f\n",modelstate.p01);
        //fprintf(stderr,"p10 = %f\n",modelstate.p10);
        modelstate->entropy=p_to_entropy(modelstate->p01, modelstate->p10,8, &mcv_prob, &mcv);
    } else if ((modelstate->gotp01 == 1) || (modelstate->gotp10==1))  {
        if (modelstate->gotp01==0) modelstate->p01 = 0.5;
        if (modelstate->gotp10==0) modelstate->p10 = 0.5;
        modelstate->correlation = 1.0 - modelstate->p10 - modelstate->p01;
        modelstate->bias  = modelstate->p10/(modelstate->p10+modelstate->p01);
        modelstate->entropy=p_to_entropy(modelstate->p01, modelstate->p10,8, &mcv_prob, &mcv);
    } else if (modelstate->gotentropy==0) {
        modelstate->entropy = 1.0;
        modelstate->p01 = 0.5;
        modelstate->p10 = 0.5;
        modelstate->bias = 0.5;
        modelstate->correlation = 0.0;
    }

    p01_dthreshold = 65536.0*(modelstate->p01);
    p10_dthreshold = 65536.0*(modelstate->p10);
    modelstate->p01_threshold = (int)p01_dthreshold;
    modelstate->p10_threshold = (int)p10_dthreshold;
   
    // Make lookup table

    if (verbose_mode==1) fprintf(stderr,"Calling make_sample_table()\n");

    make_sample_table(modelstate->p01, modelstate->p10, modelstate->bitwidth, &sampletable0, &sampletable1);

    modelstate->sampletable0 = sampletable0;
    modelstate->sampletable1 = sampletable1;

    if (verbose_mode > 0)
    {
        fprintf(stderr,"model=markov_2_param\n");
        fprintf(stderr,"  bias            = %f\n",modelstate->bias);
        fprintf(stderr,"  correlation     = %f\n",modelstate->correlation);
        fprintf(stderr,"  p01             = %f\n",modelstate->p01);
        fprintf(stderr,"  p10             = %f\n",modelstate->p10);
        fprintf(stderr,"  entropy         = %f\n",modelstate->entropy);
        //fprintf(stderr,"  MCV Prob        = %f\n",lmcv_prob);
        fprintf(stderr,"  Bits per symbol = %d\n",modelstate->bitwidth);
    }

    if (modelstate->using_json > 0)
    {
        fprintf(modelstate->json_file,"{\n");
        fprintf(modelstate->json_file,"  \"model\"           : \"markov_2_param\"\n");
        if (modelstate->using_ofile==1)
            fprintf(modelstate->json_file,"  \"filename\"        : \"%s\"\n",modelstate->filename);
        fprintf(modelstate->json_file,"  \"bias\"            : \"%f\"\n",modelstate->bias);
        fprintf(modelstate->json_file,"  \"scc\"             : \"%f\"\n",modelstate->correlation);
        fprintf(modelstate->json_file,"  \"p01\"             : \"%f\"\n",modelstate->p01);
        fprintf(modelstate->json_file,"  \"p10\"             : \"%f\"\n",modelstate->p10);
        fprintf(modelstate->json_file,"  \"minentropy\"      : \"%f\"\n",modelstate->entropy);
        fprintf(modelstate->json_file,"  \"bits_per_symbol\" : \"%d\"\n",modelstate->bitwidth);
        fprintf(modelstate->json_file,"}\n");
    }

    if (modelstate->using_yaml > 0)
    {
        fprintf(modelstate->yaml_file,"model:markov_2_param\n");
        if (modelstate->using_ofile==1)
            fprintf(modelstate->yaml_file,"filename:%s\n",modelstate->filename);
        fprintf(modelstate->yaml_file,"bias:%f\n",modelstate->bias);
        fprintf(modelstate->yaml_file,"scc:%f\n",modelstate->correlation);
        fprintf(modelstate->yaml_file,"p01:%f\n",modelstate->p01);
        fprintf(modelstate->yaml_file,"p10:%f\n",modelstate->p10);
        fprintf(modelstate->yaml_file,"minentropy:%f\n",modelstate->entropy);
        fprintf(modelstate->yaml_file,"bits_per_symbol:%d\n",modelstate->bitwidth);
    }

}

void markovsigmoidinit(t_modelstate *modelstate, t_rngstate *rngstate)
{
    int i;
    unsigned char out[16];
    double low;
    double high;
    double states;
    double range;
    double stepsize;
    double x;
    double themin;
    double themax;
    double h;
    
    //smoothinit(modelstate, rngstate);
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

    // Starting state
    modelstate->sigmoid_state = (modelstate->states) >> 1;
    

    //
    // Allocate space for Markov Chain
    //
    modelstate->chain = malloc((modelstate->states)*sizeof(double));
        
    //
    // Set up the chain transition probabilities
    // based on the curve chosen and the range.
    //
    low = modelstate->min_range;
    high = modelstate->max_range;
    states = modelstate->states;
    range = high-low;
    
    // Compute stepsize
    stepsize = (range/(states-1));
    
    // Compute x distance between each state
    stepsize = (range/(states-1));    
    
    if (modelstate->curve == CURVE_FLAT) {
        for (i=0;i<modelstate->states;i++) {
            modelstate->chain[i] = 0.5;
        }
    }
    
    if (modelstate->curve == CURVE_LINEAR) {
        for (i=0;i<modelstate->states;i++) {
            x = low+(stepsize*i);            
            modelstate->chain[i] = x;
        }
    }
    
    if (modelstate->curve == CURVE_AGEBRAIC) {
        for (i=0;i<modelstate->states;i++) {
            x = low+(stepsize*i);
            modelstate->chain[i] = x/sqrt(1.0+(x*x));
        }
    }    

    if (modelstate->curve == CURVE_ATAN) {
        for (i=0;i<modelstate->states;i++) {
            x = low+(stepsize*i);
            modelstate->chain[i] = atan(x);
        }
    } 

    if (modelstate->curve == CURVE_TANH) {
        for (i=0;i<modelstate->states;i++) {
            x = low+(stepsize*i);
            modelstate->chain[i] = tanh(x);
        }
    }
    
    if (modelstate->curve == CURVE_ERF) {
        for (i=0;i<modelstate->states;i++) {
            x = low+(stepsize*i);
            modelstate->chain[i] = erf(x);
        }
    }
    
    if (modelstate->curve == CURVE_GUDERMANN) {
        for (i=0;i<modelstate->states;i++) {
            x = low+(stepsize*i);
            modelstate->chain[i] = 2.0*atan(tanh(x/2.0));
        }
    }
    
    if (modelstate->curve == CURVE_LOGISTIC) {
        for (i=0;i<modelstate->states;i++) {
            x = low+(stepsize*i);
            modelstate->chain[i] = 1.0/(1.0+exp(-x));
        }
    }
    
    // Print out the Markov Chain
    if (verbose_mode==3) {
        printf("  MARKOV Chain\n");
        for (i=0;i<modelstate->states;i++) {
            fprintf(stderr,"  state=%03d PLeft=%0.4f\n",i,modelstate->chain[i]);  
        }
    }
    
    // Scale the curve to be between 0 and 1
    // First find the highest and lowest point.
    // Then sqish it.
    if (modelstate->curve != CURVE_FLAT) {
        themin = 10000.0;
        themax = -10000.0;
    
        for (i=0;i<modelstate->states;i++) {
            if (modelstate->chain[i] > themax) themax=modelstate->chain[i];
        }
    
        for (i=0;i<modelstate->states;i++) {
            if (modelstate->chain[i] < themin) themin=modelstate->chain[i];
        }
    
        h = themax - themin;    
    
        for (i=0;i<modelstate->states;i++) {
            x= modelstate->chain[i];
            x = (x-themin)/h;
            modelstate->chain[i] = x;
        }
    }
}


void sinbiasinit(t_modelstate *modelstate, t_rngstate *rngstate)
{
    unsigned char out[16];

    //smoothinit(modelstate, rngstate);
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

    modelstate->time = 0;
} 

void lcginit(t_modelstate *modelstate, t_rngstate *rngstate)
{
    if (rngstate->randseed==0) {
        modelstate->lcg_x = 0x63a3d28a2682b002ULL % modelstate->lcg_m;
        modelstate->lcg_index = 0;
    }
    else
    {
        nondeterministic_bytes(sizeof(unsigned long long), &(modelstate->lcg_x), rngstate);
        modelstate->lcg_x = modelstate->lcg_x % modelstate->lcg_m;
        modelstate->lcg_index = 0;
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
    //smoothinit(modelstate, rngstate);
}

void normalinit(t_modelstate *modelstate, t_rngstate *rngstate)
{
    unsigned char out[16];

    //smoothinit(modelstate, rngstate);
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

