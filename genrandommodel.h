
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
#include <stdlib.h>
#include <math.h>

typedef struct {
		double t;
		int lastbit;
		double left_stepsize;
		double right_stepsize;
		double sums_bias;
		double correlation;
		double bias;
		double mean;
		double variance;
		
		unsigned long long lcg_a; /* LCG Model */
		unsigned long long lcg_c;
		unsigned long long lcg_m;
		unsigned long long lcg_x;
		unsigned int lcg_truncate;
		unsigned int lcg_outbits;	
			
		int using_stepnoise;
		double stepnoise;
		int using_jfile;
		FILE *jfile;
		int using_infile;
		FILE *infile;
	} t_modelstate;

typedef struct {
		int c_max;
		int c;
		int input_format;
		unsigned char pool0[16];
		unsigned char pool1[16];
		unsigned char filechar;
		unsigned int  fileindex;
		unsigned char k[16];
		unsigned char kprime[16];
		unsigned char v[16];
		int temp;
		unsigned char rngbits[16];
		int randseed;
		int devurandom_available;
		int rdrand_available;
		FILE* devrandom;
		int reached_eof;
        int windowsrng;
	} t_rngstate;

int getrand16(t_rngstate* rngstate);
double getNormal(t_modelstate *modelstate, t_rngstate* rngstate);

/* entropysource(j,stepsize,v,seed); */
/* compute next bit */
/* j is the current step position */
/* stepsize is the stepsize */
/* k is a pointer to a 16 byte key for the simulation RNG */
/* v is the current CTR vector for the simulation RNG     */
/* It returns the next state for j */
int smoothsource(     t_modelstate* modelstate, t_rngstate* rngstate);
int puresource(       t_modelstate *modelstate, t_rngstate* rngstate);
int biasedsource(     t_modelstate *modelstate, t_rngstate* rngstate);
int correlatedsource( t_modelstate *modelstate, t_rngstate* rngstate);
int lcgsource( t_modelstate *modelstate, t_rngstate* rngstate);
int filesource(       t_modelstate *modelstate, t_rngstate* rngstate);
int filesourcehex(       t_modelstate *modelstate, t_rngstate* rngstate);

double normalsource(t_modelstate *modelstate, t_rngstate* rngstate);

/* initialize_sim(k,v,seed);                                      */
/* set aside two 16 byte vectors k and v to keep the simulation's */
/* PRNG state. Also create a 16 byte vector 'seed' and set it to  */
/* some value. The simulation is deterministic but will take a    */
/* different path for each value of seed.                         */
/* pass k,v and seed as pointers to unsigned char.                */
void smoothinit(     t_modelstate* modelstate, t_rngstate* rngstate);
void pureinit(       t_modelstate* modelstate, t_rngstate* rngstate);
void biasedinit(     t_modelstate* modelstate, t_rngstate* rngstate);
void correlatedinit( t_modelstate* modelstate, t_rngstate* rngstate);
void lcginit( t_modelstate* modelstate, t_rngstate* rngstate);
void normalinit(     t_modelstate *modelstate, t_rngstate* rngstate);
void fileinit(       t_modelstate *modelstate, t_rngstate* rngstate);

