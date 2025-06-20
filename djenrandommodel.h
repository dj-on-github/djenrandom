
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
#include <stdint.h>

#define CURVE_FLAT      0
#define CURVE_LINEAR    1
#define CURVE_LOGISTIC  2
#define CURVE_TANH      3
#define CURVE_ATAN      4
#define CURVE_GUDERMANN 5
#define CURVE_ERF       6
#define CURVE_AGEBRAIC  7

#define PUNCTURING_LEVEL_LOW  0
#define PUNCTURING_LEVEL_HIGH 1
#define PUNCTURING_LEVEL_BOTH 2

#define PUNC_STATE_STARTING  0
#define PUNC_STATE_GAPPING   1
#define PUNC_STATE_INJECTING 2
#define PUNC_STATE_FINISHED  3

#define PCG_LCG 1
#define PCG_MCG 2

#define XSH_RS 1
#define XSH_RR 2

typedef struct {
        double t;
        int lastbit;
        double left_stepsize;
        double right_stepsize;
        double sigmoid_bias;
        double correlation;
        double entropy;
        double averageentropy;
        int    n;
        double bias;
        double mean;
        double variance;
        int    states;
        int    sigmoid_state;
        double min_range;
        double max_range;
        double *chain;
        int    curve;
        char   curvestr[25];
        double p01;
        double p10;
        int    markov2p_phase;
        int    markov2p_symbol;
        int    *sampletable0;
        int    *sampletable1;
        int    p01_threshold;
        int    p10_threshold;
        int    bitwidth;
        int    gotentropy;
        int    gotp01;
        int    gotp10;
        int    gotbias;
        int    gotcorrelation;
        int    fast_m2p;        
        unsigned long long lcg_a; /* LCG Model */
        unsigned long long lcg_c;
        unsigned long long lcg_m;
        unsigned long long lcg_mask;
        unsigned long long lcg_x;
        unsigned int lcg_truncate;
        unsigned int lcg_outbits;
        unsigned int lcg_index;
        unsigned long long lcg_output;
        
        /* PCG States */
        unsigned int pcg_state_size;
        unsigned int pcg_index;
        int pcg_alg;
        int pcg_of;
        uint64_t pcg_output;
        
        uint16_t pcg16_state;
        uint32_t pcg32_state;
        uint64_t pcg64_state;
        uint64_t pcg128_state[2];
        
        uint16_t  pcg16_multiplier;
        uint32_t  pcg32_multiplier;
        uint64_t  pcg64_multiplier;
        uint64_t  pcg128_multiplier[2];
        
        uint16_t   pcg16_adder;
        uint32_t  pcg32_adder;
        uint64_t  pcg64_adder;
        uint64_t  pcg128_adder[2];
        
        /* sin bias states */
        uint64_t sinbias_period;
        double sinbias_amplitude;
        double sinbias_offset;
        int    sinbias_bias;
        uint64_t time;
    
        /* XORSHIFT States */
        
        uint32_t xorshift_size;
        uint32_t xorshift_state_a;
        uint32_t xorshift_state_b;
        uint32_t xorshift_state_c;
        uint32_t xorshift_state_d;
       
        /* Punctured States */

        uint32_t puncturing_level;
        uint32_t puncturing_length;
        uint32_t puncturing_start;
        uint32_t puncturing_gap;
        uint32_t punc_state;
        uint32_t punc_counter;
        uint32_t punc_both_level;
        uint32_t punc_event_count;
        uint32_t punc_event_limit;
        uint32_t punc_event_dis;

        /* General states */
        int using_stepnoise;
        double stepnoise;
        int using_jfile;
        FILE *jfile;
        int using_infile;
        FILE *infile;
        int using_json;
        int using_yaml;
        FILE *json_file;
        FILE *yaml_file;
        char filename[1000];
        int using_ofile;

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
        int got_detseed;
        unsigned char detseed[1024];
    } t_rngstate;

void init_rng(t_rngstate* rngstate);
int getrand16(t_rngstate* rngstate);
uint64_t getrand64(t_rngstate* rngstate);
double getNormal(t_modelstate *modelstate, t_rngstate* rngstate);
double get_rand_double(t_rngstate* rngstate);
uint64_t choose_exponent(uint64_t start, t_rngstate* rngstate);

/* entropysource(j,stepsize,v,seed); */
/* compute next bit */
/* j is the current step position */
/* stepsize is the stepsize */
/* k is a pointer to a 16 byte key for the simulation RNG */
/* v is the current CTR vector for the simulation RNG     */
/* It returns the next state for j */
int smoothsource(     t_modelstate* modelstate, t_rngstate* rngstate);
int puncturingsource(     t_modelstate* modelstate, t_rngstate* rngstate);
int puresource(       t_modelstate *modelstate, t_rngstate* rngstate);
int biasedsource(     t_modelstate *modelstate, t_rngstate* rngstate);
int correlatedsource( t_modelstate *modelstate, t_rngstate* rngstate);
int markov2psource( t_modelstate *modelstate, t_rngstate* rngstate);
int markov2pfastsource( t_modelstate *modelstate, t_rngstate* rngstate);
int markovsigmoidsource( t_modelstate *modelstate, t_rngstate* rngstate);
int sinbiassource( t_modelstate *modelstate, t_rngstate* rngstate);
int lcgsource( t_modelstate *modelstate, t_rngstate* rngstate);
int pcgsource( t_modelstate *modelstate, t_rngstate* rngstate);
int xorshiftsource( t_modelstate *modelstate, t_rngstate* rngstate);
int filesource(       t_modelstate *modelstate, t_rngstate* rngstate);
int filesourcehex(       t_modelstate *modelstate, t_rngstate* rngstate);
int filesourcebinary(t_modelstate* modelstate, t_rngstate* rngstate);
double normalsource(t_modelstate *modelstate, t_rngstate* rngstate);
int normalintegersource(t_modelstate *modelstate, t_rngstate* rngstate);

/* initialize_sim(k,v,seed);                                      */
/* set aside two 16 byte vectors k and v to keep the simulation's */
/* PRNG state. Also create a 16 byte vector 'seed' and set it to  */
/* some value. The simulation is deterministic but will take a    */
/* different path for each value of seed.                         */
/* pass k,v and seed as pointers to unsigned char.                */
void smoothinit(     t_modelstate* modelstate, t_rngstate* rngstate);
void puncturinginit(     t_modelstate* modelstate, t_rngstate* rngstate);
void pureinit(       t_modelstate* modelstate, t_rngstate* rngstate);
void biasedinit(     t_modelstate* modelstate, t_rngstate* rngstate);
void correlatedinit( t_modelstate* modelstate, t_rngstate* rngstate);
void markov2pinit(   t_modelstate* modelstate, t_rngstate* rngstate);
void markov2pfastinit(t_modelstate* modelstate, t_rngstate* rngstate);
void markovsigmoidinit(t_modelstate* modelstate, t_rngstate* rngstate);
void sinbiasinit(    t_modelstate* modelstate, t_rngstate* rngstate);
void lcginit(        t_modelstate* modelstate, t_rngstate* rngstate);
void pcginit(        t_modelstate* modelstate, t_rngstate* rngstate);
void xorshiftinit(   t_modelstate* modelstate, t_rngstate* rngstate);
void normalintegerinit(t_modelstate *modelstate, t_rngstate* rngstate);
void normalinit(     t_modelstate *modelstate, t_rngstate* rngstate);
void fileinit(       t_modelstate *modelstate, t_rngstate* rngstate);


