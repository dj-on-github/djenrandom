
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
#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>

#include "aes128k128d.h"
#include "djenrandommodel.h"
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include "rdrand_stdint/rdrand_stdint.h"
#include "markov2p.h"

#define MODEL_SUMS 0
#define MODEL_PURE 1
#define MODEL_BIASED 2
#define MODEL_CORRELATED 3
#define MODEL_NORMAL 4 
#define MODEL_FILE 5
#define MODEL_LCG 6
#define MODEL_PCG 7
#define MODEL_XORSHIFT 8
#define MODEL_SINBIAS 9
#define MODEL_MARKOV2P 10
#define MODEL_MARKOV_SIGMOID 11
#define MODEL_PUNCTURING 12
#define MODEL_NORMALINTEGER 13 
                    
#define INFORMAT_01 0
#define INFORMAT_HEX 1
#define INFORMAT_BINARY 2

#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"    
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

int aesni_supported;
int verbose_mode;
int puncture_count=0;
int puncture=0;
int puncture_state = PUNC_STATE_STARTING;
int reset_shiftreg = 1;  // reset xor feedback shiftregister every 512 bits

int puncture_this_bit(t_modelstate* modelstate);


void display_usage() {
fprintf(stderr,"Usage: djrandom [-bsvhn] [-x <bits>] [-y <bits>] [-z <bits>] [-c <generate length>]\n");
fprintf(stderr,"       [-m <|pure(default)|sums|biased|correlated|normal|normalinteger|sinbias|markov_2_param|puncturing|file>] [-l <left_stepsize>]\n"); 
fprintf(stderr,"       [-r <right_stepsize>] [--stepnoise=<noise on step>] [--bias=<bias>]\n");
fprintf(stderr,"       [--correlation=<correlation>] [--mean=<normal mean>] [--variance=<normal variance>]\n");
fprintf(stderr,"       [--pcg_state_16=<16|32|64>] [--pcg_generator=<LCG|MCG>] [--pcg_of=<XSH_RS|XSH|RR]\n");
fprintf(stderr,"       [--sinbias_offset=<0.0 to 1.0>] [--sinbias_amplitude=<0.0 to 1.0>] [--sinbias_period=<samples per cycle>]\n");
fprintf(stderr,"       [--p10=<probability of 10 transition] [--p01=<probability of 01 transition>]\n");
fprintf(stderr,"       [--states=<integer of number of states in the markov chain>]\n");
fprintf(stderr,"       [--sigmoid=<flat|linear|sums|logistic|tanh|atan|gudermann|erf|algebraic]\n");
fprintf(stderr,"       [--min_range=<float less than max_range>][--max_range=<float greater than min_range>]\n");
fprintf(stderr,"       [-o <output_filename>] [-j <j filename>] [-i <input filename>] [-f <hex|binary|01>]\n");
fprintf(stderr,"       [-J <json_filename>] [-Y <yaml_filename>]\n");
fprintf(stderr,"       [--bpb=<binary bits per byte>]\n");
fprintf(stderr,"       [-k <1K_Blocks>] [-w [1..256]]\n");
fprintf(stderr,"       [-D <deterministic seed string>]\n");
fprintf(stderr,"\n");
fprintf(stderr,"Generate random bits with configurable non-uniformities.\n");
fprintf(stderr,"  Author: David Johnston, dj@deadhat.com\n");
fprintf(stderr,"\n");

fprintf(stderr,"  -m, --model=<pure(default)|sums|biased|correlated|lcg|pcg|xorshift|normal|normalinteger|file>\n");
fprintf(stderr,"              Select random source model\n");

fprintf(stderr,"\nStep Update Metastable Source model (-m sums) Options\n\n");
fprintf(stderr,"  -l, --left=<left_stepsize>     stepsize when moving left as a fraction of sigma_m.\n");
fprintf(stderr,"  -r, --right=<right_stepsize>   stepsize when moving right as a fraction of sigma_m.\n");
fprintf(stderr,"  --stepnoise=<noise on step>    variance of the noise on stepsize. e.g. 0.00001.\n");

fprintf(stderr,"\nBiased model (-m biased) Options\n\n");
fprintf(stderr,"  --bias=<bias>                  bias as a number between 0.0 and 1.0.\n");
fprintf(stderr,"                                 Only for biased or markov model\n");
fprintf(stderr,"\nCorrelated model (-m correlated) Options\n\n");
fprintf(stderr,"  --correlation=<correlation>    Correlation with previous bit as a number between -1.0 and 1.0.\n");
fprintf(stderr,"                                 Only for correlation or markov model\n");

fprintf(stderr,"\nSinusoidally Varying Bias model (-m sinbias) Options\n\n");
fprintf(stderr,"  --sinbias_amplitude=<0.0 to 1.0>     Amplitude of the variation of the bias between 0.0 and 1.\n");
fprintf(stderr,"                                      Only for sinbias model\n");
fprintf(stderr,"  --sinbias_offset=<0.0 to 1.0>        Midpoint Offset of the varying bias between 0.0 and 1.0.\n");
fprintf(stderr,"                                       Only for sinbias model\n");
fprintf(stderr,"  --sinbias_period=<samples per cycle> Number of samples for a full cycle of the sinusoidally\n");
fprintf(stderr,"                                       varying bias. Only for sinbias model\n");

fprintf(stderr,"\nTwo Parameter Markov model (-m markov_2_param) Options\n\n");
fprintf(stderr,"  --fast                    Use a fast version on the generator.\n");
fprintf(stderr,"         and one set of:\n");
fprintf(stderr,"  --p10=<0.0 to 1.0>        The probability of a 1 following a 0, default 0.5\n");
fprintf(stderr,"  --p01=<0.0 to 1.0>        The probability of a 0 following a 1, default 0.5\n");
fprintf(stderr,"         or\n");
fprintf(stderr,"  --bias=<0.0 to 1.0>               The ones probability, default 0.5\n");
fprintf(stderr,"  --correlation=<-1.0 to 1.0>       The serial correlation coefficient, default 0.0\n");
fprintf(stderr,"         or\n");
fprintf(stderr,"  --entropy=<0.0 to 1.0>    The per bit entropy, default 1.0\n");
fprintf(stderr,"  --bitwidth=<3 to 64>      The number of bits per symbol\n");

fprintf(stderr,"\nSigmoid Markov model (-m markov_sigmoid) Options\n\n");
fprintf(stderr,"  --states=<n>              The number of states in the Markov Chain\n");
fprintf(stderr,"  --sigmoid=<curve>         Curve name, one of: flat, linear, sums, logistic, tah, atan,\n");
fprintf(stderr,"                            gudermann, erf or algebraic, default linear\n");
fprintf(stderr,"  --min_range=<float>       The start of the range of the curve. Usually between -5.0 and -2.0\n");
fprintf(stderr,"  --max_range=<float>       The end of the range of the curve. Usually between 2.0 and 5.0\n");

fprintf(stderr,"\nNormal model (-m normal) Options\n\n");
fprintf(stderr,"  --mean=<normal mean>           mean of the normally distributed data. Only for normal model\n");
fprintf(stderr,"  --variance=<normal variance>   variance of the normally distributed data\n");

fprintf(stderr,"\nNormal Integer model (-m normalinteger) Options\n\n");
fprintf(stderr,"    Pass in floating point mean and variance, but set the mean to the middle of where you want it in the 32 bit integer space.\n");
fprintf(stderr,"    A floating point Gaussian generator is used and then sampled into the int32 values.\n");
fprintf(stderr,"  --mean=<normal mean>           mean of the normally distributed data. Only for normal model\n");
fprintf(stderr,"  --variance=<normal variance>   variance of the normally distributed data\n");

fprintf(stderr,"\nLinear Congruential Generator model (-m lcg) Options\n\n");
fprintf(stderr,"  --lcg_a=<LCG multipler term>  Positive integer less than lcg_m\n");
fprintf(stderr,"  --lcg_c=<LCG additive term>   Positive integer less than lcg_m\n");
fprintf(stderr,"  --lcg_m=<LCG modulo term>     Positive integer defining size of the group\n");
fprintf(stderr,"  --lcg_truncate=<lower bits to truncate>     Positive integer\n");
fprintf(stderr,"  --lcg_outbits=<Number of bits per output>     Positive integer\n");

fprintf(stderr,"\nPermuted Congruential Generator model (-m pcg) Options\n\n");
fprintf(stderr,"  --pcg_state_size=<state size of PCG>  16 ,32 or 64\n");
fprintf(stderr,"  --pcg_generator=<Generator Algorithm> MCG or LCG\n");
fprintf(stderr,"  --pcg_of=<Output Function>            XSH_RS or XSH_RR\n");

fprintf(stderr,"\nXorShift model (-m xorshift) Options\n\n");
fprintf(stderr,"  --xorshift_size=[state size of xorshift]  32 or 128\n");

fprintf(stderr,"\nPuncturing model (-m puncturing) or other models (--puncture) Options\n\n");
fprintf(stderr,"  Applies to uniform data with -m puncturing or any other model with --puncture\n");
fprintf(stderr,"  --puncturing_start=[bit position to start injecting constant data]\n");
fprintf(stderr,"  --puncturing_length=[length in bits of injected data]\n");
fprintf(stderr,"  --puncturing_level=[value to inject] 0, 1 or 2=randomize between high or low\n");
fprintf(stderr,"  --puncturing_gap=[length in bits of the gap between injection points]\n");
fprintf(stderr,"  --puncturing_limit=[Maximum number of puncturing events to perform. 0 for no limit.]\n");

fprintf(stderr,"\nGeneral Options\n\n");
fprintf(stderr,"  -x, --xor=<bits>               XOR 'bits' of entropy together for each output bit\n");
fprintf(stderr,"  -y, --xmin=<bits>              Provides the start of a range of XOR ratios\n");
fprintf(stderr,"                                 to be chosen at random per sample\n");
fprintf(stderr,"  -z, --xmax=<bits>              Provides the end of a range of XOR ratios to be\n");
fprintf(stderr,"                                 chosen at random per sample\n");
fprintf(stderr,"  --xor4                         Enable 4 bit xor feedback digitizer\n");
fprintf(stderr,"  --xor11                        Enable 11 bit xor feedback digitizer\n");
fprintf(stderr,"  --noresetxor                   xor feedback is reset every 512 bits. This disabled that.\n");
fprintf(stderr,"  -s, --seed                     Nondeterministically seed the internal RNG with /dev/random\n");
fprintf(stderr,"  -D, --detseed <seed string>    Deterministically seed the internal RNG with the given string\n");
fprintf(stderr,"  -n, --noaesni                  Don't use AESNI instruction.\n");
fprintf(stderr,"  -c, --cmax=<generate length>   number of PRNG generates before a reseed\n");
fprintf(stderr,"  -v, --verbose                  output the parameters\n");

fprintf(stderr,"\nFile Options\n\n");
fprintf(stderr,"  -o <output_filename>             output file\n");
fprintf(stderr,"  -j, --jfile=<j filename>         filename to push source model internal state to\n");
fprintf(stderr,"  -i, --infile=<input filename>    filename of entropy file for file model\n");
fprintf(stderr,"  -f, --informat=<hex|binary|01>   Format of input file.:\n");
fprintf(stderr,"                                       hex=Ascii hex(default),\n");
fprintf(stderr,"                                       4 bit per hex character.\n");
fprintf(stderr,"                                       binary=raw binary.\n");
fprintf(stderr,"                                       01=ascii binary.\n");
fprintf(stderr,"                                       Non valid characters are ignored\n");
fprintf(stderr,"  -J, --json=<JSON filename>       filename to output JSON information of the data to\n");
fprintf(stderr,"  -Y, --yaml=<YAML filename>       filename to output YAML information of the data to\n");
fprintf(stderr,"  -k, --blocks=<1K_Blocks>         Size of output in kilobytes\n");

fprintf(stderr,"\nOutput Format Options\n\n");
fprintf(stderr,"  -b, --binary                output in raw binary format\n");
fprintf(stderr,"  --nistoddball               output in nistoddball format, 1 bit per byte\n");
fprintf(stderr,"  --bpb                       Number of bits per byte to output in binary output mode.\n");
fprintf(stderr,"                              Default 8 bits.\n");
fprintf(stderr,"  -w, --width=[1...256]       Byte per line of output\n");

fprintf(stderr,"\nThe most important option of all\n\n");
fprintf(stderr,"  -h, --help                     print this help and exit\n");
}

void printsample(unsigned char *thesample)
{
    int tempindex;
    int j;
    int i;
   tempindex = 0;
    for (j=0;j<16;j++)
    {
            for (i=0;i<16;i++) printf("%02X",thesample[tempindex++]);
            printf("\n");
    }
}

double logtwo(double x)
{
    double result;
        //result = log(x)/log(2);
        result = log2l(x);
        return(result);
}

/* Return the entropy in a single biased bit */
double bias2entropy(double bias)
{
    double result;

    if (bias == 0.0) result = 0.0;
    else if (bias == 1.0) result = 0.0;
    else
    {
        result = -(bias)*logtwo(bias) -((1.0L-bias)*logtwo(1.0L-bias));
    
        if (isnan(result)) fprintf(stderr," Got NAN in bias2entropy bias=%f\n",bias);
    }
    return(result);
}

/* return bias on two biased inputs to an xor */
double xorbias_2bit(double pa1, double pb1)
{
    double result;
    double pa0;
    double pb0;

    pa0 = 1.0-pa1;
    pb0 = 1.0-pb1;

    result = (pa0 * pb1) + (pa1 * pb0);

    if (isnan(result)) fprintf(stderr," Got NAN in xorbias_2bit\n");
    return(result);
}

/* return bias on three biased inputs to an xor */
double xorbias_3bit(double pa1, double pb1, double pc1)
{
    double result;
    double pa0;
    double pb0;
    double pc0;

    pa0 = 1.0-pa1;
    pb0 = 1.0-pb1;
    pc0 = 1.0-pc1;

    /* Sum the joint probabilities of all the odd parity patterns. */
    result = (pa0 * pb0 * pc1) + (pa0 * pb1 * pc0) + (pa1 * pb0 * pc0) + (pa1 * pb1 * pc1);

    if (isnan(result)) fprintf(stderr," Got NAN in xorbias_3bit\n");
    return(result);
}

/* Select between the different
 * entropy source models.
 */

double dbl_entropysource(int model, t_modelstate* modelstate, t_rngstate* rngstate)
{
    double result;
    
    if (model==MODEL_NORMAL)
    {
        result = normalsource(modelstate, rngstate);
        return result;
    }
    else return(0.0);
}

int entropysource(int model, t_modelstate* modelstate, t_rngstate* rngstate)
{
    int result;
    if (model==MODEL_NORMALINTEGER)
    {
        result = normalintegersource(modelstate, rngstate);
    }
    else if (model==MODEL_SUMS)
    {
        result = smoothsource(modelstate, rngstate);
    }
    else if (model==MODEL_PURE)
    {
        result = puresource(modelstate, rngstate);
    }
    else if (model==MODEL_PUNCTURING)
    {
        result = puncturingsource(modelstate, rngstate);
    }
    else if (model==MODEL_BIASED)
    {
        result = biasedsource(modelstate, rngstate);
        return result;
    }
    else if (model==MODEL_CORRELATED)
    {
        result = correlatedsource(modelstate, rngstate);
    }
    else if (model==MODEL_MARKOV2P)
    {
        if (modelstate->fast_m2p==1) {
            result = markov2pfastsource(modelstate, rngstate);
        } else {
            result = markov2psource(modelstate, rngstate);
        }
    }
    else if (model==MODEL_MARKOV_SIGMOID)
    {
        result = markovsigmoidsource(modelstate, rngstate);
    }
    else if (model==MODEL_SINBIAS)
    {
        result = sinbiassource(modelstate, rngstate);
    }
    else if (model==MODEL_LCG)
    {
        result = lcgsource(modelstate, rngstate);
    }
    else if (model==MODEL_PCG)
    {
        result = pcgsource(modelstate, rngstate);
    }
    else if (model==MODEL_XORSHIFT)
    {
        result = xorshiftsource(modelstate, rngstate);
    }
    else if (model==MODEL_FILE)
    {
        if (rngstate->input_format==INFORMAT_HEX)
            result = filesourcehex(modelstate,rngstate);
        else if (rngstate->input_format==INFORMAT_01)
            result = filesource(modelstate,rngstate);
        else
            result = filesourcebinary(modelstate, rngstate);
    }
    else {
        result = 0;
    }
    
    // Apply the puncturing

    if (puncture==1) {
        if (puncture_this_bit(modelstate) == 1){
            if (modelstate->puncturing_level==0) result = 0;
            if (modelstate->puncturing_level==1) result = 1;
            if (modelstate->puncturing_level==2) result = 0; // TBD have not implemented this yet
        }
    }

    return result;
}

int puncture_this_bit(t_modelstate* modelstate) {
    if (puncture_state == PUNC_STATE_STARTING) {
        if (puncture_count < (modelstate->puncturing_start-1)) {
            puncture_count++;
            //fprintf(stderr,"STARTING %d, start=%d\n",puncture_count,(modelstate->puncturing_start-1));
            return 0;
        } else {
            modelstate->punc_event_count = 0;
            puncture_state = PUNC_STATE_INJECTING;
            puncture_count=0;
            return 0;
        }
    } else if (puncture_state == PUNC_STATE_INJECTING) {
        if (puncture_count < (modelstate->puncturing_length-1)) {
            puncture_count++;
            //fprintf(stderr,"INJECTING %d\n",puncture_count);
            return 1;
        } else {
            modelstate->punc_event_count++;
            puncture_state = PUNC_STATE_GAPPING;
            puncture_count=0;
            return 1;
        }
        
    } else if (puncture_state == PUNC_STATE_GAPPING) {
        if (puncture_count < (modelstate->puncturing_gap-1)) {
            puncture_count++;
            //fprintf(stderr,"GAPPING %d\n",puncture_count);
            return 0;
        } else if ((modelstate->punc_event_count < modelstate->punc_event_limit) || (modelstate->punc_event_dis==1)) {
            puncture_state = PUNC_STATE_INJECTING;
            puncture_count=0;
            return 0;
        } else {
            puncture_state = PUNC_STATE_FINISHED;
            return 0;
        }
    } else if (puncture_state == PUNC_STATE_FINISHED) {
        puncture_state = PUNC_STATE_FINISHED;
        return 0;
    } else {
        fprintf(stderr,"ERROR : puncture_state wrong %d\n",puncture_state);
        exit(1);
        return 0;
    }

}

void initialize_sim(int model, t_modelstate* modelstate, t_rngstate* rngstate)
{
    if (model==MODEL_SUMS)
    {
        smoothinit(modelstate, rngstate);
    }
    else if (model==MODEL_PURE)
    {
        pureinit(modelstate, rngstate);
    }
    else if (model==MODEL_PUNCTURING)
    {
        puncturinginit(modelstate, rngstate);
    }
    else if (model==MODEL_BIASED)
    {
        biasedinit(modelstate, rngstate);
    }
    else if (model==MODEL_CORRELATED)
    {
        correlatedinit(modelstate, rngstate);
    }
    else if (model==MODEL_MARKOV2P)
    {
        if (modelstate->fast_m2p==1) {
            markov2pfastinit(modelstate, rngstate);
        } else {
            markov2pinit(modelstate, rngstate);
        }
    }
    else if (model==MODEL_MARKOV_SIGMOID)
    {
        markovsigmoidinit(modelstate, rngstate);
    }
    else if (model==MODEL_SINBIAS)
    {
        sinbiasinit(modelstate, rngstate);
    }
    else if (model==MODEL_NORMAL)
    {
        normalinit(modelstate, rngstate);
    }
    else if (model==MODEL_NORMALINTEGER)
    {
        normalintegerinit(modelstate, rngstate);
    }
    else if (model==MODEL_LCG)
    {
        lcginit(modelstate, rngstate);
    }
    else if (model==MODEL_PCG)
    {
        pcginit(modelstate, rngstate);
    }
    else if (model==MODEL_XORSHIFT)
    {
        xorshiftinit(modelstate, rngstate);
    }
    else if (model==MODEL_FILE)
    {
        fileinit(modelstate, rngstate);
    }
}

/********
* main() is mostly about parsing and qualifying the command line options. The argtable2 linrary is used to help with this.
*/

int main(int argc, char** argv)
{
    int opt;
    int i;
    int j;
    unsigned char abyte;
    int onek;
    int tempindex;
    int xor4bit;
    int xorcount;
    int xor11bit;
    int downsample;
    int noresetxor=0;

    //int inputbits;
    int shiftreg;
    int newbit;
    int xoriter;
    int abort;
    double thevalue;
    int thebit;
    unsigned char thebyte;
    int binary_mode;
    int hex_mode;
    int nistoddball_mode;
    int bits_per_byte;
    int xormode;
    int xorbits;
    //int verbose_mode;
    int kilobytes;
    int no_k=1;
    int ofile;
    char errstr[2000];
    char filename[1000];
    char jfilename[1000];
    char json_filename[1000];
    char yaml_filename[1000];
    char infilename[1000];
    int model;
    int using_xor_range;
    int xmin;
    int xmax;
    
    int input_format;
    int linewidth;
    int lineindex;
    int width;
   
    /* unsigned char entropy;*/
    //double sums_entropy;
    double postxor_entropy;
    double total_entropy;
    double prob;
    
    int samplenum;
    int simrun;
    unsigned char thesample[256];
    unsigned char thebpbsample[2048];
    double floatingpointsamples[256];
    int32_t normalintegersamples[32];
    t_rngstate rngstate;
    t_modelstate modelstate;

    int gotcorrelation;
    int gotbias;
    //int gotmean;
    int gotp01;
    int gotp10;
    int gotentropy;
    //int gotbitwidth;
    FILE *fp=NULL;

    //double epsilon;
    
    /* Defaults */
    xor4bit = 0;    // 0=Do not use the 4 bit xor decorrelator, 1 = do.
    xorcount = 0;
    xor11bit = 0;    // 0=Do not use the 11 bit xor decorrelator, 1 = do.
    downsample = 0;    // 0=Do not use the downsampler of the decorellator, 1 = do.
    shiftreg = 0;  // Preset the decorrelator shiftregister to 0.
    puncture = 0; // Default to not puncture the source data

    binary_mode = 0; /* binary when 1 */
    nistoddball_mode = 0;
    hex_mode = 1;

    bits_per_byte = 8; /* default 8 bits per byte */
    ofile = 0;       /* use stdout instead of outputfile*/
    model = MODEL_PURE;
    xormode = 0;  /* do xor when 1, else don't do xor */
    xorbits = 3;  /* the number of bits to xor together when xormode=1 */
    verbose_mode = 0;
    kilobytes = 1;
    modelstate.using_jfile = 0;
    modelstate.using_infile = 0;
    modelstate.using_ofile = 0;
    modelstate.using_yaml = 0;
    modelstate.using_json = 0;
    input_format = INFORMAT_HEX;
    rngstate.got_detseed=0;
    rngstate.detseed[0]=0;
    rngstate.input_format = INFORMAT_HEX;
    rngstate.randseed=0;
    rngstate.rdrand_available=0;
    rngstate.devurandom_available=0;
    using_xor_range=0;
    xmin=0;
    xmax=0;
    aesni_supported = 0;

    modelstate.fast_m2p=0;  
    gotcorrelation = 0;
    gotbias = 0;
    //gotmean = 0;
    gotp01 = 0;
    gotp10 = 0;
    gotentropy = 0;
    //gotbitwidth = 0;
    
    modelstate.lcg_a = 0x05DEECE66DULL;  /* Posix RAND48 default */
    modelstate.lcg_c = 11ULL;
    modelstate.lcg_m = 0x0001000000000000ULL; /* 2**48 */
    modelstate.lcg_truncate = 15; /* Lower bits to truncate */
    modelstate.lcg_x = 0x63a3d28a2682b002ULL; /* Initial state */
    modelstate.lcg_outbits= 33; /* number of bits in output */
    modelstate.lcg_index = 0;

    modelstate.pcg_state_size=32;
    modelstate.pcg_index=0;
    modelstate.pcg_alg = PCG_LCG;
    modelstate.pcg_of = XSH_RR;  
    
    modelstate.sinbias_amplitude = 0.5;
    modelstate.sinbias_offset = 0.5;
    modelstate.sinbias_period = 1000;
    
    modelstate.time = 0;
    
    rngstate.c_max=511;
    
    modelstate.left_stepsize = 0.1;
    modelstate.right_stepsize = 0.1;
    modelstate.bias = 0.5;
    modelstate.correlation = 0.1;
    modelstate.mean = 0.0;
    modelstate.variance = 0.1;
    modelstate.using_stepnoise = 0;
    modelstate.stepnoise = 0.0;
    modelstate.p01 = 0.5;
    modelstate.p10 = 0.5;
    modelstate.bitwidth=8;
    
    modelstate.curve = CURVE_LINEAR;
    strcpy(modelstate.curvestr,"Linear");
    modelstate.states = 21;
    modelstate.sigmoid_state = 10;
    modelstate.min_range=-1.0;
    modelstate.max_range=1.0;
    modelstate.xorshift_size=32;

    modelstate.puncturing_level = PUNCTURING_LEVEL_LOW;
    modelstate.puncturing_start = 0;      // start injecting on first bit
    modelstate.puncturing_length = 167;   // ceil(100/H) , H=0.6
    modelstate.puncturing_gap = 1024-modelstate.puncturing_length; // Inject 167 bits every 1024 bits
    modelstate.punc_event_count = 0;
    modelstate.punc_event_limit = 0;
    modelstate.punc_event_dis = 1;

    //sums_entropy = 0.0;
    postxor_entropy = 0.0;
    total_entropy = 0.0;
    width = 32;
    linewidth = 32;
    
    filename[0] = (char)0;
    modelstate.filename[0] = (char)0;
    jfilename[0] = (char)0;
    infilename[0] = (char)0;
    json_filename[0] = (char)0;
    yaml_filename[0] = (char)0;

    fp = NULL;
    modelstate.jfile = NULL;
    modelstate.infile = NULL;
    modelstate.yaml_file = NULL;
    modelstate.json_file = NULL;

    aesni_supported = aesni_check_support();

    int tempa;
    int tempb;
    int xorselector;
    int xorrange;
    double pa1;
    double pb1;
    double pc1;

    /* get the options and arguments */
    int longIndex;
    int gotxmin;
    int gotxmax;

    /*
    fprintf(stderr,"  ARGC=%d\n",argc);
    for (i=0;i<argc;i++) {
        fprintf(stderr,"  ARGV[%d] = %s  ",i,argv[i]);
        for(j=0;j<strlen(argv[i]);j++) fprintf(stderr,"%02x ",(unsigned char)(argv[i][j]));
        fprintf(stderr,"\n");
    }  
    */

    gotxmin = 0;
    gotxmax = 0;
    char optString[] = "c:m:l:r:B:C:o:j:J:Y:i:f:k:V:w:D:bxpsnvh";
    static const struct option longOpts[] = {
    { "binary", no_argument, NULL, 'b' },
    { "nistoddball", no_argument, NULL, 0 },
    { "hex", no_argument, NULL, 0 },
    { "p01", required_argument, NULL, 0 },
    { "p10", required_argument, NULL, 0 },
    { "bitwidth", required_argument, NULL, 0 },
    { "entropy", required_argument, NULL, 0 },
    { "fast", no_argument, NULL, 0 },
    { "bpb", required_argument, NULL, 0 },
    { "xor4bit", no_argument, NULL, 0 },
    { "xor11bit", no_argument, NULL, 0 },
    { "downsample", no_argument, NULL, 0 },
    { "noresetxor", no_argument, NULL, 0 },
    { "xor", required_argument, NULL, 'x' },
    { "xmin", required_argument, NULL, 0 },
    { "xmax", required_argument, NULL, 0 },
    { "seed", no_argument, NULL, 's' },
    { "detseed", required_argument, NULL, 'D' },
    { "noaesni", no_argument, NULL, 'n' },
    { "cmax", required_argument, NULL, 'c' },
    { "model", required_argument, NULL, 'm' },
    { "verbose", no_argument, NULL, 'v' },
    { "verbose_level", required_argument, NULL, 'V' },
    { "left", required_argument, NULL, 'l' },
    { "right", required_argument, NULL, 'r' },
    { "stepnoise", required_argument, NULL, 0 },
    { "bias", required_argument, NULL, 'B' },
    { "correlation", required_argument, NULL, 'C' },
    { "mean", required_argument, NULL, 0 },
    { "variance", required_argument, NULL, 0 },

    { "states", required_argument, NULL, 0 },
    { "sigmoid", required_argument, NULL, 0 },
    { "min_range", required_argument, NULL, 0 },
    { "max_range", required_argument, NULL, 0 },
           
    { "puncture", no_argument, NULL, 'p' },
    { "puncturing_start", required_argument, NULL, 0},
    { "puncturing_length", required_argument, NULL, 0},
    { "puncturing_gap", required_argument, NULL, 0},
    { "puncturing_level", required_argument, NULL, 0},
    { "puncturing_event_count", required_argument, NULL, 0},
    
    { "lcg_a", required_argument, NULL, 0 },
    { "lcg_c", required_argument, NULL, 0 },
    { "lcg_m", required_argument, NULL, 0 },
    { "lcg_truncate", required_argument, NULL, 0 },
    { "lcg_outbits", required_argument, NULL, 0 },
    
    { "pcg_state_size", required_argument, NULL, 0 },
    { "pcg_generator", required_argument, NULL, 0 },
    { "pcg_of", required_argument, NULL, 0 },
    
    { "sinbias_amplitude", required_argument, NULL, 0 },
    { "sinbias_offset", required_argument, NULL, 0 },
    { "sinbias_period", required_argument, NULL, 0 },
    
    { "xorshift_size", required_argument, NULL, 0 },
    
    { "output", required_argument, NULL, 'o' },
    { "jfile", required_argument, NULL, 'j' },
    { "json", required_argument, NULL, 'J' },
    { "yaml", required_argument, NULL, 'Y' },
    { "infile", required_argument, NULL, 'i' },
    { "informat", required_argument, NULL, 'f' },
    { "blocks", required_argument, NULL, 'k' },
    { "width", required_argument, NULL, 'w' },
    { "help", no_argument, NULL, 'h' },
    { NULL, no_argument, NULL, 0 }
    };

    /* Test for nondeterministic random source */
    
    //if (rngstate.randseed==1) {
        if (rdrand_check_support()==1) {    
            rngstate.rdrand_available=1;    
        }
        else if ((rngstate.devrandom =  fopen("/dev/urandom", "r")) !=NULL) {
            rngstate.devurandom_available = 1;
        }
        else {
            rngstate.rdrand_available=0;
            rngstate.devurandom_available = 0;
            //fprintf(stderr,"Neither /dev/urandom nor RdRand Supported for nondeterministic seeding.");
            //exit(1);
        }
    //}
    
    
    opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
    while( opt != -1 ) {
        //fprintf(stderr,"OPT = %c\n",opt);
        switch( opt ) {
            case 'b':
                binary_mode = 1;
                hex_mode = 0;
                nistoddball_mode = 0;
                break;
                
            case 'x':
                xormode = 1;
                xorbits = atoi(optarg);
                break;
                
            case 's':
                rngstate.randseed = 1;
                break;
            
            case 'D':
                if (strlen(optarg) < 1024) {
                    strcpy((char *)rngstate.detseed, optarg);
                    
                } else {
                    fprintf(stderr, "Error, deterministic seed cannot be longer than 1024 characters\n");
                    exit(1);
                }
                rngstate.got_detseed = 1;
                break;
                
            case 'n':
                aesni_supported = 0;
                break;
            
            case 'v':
                verbose_mode = 1;
                break;
            
            case 'V':
                verbose_mode = atoi(optarg);
                break;
            
            case 'l':
                modelstate.left_stepsize = atof(optarg);
                break;
               
            case 'B': 
                modelstate.bias = atof(optarg);
                gotbias=1;
                modelstate.gotbias=1;
                break;

            case 'C': 
                modelstate.correlation = atof(optarg);
                gotcorrelation=1;
                modelstate.gotcorrelation=1;
                break;

            case 'r':
                modelstate.right_stepsize = atof(optarg);
                break;
                
            case 'o':
                ofile = 1;
                strcpy(filename,optarg);
                strcpy(modelstate.filename,filename);
                modelstate.using_ofile=1;
                break;
                
            case 'j':
                modelstate.using_jfile = 1;
                strcpy(jfilename,optarg);
                break;
                
            case 'J':
                modelstate.using_json = 1;
                strcpy(json_filename,optarg);
                break;

            case 'Y':
                modelstate.using_yaml = 1;
                strcpy(yaml_filename,optarg);
                break;
            
            case 'i':
                modelstate.using_infile = 1;
                strcpy(infilename,optarg);
                break;
                
            case 'f':
                if (strcmp(optarg,"binary")==0) input_format=INFORMAT_BINARY;
                else if (strcmp(optarg,"hex")==0) input_format=INFORMAT_HEX;
                else if (strcmp(optarg,"01")==0) input_format=INFORMAT_01;
                else
                {
                    fprintf(stderr,"Input file format %s not recognized. Choose from binary, hex or 01.\n",optarg);
                    exit(1);
                }
                rngstate.input_format = input_format;
                break; 

            case 'm':
                if (strcmp(optarg,"sums")==0) model=MODEL_SUMS;
                else if (strcmp(optarg,"pure")==0) model=MODEL_PURE;
                else if (strcmp(optarg,"puncturing")==0) model=MODEL_PUNCTURING;
                else if (strcmp(optarg,"biased")==0) model=MODEL_BIASED;
                else if (strcmp(optarg,"correlated")==0) model=MODEL_CORRELATED;
                else if (strcmp(optarg,"markov_2_param")==0) model=MODEL_MARKOV2P;
                else if (strcmp(optarg,"markov_sigmoid")==0) model=MODEL_MARKOV_SIGMOID;
                else if (strcmp(optarg,"sinbias")==0) model=MODEL_SINBIAS;
                else if (strcmp(optarg,"lcg")==0) model=MODEL_LCG;
                else if (strcmp(optarg,"pcg")==0) model=MODEL_PCG;
                else if (strcmp(optarg,"xorshift")==0) model=MODEL_XORSHIFT;
                else if (strcmp(optarg,"normal")==0) model=MODEL_NORMAL;
                else if (strcmp(optarg,"normalinteger")==0) model=MODEL_NORMALINTEGER;
                else if (strcmp(optarg,"file")==0) model=MODEL_FILE;
                else
                {
                    fprintf(stderr,"model type %s not recognized. Choose from sums, pure, puncturing, biased, correlated, markov_2_param, markov_sigmoid, normal, normalinteger or file.\n",optarg);
                    exit(1);
                }
                break; 

                
            case 'k':
                kilobytes = atoi(optarg);
                no_k = 0;
                //fprintf(stderr,"atoi(optarg) = %d,  optarg = %s",kilobytes,optarg); 
                break;
                
            case 'w':
                width = atoi(optarg);
                if ((width > 0) && (width < 257)) {
                    linewidth=width;
                    //if (width == 8) linewidth=1;
                    //else if (width == 16) linewidth=2;
                    //else if (width == 32) linewidth=4;
                    //else if (width == 64) linewidth=8;
                    //else if (width == 128) linewidth=16;
                    //else if (width == 256) linewidth=32;
                }
                break;
            case 'p':
                //fprintf(stderr,"SETTING puncture=1 from options\n");
                puncture = 1;
                break;                
               
            case 0:     /* long option without a short arg */
                //fprintf(stderr," LONGOPT = %s\n",longOpts[longIndex].name);
                if( strcmp( "hex", longOpts[longIndex].name ) == 0 ) {
                    hex_mode = 1;
                    binary_mode = 0;
                    nistoddball_mode = 0;
                }
                if( strcmp( "binary", longOpts[longIndex].name ) == 0 ) {
                    hex_mode = 0;
                    binary_mode = 1;
                    nistoddball_mode = 0;
                }
                if( strcmp( "puncturing_start", longOpts[longIndex].name ) == 0 ) {
                    modelstate.puncturing_start = atoi(optarg);
                    if (modelstate.puncturing_start == 0) puncture_state = PUNC_STATE_INJECTING;
                }
                if( strcmp( "puncturing_length", longOpts[longIndex].name ) == 0 ) {
                    modelstate.puncturing_length = atoi(optarg);
                }
                if( strcmp( "puncturing_gap", longOpts[longIndex].name ) == 0 ) {
                    modelstate.puncturing_gap = atoi(optarg);
                }
                if( strcmp( "puncturing_level", longOpts[longIndex].name ) == 0 ) {
                    modelstate.puncturing_level = atoi(optarg);
                }
                if( strcmp( "puncturing_event_count", longOpts[longIndex].name ) == 0 ) {
                    modelstate.punc_event_limit = atoi(optarg);
                }
                if( strcmp( "nistoddball", longOpts[longIndex].name ) == 0 ) {
                    hex_mode = 0;
                    binary_mode = 0;
                    nistoddball_mode = 1;
                }
                if( strcmp( "bpb", longOpts[longIndex].name ) == 0 ) {
                    bits_per_byte = atoi(optarg);
                }

                if( strcmp( "xor4bit", longOpts[longIndex].name ) == 0 ) {
                    xor4bit = 1;
                }

                if( strcmp( "xor11bit", longOpts[longIndex].name ) == 0 ) {
                    xor11bit = 1;
                }

                if( strcmp( "downsample", longOpts[longIndex].name ) == 0 ) {
                    downsample = 1;
                }

                if ( strcmp( "noresetxor", longOpts[longIndex].name ) == 0 ) {
                    noresetxor = 1;
                }

                if( strcmp( "xmin", longOpts[longIndex].name ) == 0 ) {
                    gotxmin=1;
                    tempa = atoi(optarg);
                }
                if( strcmp( "xmax", longOpts[longIndex].name ) == 0 ) {
                    gotxmax=1;
                    tempb = atoi(optarg);
                }   
                if( strcmp( "stepnoise", longOpts[longIndex].name ) == 0 ) {
                    modelstate.using_stepnoise = 1;
                    modelstate.stepnoise = atof(optarg);
                }
                if( strcmp( "fast", longOpts[longIndex].name ) == 0 ) {
                    //fprintf(stderr,"FAST OPTION Selected\n");
                    modelstate.fast_m2p=1;
                }
                //if( strcmp( "correlation", longOpts[longIndex].name ) == 0 ) {
                //    modelstate.correlation = atof(optarg);
                //    gotcorrelation=1;
                //}
                if( strcmp( "entropy", longOpts[longIndex].name ) == 0 ) {
                    modelstate.entropy = atof(optarg);
                    gotentropy=1;
                    modelstate.gotentropy=1;
                }
                if( strcmp( "mean", longOpts[longIndex].name ) == 0 ) {
                    modelstate.mean = atof(optarg);
                    //gotmean=1;
                }
                if( strcmp( "variance", longOpts[longIndex].name ) == 0 ) {
                    modelstate.variance = atof(optarg);
                }
                if( strcmp( "states", longOpts[longIndex].name ) == 0 ) {
                    modelstate.states = atoi(optarg);
                }
                if( strcmp( "min_range", longOpts[longIndex].name ) == 0 ) {
                    modelstate.min_range = atof(optarg);
                }
                if( strcmp( "max_range", longOpts[longIndex].name ) == 0 ) {
                    modelstate.max_range = atof(optarg);
                }
                if( strcmp( "sigmoid", longOpts[longIndex].name ) == 0 ) {
                    if (strcmp(optarg,"flat")==0) {
                        modelstate.curve = CURVE_FLAT;
                    }
                    if (strcmp(optarg,"linear")==0) {
                        modelstate.curve = CURVE_LINEAR;
                        strcpy(modelstate.curvestr,"Linear");
                    }
                    if (strcmp(optarg,"logistic")==0) {
                        modelstate.curve = CURVE_LOGISTIC;
                        strcpy(modelstate.curvestr,"Logistic");
                    }
                    
                    if (strcmp(optarg,"tanh")==0) {
                        modelstate.curve = CURVE_TANH;
                        strcpy(modelstate.curvestr,"Hyperbolic Tangent");
                    }
                    
                    if (strcmp(optarg,"atan")==0) {
                        modelstate.curve = CURVE_ATAN;
                        strcpy(modelstate.curvestr,"Arctangent");
                    }
                    
                    if (strcmp(optarg,"gudermann")==0) {
                        modelstate.curve = CURVE_GUDERMANN;
                        strcpy(modelstate.curvestr,"Gudermannian");
                    }
                    
                    if (strcmp(optarg,"erf")==0) {
                        modelstate.curve = CURVE_ERF;
                        strcpy(modelstate.curvestr,"Error Function");
                    }
                    if (strcmp(optarg,"algebraic")==0) {
                        modelstate.curve = CURVE_AGEBRAIC;
                        strcpy(modelstate.curvestr,"Algebraic");
                    }
                }                                  
                if( strcmp( "lcg_a", longOpts[longIndex].name ) == 0 ) {
                    modelstate.lcg_a = strtoull(optarg,NULL,0);
                }
                if( strcmp( "lcg_c", longOpts[longIndex].name ) == 0 ) {
                    modelstate.lcg_c = strtoull(optarg,NULL,0);
                }
                if( strcmp( "lcg_m", longOpts[longIndex].name ) == 0 ) {
                    modelstate.lcg_m = strtoull(optarg,NULL,0);
                }
                if( strcmp( "lcg_truncate", longOpts[longIndex].name ) == 0 ) {
                    modelstate.lcg_truncate = atoi(optarg);
                }
                if( strcmp( "lcg_outbits", longOpts[longIndex].name ) == 0 ) {
                    modelstate.lcg_outbits = atoi(optarg);
                }
                if( strcmp( "pcg_state_size", longOpts[longIndex].name ) == 0 ) {
                    modelstate.pcg_state_size = atoi(optarg);
                }
                if( strcmp( "pcg_generator", longOpts[longIndex].name ) == 0 ) {
                    if (strcmp(optarg, "LCG")==0) {
                        modelstate.pcg_alg = PCG_LCG;
                    }
                    else if (strcmp(optarg, "MCG")==0) {
                        modelstate.pcg_alg = PCG_MCG;
                    }
                }
                if( strcmp( "pcg_of", longOpts[longIndex].name ) == 0 ) {
                    if (strcmp(optarg, "XSH_RR")==0) {
                        modelstate.pcg_of = XSH_RR;
                    }
                    else if (strcmp(optarg, "XSH_RS")==0) {
                        modelstate.pcg_of = XSH_RS;
                    }
                }
                if( strcmp( "xorshift_size", longOpts[longIndex].name ) == 0 ) {
                    modelstate.xorshift_size = atoi(optarg);
                }
                
                if( strcmp( "sinbias_amplitude", longOpts[longIndex].name ) == 0 ) {
                    modelstate.sinbias_amplitude = atof(optarg);
                }
                if( strcmp( "sinbias_offset", longOpts[longIndex].name ) == 0 ) {
                    modelstate.sinbias_offset = atof(optarg);
                }
                if( strcmp( "sinbias_period", longOpts[longIndex].name ) == 0 ) {
                    modelstate.sinbias_period = atoi(optarg);
                }
                
                if( strcmp( "p01", longOpts[longIndex].name ) == 0 ) {
                    modelstate.p01 = atof(optarg);
                    modelstate.gotp01 = 1;
                    gotp01 = 1;
                }
                if( strcmp( "p10", longOpts[longIndex].name ) == 0 ) {
                    modelstate.p10 = atof(optarg);
                    modelstate.gotp10 = 1;
                    gotp10 = 1;
                }
                if( strcmp( "bitwidth", longOpts[longIndex].name ) == 0 ) {
                    modelstate.bitwidth = atoi(optarg);
                    //gotbitwidth = 1;
                }
                                
                break;
            case 'h':   /* fall-through is intentional */
            case '?':
                display_usage();
                exit(0);
                 
            default:
                /* You won't actually get here. */
                break;
        }
         
        opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
    } // end while
    
    /* Sort xmin and xmax */ 
    if ((gotxmin==1) && (gotxmax==1))
    {
        using_xor_range = 1;
        
        if (tempa > tempb)
        {
            xmin = tempb;
            xmax = tempa;
        }
        else
        {
            xmin = tempa;
            xmax = tempb;
        }
    }

    if ((rngstate.randseed==1) && (rngstate.got_detseed)) {
        fprintf(stderr,"Error - Can't have both deterministic and nondeterministic seeding (-s with -D).");
        exit(1);
    }
    
    if (rngstate.randseed==1) {
        if ((rngstate.rdrand_available==0) && (rngstate.devurandom_available==0)){
            fprintf(stderr,"Neither /dev/urandom nor RdRand Supported for nondeterministic seeding.");
            exit(1);
        }
    }
    
    /* start the RNG */
    init_rng(&rngstate);
    
    /* Range check the var args */
    abort = 0;

    
    modelstate.punc_event_dis=0;
    if (modelstate.punc_event_limit==0) modelstate.punc_event_dis=1;

    if ((puncture==1) || (model==MODEL_PUNCTURING)) {
        if (modelstate.puncturing_start < 0) {
            fprintf(stderr,"Error: puncturing_start cannot be negative\n");
            abort = 1;
        }

        if (modelstate.puncturing_length < 1) {
            fprintf(stderr,"Error: puncturing_length must be positive and greater than 0\n");
            abort = 1;
        }

        if (modelstate.puncturing_gap < 0) {
            fprintf(stderr,"Error: puncturing_gap must be positive and greater than 0\n");
            abort = 1;
        }

        if (modelstate.puncturing_level < 0) {
            fprintf(stderr,"Error: puncturing_level must be 0 (low), 1 (high) or 2 (randomized)\n");
            abort = 1;
        }

        if (model == MODEL_NORMAL) {
            fprintf(stderr,"Error: puncturing (--puncture) cannot be used with the normal model (-m MODEL_NORMAL)\n");
            abort = 1;
        }
    }

    if (using_xor_range == 1)
    {
        if (gotxmin != gotxmax)
        {
            fprintf(stderr,"Error: -y and -z (xor_min and xor_max) must be used together or not at all\n");
            abort = 1;
        }
        if (xmin == xmax)
        {
            fprintf(stderr,"Error: -y and -z (xor_min and xor_max) must be different and positive\n");
            abort = 1;
        }
        if (xmin < 1)
        {
            fprintf(stderr,"Error: -y (xor min) must be 1 or greater. Provided value = %d\n",xmin);
            abort = 1;
        }
        if (xmax < 1)
        {
            fprintf(stderr,"Error: -z (xor max) must be 1 or greater. Provided value = %d\n",xmax);
            abort = 1;
        }

        if (xormode==1)
        {
            fprintf(stderr,"Error: -x (xor ratio) cannot be used with the XOR range (-y, -z) options");
            abort = 1;
        }

    }

    if ((xor4bit+xor11bit+downsample) > 1) {
            fprintf(stderr,"Error: Can only use one of --xor4bit, --xor11bit or --downsample at the same time.\n");
            abort=1;
    }

    if ((bits_per_byte != 1) && (bits_per_byte != 2) && (bits_per_byte != 4) && (bits_per_byte != 8))
    {
            fprintf(stderr,"Error: -b -bpb = %d. Bits per byte must be 1, 2, 4 or 8.\n",bits_per_byte);
            abort=1;
    }
    if (rngstate.c_max < 1)
    {
            fprintf(stderr,"Error: -c n: cmax must be an integer of 1 or greater. Supplied value = %d\n",rngstate.c_max);
        
            abort=1;
    }

    if (xorbits < 1)
    {
            fprintf(stderr,"Error: -x n: XOR ratio out of bounds. n must be an integer of 1 or greater. Supplied value = %d\n",xorbits);
            abort=1;
        }

    if (modelstate.sinbias_period < 4)
    {
            fprintf(stderr,"Error: --sinbiad_period=%" PRId64 ": Sinusoid period must be 4 or more\n",modelstate.sinbias_period);
            abort=1;
    }
    
    if (modelstate.left_stepsize == 0.0)
    {
            fprintf(stderr,"Error: -l n: left stepsize cannot be 0.0. Supplied value = %f\n",modelstate.left_stepsize);
            abort=1;
    }
    if (modelstate.right_stepsize == 0.0)
    {
            fprintf(stderr,"Error: -l n: right stepsize cannot be 0.0 . Supplied value = %f\n",modelstate.right_stepsize);
            abort=1;
        }
    if (modelstate.left_stepsize < 0.0)
    {
            fprintf(stderr,"Error: -l n: left stepsize cannot be negative. Supplied value = %f\n",modelstate.left_stepsize);
            abort=1;
        }
    if (modelstate.right_stepsize < 0.0)
    {
            fprintf(stderr,"Error: -l n: right stepsize cannot be negative . Supplied value = %f\n",modelstate.right_stepsize);
            abort=1;
        }
    if (modelstate.stepnoise < 0.0)
    {
            fprintf(stderr,"Error: --stepnoise n: stepnoise cannot be negative . Supplied value = %f\n",modelstate.stepnoise);
            abort=1;
    }
    if ((modelstate.bias < 0.0) || (modelstate.bias >1.0))
    {
        fprintf(stderr,"Error: --bias n: bias must be between 0.0 and 1.0. Supplied value = %f\n",modelstate.bias);
        abort=1;
    }

    if ((modelstate.correlation < -1.0) || (modelstate.correlation >1.0))
    {
        fprintf(stderr,"Error: --correlation n: correlation must be between -1.0 and 1.0. Supplied value = %f\n",modelstate.correlation);
        abort=1;
    }
    
    if ((modelstate.p01 < 0.0) || (modelstate.p01 >1.0))
    {
        fprintf(stderr,"Error: --p01 n: P01 must be between 0.0 and 1.0. Supplied value = %f\n",modelstate.p01);
        abort=1;
    }
    if ((modelstate.p10 < 0.0) || (modelstate.p10 >1.0))
    {
        fprintf(stderr,"Error: --p10 n: P01 must be between 0.0 and 1.0. Supplied value = %f\n",modelstate.p10);
        abort=1;
    }
        
    if (kilobytes < 1)
    {
        fprintf(stderr,"Error: -k n: Output size must be 1 or more kilobytes. n must be an integer of 1 or greater. Supplied value = %d\n",kilobytes);
        abort=1;
    }

    if ((width >256) || (width < 1))
    {
        fprintf(stderr,"Error: Width must be from 1 to 256\n");
        abort = 1;  
    } 
    if (model==MODEL_PCG) {
        if ((modelstate.pcg_state_size != 16) && (modelstate.pcg_state_size != 32)
              && (modelstate.pcg_state_size != 64) && (modelstate.pcg_state_size != 128)) {
            fprintf(stderr,"Error: A pcg_size must be one of 16, 32, 64 or 128\n");
            exit(1);
        }    
    }
    if (model==MODEL_XORSHIFT) {
        if ((modelstate.xorshift_size != 32) &&  (modelstate.xorshift_size != 128)) {
            fprintf(stderr,"Error: A xorshift_size must be one of 32 or 128\n");
            exit(1);
        }    
    }
    if (model==MODEL_MARKOV2P) {
        if (((gotcorrelation==1) || (gotbias==1)) &&  ((gotp01) || (gotp10))) {
            fprintf(stderr,"Error: Cannot give both correlation,bias and p01,p10 parameters with Markov model\n");
            exit(1);
        }    
        if (((gotcorrelation==1) || (gotbias==1)) &&  (gotentropy)) {
            fprintf(stderr,"Error: Cannot give both correlation,bias and entropy parameters with Markov model\n");
            exit(1);
        }  
        if (((gotp01==1) || (gotp10==1)) &&  (gotentropy)) {
            fprintf(stderr,"Error: Cannot give both p01,p10 and entropy parameters with Markov model\n");
            exit(1);
        }  
        
        // Deal with the 3 parameter types
        /*
        if (gotentropy==1) {
            epsilon = pow(2.0,-50);
            pick_point(&(modelstate.p01),&(modelstate.p10),modelstate.entropy,epsilon,modelstate.bitwidth,&rngstate);
            modelstate.bias = modelstate.p01/(modelstate.p10+modelstate.p01);
            modelstate.correlation = 1.0 - modelstate.p01 - modelstate.p10;   
        }
        else if ((gotbias == 1) || (gotcorrelation==1)) {
            if (gotbias==0){
                 modelstate.bias = 0.5;
            }
            
            if (gotcorrelation==0) {
                modelstate.correlation = 0.0;
            }
            modelstate.p01 = modelstate.bias * (1.0 - modelstate.correlation);
            modelstate.p10 = (1.0-modelstate.bias)*(1.0-modelstate.correlation);
            //fprintf(stderr,"bias = %f\n",modelstate.bias);
            //fprintf(stderr,"correlation = %f\n",modelstate.correlation);
            //fprintf(stderr,"p01 = %f\n",modelstate.p01);
            //fprintf(stderr,"p10 = %f\n",modelstate.p10);
        } else if ((gotp01 == 1) || (gotp10==1))  {
            if (gotp01==0) modelstate.p01 = 0.5;
            if (gotp10==0) modelstate.p10 = 0.5;
            modelstate.correlation = 1.0 - modelstate.p10 - modelstate.p01;
            modelstate.bias  = modelstate.p10/(modelstate.p10+modelstate.p01);
        } else if (gotentropy==0) {
            modelstate.entropy = 1.0;
            modelstate.p01 = 0.5;
            modelstate.p10 = 0.5;
            modelstate.bias = 0.5;
            modelstate.correlation = 0.0;
        }
        */
        
        
    }
    if (model==MODEL_FILE) {
        if (modelstate.using_infile == 0) {
            fprintf(stderr,"Error: A file must be provided for the file input model using -i <filename> or --infile=<filename>\n");
        abort = 1;
        } 
    }
    if (abort==1) {
        exit(1);
    }


    
    /* Print out the job parameters */
    if (verbose_mode>0)
    {
        if (aesni_check_support() == 1)
            fprintf(stderr,"AESNI Supported in instruction set\n");
        else
            fprintf(stderr,"AESNI Not supported in instruction set\n");
 
        if (hex_mode == 1)
            fprintf(stderr,"Format=Hex\n");
        else if (binary_mode == 1)
            fprintf(stderr,"Format=Binary\n");
        else if (nistoddball_mode == 1)
            fprintf(stderr,"Format=Nist Oddball");

        if (model == MODEL_SUMS)
        {
            fprintf(stderr,"model=sums\n");
            fprintf(stderr,"  left stepsize  = %f\n",modelstate.left_stepsize);
            fprintf(stderr,"  right stepsize = %f\n",modelstate.right_stepsize);
            if (modelstate.using_stepnoise==0)
            {
                fprintf(stderr,"  Not adding noise to stepsize\n");
            }
            else
                fprintf(stderr,"  Adding noise of variance %f to stepsize\n",modelstate.stepnoise);
        }
        if ((model == MODEL_PUNCTURING) || (puncture == 1))
        {
            if (puncture == 1) fprintf(stderr,"Puncturing enabled\n");
            else fprintf(stderr,"model=puncturing\n");
            if (modelstate.puncturing_level == PUNCTURING_LEVEL_LOW) {
                fprintf(stderr,"  level  = zeroes\n");
            } else if (modelstate.puncturing_level == PUNCTURING_LEVEL_HIGH) {
                fprintf(stderr,"  level  = ones\n");
            } else if (modelstate.puncturing_level == PUNCTURING_LEVEL_BOTH) {
                fprintf(stderr,"  level  = zeroes or ones randomly\n");
            }
            
            fprintf(stderr,"  start bit position         = %d\n", modelstate.puncturing_start);
            fprintf(stderr,"  length of injected string  = %d\n", modelstate.puncturing_length);
            fprintf(stderr,"  gap between injections     = %d\n", modelstate.puncturing_gap);
        }
        if (model == MODEL_PURE)
        {
            fprintf(stderr,"model=pure\n");
        }

        if (model == MODEL_BIASED)
        {
            fprintf(stderr,"model=biased\n");
            fprintf(stderr,"  bias  = %f\n",modelstate.bias);
        }

        if (model == MODEL_CORRELATED)
        {
            fprintf(stderr,"model=correlated\n");
            fprintf(stderr,"  correlation  = %f\n",modelstate.correlation);
        }

        //if (model == MODEL_MARKOV2P)
        //{
        //    double entropy;
        //    double lmcv_prob;
        //    uint64_t lmcv;
        //
        //    entropy = p_to_entropy(modelstate.p01, modelstate.p10,modelstate.bitwidth, &lmcv_prob, &lmcv) ;
        //    
        //    fprintf(stderr,"model=markov_2_param\n");
        //    fprintf(stderr,"  bias            = %f\n",modelstate.bias);
        //    fprintf(stderr,"  correlation     = %f\n",modelstate.correlation);
        //    fprintf(stderr,"  p01             = %f\n",modelstate.p01);
        //    fprintf(stderr,"  p10             = %f\n",modelstate.p10);
        //    fprintf(stderr,"  entropy         = %f\n",entropy);
        //    fprintf(stderr,"  MCV Prob        = %f\n",lmcv_prob);
        //    
        //}
        if (model == MODEL_MARKOV_SIGMOID)
        {
            fprintf(stderr,"model=markov_sigmoid\n");
            fprintf(stderr,"  Chain Length    = %d\n",modelstate.states);
            fprintf(stderr,"  curve           = %s\n",modelstate.curvestr);
            fprintf(stderr,"  range from      = %f\n",modelstate.min_range);
            fprintf(stderr,"          to      = %f\n",modelstate.max_range);
        }
                
        if (model == MODEL_LCG)
        {
            fprintf(stderr,"model=linear congruential generator\n");
            fprintf(stderr,"  a  = 0x%llx\n",modelstate.lcg_a);
            fprintf(stderr,"  c  = 0x%llx\n",modelstate.lcg_c);
            fprintf(stderr,"  m  = 0x%llx\n",modelstate.lcg_m);
            fprintf(stderr,"  start x = 0x%llx\n",modelstate.lcg_x % modelstate.lcg_m);
            fprintf(stderr,"  Output bit field = %d:%d\n",
                    (modelstate.lcg_truncate)+(modelstate.lcg_outbits)-1,modelstate.lcg_truncate);
        }
        if (model == MODEL_PCG)
        {
            fprintf(stderr,"model=permuted congruential generator\n");
            fprintf(stderr,"  state size = %d\n",modelstate.pcg_state_size);
            if (modelstate.pcg_alg == PCG_MCG)
                fprintf(stderr,"  State update algorithm MCG\n");
            else if (modelstate.pcg_alg == PCG_LCG)
                fprintf(stderr,"  State update algorithm LCG\n");
            else
                fprintf(stderr,"  Unknown State Update Function %d\n", modelstate.pcg_alg);
                    
            if (modelstate.pcg_of == XSH_RS)
                fprintf(stderr,"  Output function XSH_RS\n");
            else if (modelstate.pcg_of == XSH_RR)
                fprintf(stderr,"  Output function XSH_RR\n");
            else
                fprintf(stderr,"  Unknown Output function %d\n", modelstate.pcg_of);
        }
        if (model == MODEL_XORSHIFT)    
        {
            fprintf(stderr,"model=XORSHIFT\n");
        }   
        if (model == MODEL_NORMAL)
        {
            fprintf(stderr,"model=normal\n");
            fprintf(stderr,"  mean  = %f\n",modelstate.mean);
            fprintf(stderr,"  variance  = %f\n",modelstate.variance);
        }
        if (model == MODEL_NORMALINTEGER)
        {
            fprintf(stderr,"model=normalinteger\n");
            fprintf(stderr,"  mean  = %f\n",modelstate.mean);
            fprintf(stderr,"  variance  = %f\n",modelstate.variance);
        }

        if (model == MODEL_FILE)
        {
            fprintf(stderr,"model=file\n");
            fprintf(stderr,"  filename  = %s\n",infilename);
            if (rngstate.input_format == INFORMAT_HEX)
                fprintf(stderr,"  File input format Hex\n");
            else if (rngstate.input_format == INFORMAT_BINARY)
                fprintf(stderr,"  File input format Binary\n");
            else if (rngstate.input_format == INFORMAT_01 )
                fprintf(stderr,"  File input format ASCII Binary\n");
            else
                fprintf(stderr,"  Unknown input format :%d\n",rngstate.input_format);
        }

        fprintf(stderr,"size = %d kilobytes\n", kilobytes);

        if (xor4bit==1)
            fprintf(stderr,"4 bit decorrelator/decimator being used");
        
        if (xor11bit==1)
            fprintf(stderr,"11 bit decorrelator/decimator being used");
        
        if (downsample==1)
            fprintf(stderr,"16-4 downsampler being used");

        if (xormode == 1)
            fprintf(stderr,"XOR mode on, fixed ratio=%d:1\n",xorbits);
        else
            fprintf(stderr,"XOR mode off\n");

        if (using_xor_range == 1)
            fprintf(stderr,"XOR range mode on, ratios between %d:1 and %d:1, chosen randomly\n",xmin, xmax);
        else
            fprintf(stderr,"XOR range mode off\n");

        if (ofile == 0)
            fprintf(stderr,"Output to STDOUT\n");
        else
            fprintf(stderr,"Output to file %s\n",filename);

        if (rngstate.randseed==1)
        {
            fprintf(stderr,"Hardware Random Seeding on. Non deterministic mode\n");
            fprintf(stderr,"  Reseed c_max=%d\n",rngstate.c_max);
            if (rngstate.rdrand_available) {
                fprintf(stderr,"  Using RdRand as nondeterministic source\n");
            }
            else if(rngstate.devurandom_available) {
                fprintf(stderr,"  Using /dev/urandom as nondeterministic source");
            }   
        }
        else
        {
            fprintf(stderr,"Hardware Random Seeding off. Deterministic mode.\n");
            fprintf(stderr,"  Restir c_max=%d\n",rngstate.c_max);
        }
        
        if (modelstate.using_jfile==1)
        {
            fprintf(stderr,"Outputting internal per bit bias to file %s\n",jfilename);
        }
    }

    /* open the output file if needed */

    if (ofile==1)
    {
        fp = fopen(filename, "wb");
        if (fp == NULL) {
            sprintf(errstr,"failed to open output file %s for writing",filename);
            perror(errstr);
            exit(1);
        }
    }

    /* open the j file if needed */
    if (modelstate.using_jfile==1)
    {
        modelstate.jfile = fopen(jfilename, "wb");
        if (modelstate.jfile == NULL) {
            sprintf(errstr,"failed to open output j file %s for writing",jfilename);
            perror(errstr);
            exit(1);
        }
    }

    /* open the input file if needed */
    if (modelstate.using_infile==1)
    {
        if (rngstate.input_format==INFORMAT_BINARY)
            modelstate.infile = fopen(infilename,"rb");
        else
            modelstate.infile =  fopen(infilename, "r");
        if (modelstate.infile == NULL) {
            sprintf(errstr,"Failed to input file %s for reading",infilename);
            perror(errstr);
            exit(1);
        }
        else
        {
            if (verbose_mode>0)
                fprintf(stderr,"opened input file %s for reading\n",infilename);
        }
    }

    /* Open the JSON File if needed */
    if (modelstate.using_json==1) {
        modelstate.json_file = fopen(json_filename,"w");
        if (modelstate.json_file == NULL) {
            sprintf(errstr,"Failed to open JSON output file %s for writing\n",json_filename);
            perror(errstr);
            exit(1);
        }
        else {
            if (verbose_mode>0)
                fprintf(stderr,"opened JSON ouput file %s for writing\n",json_filename);
        }
    }

    /* Open the YAML File if needed */
    if (modelstate.using_yaml==1) {
        modelstate.yaml_file = fopen(yaml_filename,"w");
        if (modelstate.yaml_file == NULL) {
            sprintf(errstr,"Failed to open YAML output file %s for writing\n",yaml_filename);
            perror(errstr);
            exit(1);
        }
        else {
            if (verbose_mode>0)
                fprintf(stderr,"opened YAML ouput file %s for writing\n",yaml_filename);
        }
    }

    /* Initialize the RNG */
    initialize_sim(model, &modelstate, &rngstate);

    /* For each stepsize, perform the simulation over 256 samples */
    /* And do it 4 times */

    /* entropy = 0x00;*/
    /* Pull some bits to let it settle. But not in the following modes since there is not settling time needed */
    if (!((model==MODEL_PURE) || (model==MODEL_PUNCTURING) || (model==MODEL_FILE) || (model==MODEL_NORMAL) || (model==MODEL_NORMALINTEGER) || (model==MODEL_LCG) || (model==MODEL_PCG) || (model==MODEL_MARKOV2P) || (model==MODEL_MARKOV_SIGMOID)))
    for(i=0; i<128; i++)
    {
        thebit = entropysource(model, &modelstate, &rngstate);
        modelstate.lastbit = thebit;
    }

    /* Start with pulling samples and testing them */
    int32_t normintval; //2s comp value from normal integer model

    if (model==MODEL_NORMALINTEGER) /* or any other floating point model that is added */
    {
        
        for (simrun =0; simrun < kilobytes; simrun++)
        {
            for (onek=0;onek<32;onek++)
            {
                normintval = entropysource(model, &modelstate, &rngstate);
                normalintegersamples[onek]=normintval;
            }
            /* Output the 32 value block */

            if (ofile == 1 && binary_mode==1) /* binary to a file */
            {
                fwrite(normalintegersamples, 128, 1, fp);
            }
            else if (ofile == 0 && binary_mode == 0) /* ascii text to stdout */
            {
                for (j=0;j<32;j++)
                {
                    fprintf(stdout,"%08X\n",normalintegersamples[j]);
                }
            }
            else if (ofile == 1 && binary_mode == 0) /* Floatingpoint text to a file */
            {
                for (j=0;j<32;j++)
                {
                    fprintf(fp,"%08X\n",normalintegersamples[j]);
                }

            }
            else /* binary to stdout */
            {
                fwrite(normalintegersamples, 32, 1, stdout);
            }
        }
    } else if (model==MODEL_NORMAL) {
        for (simrun =0; simrun < kilobytes; simrun++)
        {
            for (onek=0;onek<4;onek++)
            {
                for(samplenum=0;samplenum<256;samplenum++)
                {
                    thevalue = entropysource(model, &modelstate, &rngstate);
                    floatingpointsamples[samplenum]=thevalue;
                }

                /* Output the 256 value block */

                if (ofile == 1 && binary_mode==1) /* binary to a file */
                {
                    fwrite(floatingpointsamples, 256*sizeof(double), 1, fp);
                }
                else if (ofile == 0 && binary_mode == 0) /* floatingpoint text to stdout */
                {
                    for (j=0;j<256;j++)
                    {
                        fprintf(stdout,"%0.8f\n",floatingpointsamples[j]);
                    }
                }
                else if (ofile == 1 && binary_mode == 0) /* Floatingpoint text to a file */
                {
                    for (j=0;j<256;j++)
                    {
                        fprintf(fp,"%0.8f\n",floatingpointsamples[j]);
                    }

                }
                else /* binary to stdout */
                {
                    fwrite(floatingpointsamples, 256*sizeof(double), 1, stdout);
                }
            }

        }

    } else /* model = sums, pure, puncturing, biased, correlated, lcg, pcg, xorshift, markov_2_param, markov_sigmoid or file*/
    {
        lineindex = 0;
        
        for (simrun =0; ((simrun < kilobytes) || ((model==MODEL_FILE) && (no_k==1))); simrun++)
        {
            for (onek=0;onek<4;onek++)
            {
                for(samplenum=0;samplenum<256;samplenum++) {
                    thebyte = (unsigned char)0x00;
                    
                    // Reset the shift register when OSTE buffers fill.
                    if (noresetxor == 0) {
                        if ((xorcount*8) >= 512) {
                            xorcount = 0;
                            shiftreg = 0;
                        }
                        xorcount++;
                    }

                    /* Pull 8 bits */
                    if (xor4bit==1) {
                        // Every 512 bits, clear the shiftregister.
                        // This models the point where the OSTE registers have been
                        // consumed by the conditioner and the conditioner is busy
                        // So the ES is paused and reset.

                        
                        // The first 4 bits of the byte
                        for (xoriter=0; xoriter < 16; xoriter++) {
                            newbit = entropysource(model, &modelstate, &rngstate);
                            modelstate.lastbit = newbit;
                            //shiftreg = (((shiftreg & 0x1)^newbit)<<3) | ((shiftreg>>1) & 0x7);
                            shiftreg = ((((shiftreg & 0x8)>>3)^newbit) & 0x1) | ((shiftreg<<1) & 0xe);
                        }
                        thebyte = (shiftreg & 0xf);
                        
                        // The second 4 bits of the byte
                        for (xoriter=0; xoriter < 16; xoriter++) {
                            newbit = entropysource(model, &modelstate, &rngstate);
                            modelstate.lastbit = newbit;
                            shiftreg = ((((shiftreg & 0x8)>>3)^newbit) & 0x1) | ((shiftreg<<1) & 0xe);
                            //shiftreg = (((shiftreg & 0x1)^newbit)<<3) | ((shiftreg>>1) & 0x7); // new bit and bit 0 xored into bit 3. bits 3-1 shifted to 2-0.
                        }
                        //thebyte = thebyte | ((shiftreg & 0xf) << 4);
                        thebyte = (thebyte << 4) | (shiftreg & 0xf);
                    } else if (xor11bit==1) {
                        // The first 4 bits of the byte
                        for (xoriter=0; xoriter < 16; xoriter++) {
                            newbit = entropysource(model, &modelstate, &rngstate);
                            modelstate.lastbit = newbit;
                            shiftreg = ((((shiftreg & 0x400)>>10)^newbit) & 0x1) | ((shiftreg<<1) & 0x7fe);
                            //shiftreg = (((shiftreg & 0x1)^newbit)<<10) | ((shiftreg>>1) & 0x3ff); // new bit and bit 0 xored into bit 10. bits 10-1 shifted to 9-0. 
                        }
                        thebyte = (shiftreg & 0xf);
                        
                        // The second 4 bits of the byte
                        for (xoriter=0; xoriter < 16; xoriter++) {
                            newbit = entropysource(model, &modelstate, &rngstate);
                            modelstate.lastbit = newbit;
                            shiftreg = ((((shiftreg & 0x400)>>10)^newbit) & 0x1) | ((shiftreg<<1) & 0x7fe);
                            //shiftreg = (((shiftreg & 0x1)^newbit)<<10) | ((shiftreg>>1) & 0x3ff); // new bit and bit 0 xored into bit 10. bits 10-1 shifted to 9-0. 
                        }
                        //thebyte = thebyte | ((shiftreg & 0xf) << 4);
                        thebyte = (thebyte << 4) | (shiftreg & 0xf);
                    } else if (downsample==1) {
                        // The first 4 bits of the byte
                        //inputbits = 0;
                        for (xoriter=0; xoriter < 16; xoriter++) {
                            newbit = entropysource(model, &modelstate, &rngstate);
                            modelstate.lastbit = newbit;
                            //inputbits = (inputbits << 1) | (newbit & 0x1);
                            shiftreg = ((newbit & 0x1) | ((shiftreg<<1) & 0x0e)); // new bit to bit 3. bits 3-1 shifted to 2-0.
                            
                            //shiftreg = ((newbit<<3) | ((shiftreg>>1) & 0x7)); // new bit to bit 3. bits 3-1 shifted to 2-0. 
                        }
                        thebyte = (shiftreg & 0xf);
                        //fprintf(stderr, "\ndownsample got %04x in.",inputbits);
                        //fprintf(stderr, " Returned %x.\n",thebyte);
                        
                        // The second 4 bits of the byte
                        //inputbits = 0;
                        for (xoriter=0; xoriter < 16; xoriter++) {
                            newbit = entropysource(model, &modelstate, &rngstate);
                            modelstate.lastbit = newbit;
                            //inputbits = (inputbits << 1) | (newbit & 0x1);
                            shiftreg = ((newbit & 0x1) | ((shiftreg<<1) & 0x0e)); // new bit to bit 3. bits 3-1 shifted to 2-0. 
                            //shiftreg = ((newbit<<3) | ((shiftreg>>1) & 0x7)); // new bit to bit 3. bits 3-1 shifted to 2-0. 
                        }
                        //fprintf(stderr, "\ndownsample got %04x in",inputbits);
                        //thebyte = thebyte | ((shiftreg & 0xf) << 4);
                        thebyte = (thebyte << 4) | (shiftreg & 0xf);
                        //fprintf(stderr, " Returned %02x.\n",thebyte);
                    } else if ((model==MODEL_MARKOV2P) && (modelstate.fast_m2p==1)) {
                        thebyte = markov2pfastsource(&modelstate, &rngstate);
                        //if (rngstate.reached_eof == 1) {
                        //    if ((samplenum > 0) && (samplenum < 256)) {
                        //        /*fprintf(stderr,"reached EOF with samplenum > 0 = %d\n",samplenum);*/
                        //        goto eof_with_partial_block;
                        //    }
                        //    else {
                        //        /*fprintf(stderr,"Going to EOF with full block. samplenum = %d\n",samplenum);*/
                        //        goto reached_eof;
                        //    }
                        //}
                        //modelstate.lastbit = thebit;
                        //prob = modelstate.bias;
                    } else {

                        for(i=0; i<8; i++)
                        {
                            if (xormode==1)
                            {
                                for (xoriter=0; xoriter < xorbits; xoriter++)
                                {
                                    thebit ^= entropysource(model, &modelstate, &rngstate);
                                    if (rngstate.reached_eof == 1)
                                        goto reached_eof;
                                    modelstate.lastbit = thebit;
                                    prob = modelstate.bias;
                                    if (xoriter == 0) pa1 = prob;
                                    if (xoriter == 1) pb1 = prob;
                                    if (xoriter == 2) pc1 = prob;
                                }
                                if (xorbits==2) postxor_entropy = bias2entropy(xorbias_2bit(pa1,pb1));
                                if (xorbits==3) postxor_entropy = bias2entropy(xorbias_3bit(pa1,pb1,pc1));
                                total_entropy += postxor_entropy;
                            }
                            else if (using_xor_range==1)
                            {
                                xorrange = 1+xmax-xmin;
                                xorselector = getrand16(&rngstate);
                                xorselector = xorselector % xorrange;
                                xorbits = xorselector;
        
                                for (xoriter=0; xoriter < (xmin+xorbits); xoriter++)
                                {
                                    thebit ^= entropysource(model, &modelstate, &rngstate);
                                    if (rngstate.reached_eof == 1) goto reached_eof;
                                    modelstate.lastbit = thebit;
                                    //prob = modelstate.bias;
                                    //sums_entropy = bias2entropy(prob);
                                }
                            }
                            else /* no xorring */
                            {
                                thebit = entropysource(model, &modelstate, &rngstate);
                                if (rngstate.reached_eof == 1) {
                                    if ((samplenum > 0) && (samplenum < 256)) {
                                        /*fprintf(stderr,"reached EOF with samplenum > 0 = %d\n",samplenum);*/
                                        goto eof_with_partial_block;
                                    }
                                    else {
                                        /*fprintf(stderr,"Going to EOF with full block. samplenum = %d\n",samplenum);*/
                                        goto reached_eof;
                                    }
                                }
                                modelstate.lastbit = thebit;
                                prob = modelstate.bias;
                                //sums_entropy = bias2entropy(prob);
                                //total_entropy += sums_entropy;
                            }

                            if ((thebit & 0x01)==1)
                            {
                                //fprintf(stderr,"RETURNED BIT = 1\n");
                                //thebyte = ((thebyte >> 1) & 0x7f) | 0x80;
                                thebyte = (thebyte << 1) | 0x01;
                            }
                            else
                            {
                                //fprintf(stderr,"RETURNED BIT = 0\n");
                                //thebyte = ((thebyte >> 1) & 0x7f);
                                thebyte = (thebyte << 1) & 0xfe;
                            }
                        } // end for i = 1..7
                    } // end else not fast

                    thesample[samplenum]=thebyte;
                    //fprintf(stderr, "non downsample returned %02X\n",thebyte);
                }

                /* Output the 256 byte block */
                eof_with_partial_block:
                if (ofile == 1 && binary_mode==1) /* binary to a file */
                {
                    if (bits_per_byte == 8) {
                        fwrite(thesample, samplenum, 1, fp);
                    } else if (bits_per_byte == 4) {
                        for (i=0;i<samplenum;i++) {
                            abyte = thesample[i];
                            thebpbsample[i*2] = abyte & 0x0f;
                            thebpbsample[(i*2)+1] = (abyte >> 4) & 0x0f;
                        }
                        fwrite(thebpbsample, (samplenum * 2) ,1, fp);
                    } else if (bits_per_byte == 2) {
                        for (i=0;i<samplenum;i++) {
                            abyte = thesample[i];
                            thebpbsample[i*4] = abyte & 0x03;
                            thebpbsample[(i*4)+1] = (abyte >> 2) & 0x03;
                            thebpbsample[(i*4)+2] = (abyte >> 4) & 0x03;
                            thebpbsample[(i*4)+3] = (abyte >> 6) & 0x03;
                        }
                        fwrite(thebpbsample, (samplenum * 4) ,1, fp);
                    } else if (bits_per_byte == 1) {
                        for (i=0;i<samplenum;i++) {
                            abyte = thesample[i];
                            thebpbsample[i*8] = abyte & 0x01;
                            thebpbsample[(i*8)+1] = (abyte >> 1) & 0x01;
                            thebpbsample[(i*8)+2] = (abyte >> 2) & 0x01;
                            thebpbsample[(i*8)+3] = (abyte >> 3) & 0x01;
                            thebpbsample[(i*8)+4] = (abyte >> 4) & 0x01;
                            thebpbsample[(i*8)+5] = (abyte >> 5) & 0x01;
                            thebpbsample[(i*8)+6] = (abyte >> 6) & 0x01;
                            thebpbsample[(i*8)+7] = (abyte >> 7) & 0x01;
                        }
                        fwrite(thebpbsample, (samplenum * 8) ,1, fp);
                    }
                }
                else if (ofile == 0 && binary_mode == 0) /* hex to stdout */
                {
                    tempindex = 0;
                    /* fprintf(stderr," Samplenum %d\n",samplenum);*/
                    do {
                        printf("%02X",thesample[tempindex++]);
                        lineindex++;
                        if (lineindex == linewidth) {
                            printf("\n");
                            lineindex = 0;
                        }
                    } while (tempindex <samplenum);
                    if (lineindex != 0) printf("\n");
                }
                else if (ofile == 1 && binary_mode == 0) /* Hex to a file */
                {
                    tempindex = 0;
                    lineindex = 0;
                    
                    do {
                        fprintf(fp,"%02X",thesample[tempindex++]);
                        lineindex++;
                        if (lineindex == linewidth) {
                            fprintf(fp,"\n");
                            lineindex = 0;
                        }
                    } while (tempindex < samplenum);
                    if (lineindex != 0) fprintf(fp,"\n");
                }
                else /* binary to stdout */
                {
                    if (bits_per_byte == 8) {
                        fwrite(thesample, samplenum, 1, stdout);
                    } else if (bits_per_byte == 4) {
                        for (i=0;i<samplenum;i++) {
                            abyte = thesample[i];
                            thebpbsample[i*2] = abyte & 0x0f;
                            thebpbsample[(i*2)+1] = (abyte >> 4) & 0x0f;
                        }
                        fwrite(thebpbsample, (samplenum * 2) ,1, stdout);
                    } else if (bits_per_byte == 2) {
                        for (i=0;i<samplenum;i++) {
                            abyte = thesample[i];
                            thebpbsample[i*4] = abyte & 0x03;
                            thebpbsample[(i*4)+1] = (abyte >> 2) & 0x03;
                            thebpbsample[(i*4)+2] = (abyte >> 4) & 0x03;
                            thebpbsample[(i*4)+3] = (abyte >> 6) & 0x03;
                        }
                        fwrite(thebpbsample, (samplenum * 4) ,1, stdout);
                    } else if (bits_per_byte == 1) {
                        for (i=0;i<samplenum;i++) {
                            abyte = thesample[i];
                            thebpbsample[i*8] = abyte & 0x01;
                            thebpbsample[(i*8)+1] = (abyte >> 1) & 0x01;
                            thebpbsample[(i*8)+2] = (abyte >> 2) & 0x01;
                            thebpbsample[(i*8)+3] = (abyte >> 3) & 0x01;
                            thebpbsample[(i*8)+4] = (abyte >> 4) & 0x01;
                            thebpbsample[(i*8)+5] = (abyte >> 5) & 0x01;
                            thebpbsample[(i*8)+6] = (abyte >> 6) & 0x01;
                            thebpbsample[(i*8)+7] = (abyte >> 7) & 0x01;
                        }
                        fwrite(thebpbsample, (samplenum * 8) ,1, stdout);
                    }
                    /*fwrite(thesample, samplenum, 1, stdout);*/
                }
            }
        }
        reached_eof:
        
        if (verbose_mode>0)
        {
            fprintf(stderr,"Total Entropy = %F\n",total_entropy);
            fprintf(stderr,"Per bit Entropy = %F %% \n",(100.0*(total_entropy/(8.0*kilobytes*1024.0))));
        }
    }

    if (fp                   != NULL) fclose(fp);
    if (modelstate.jfile     != NULL) fclose(modelstate.jfile);
    if (modelstate.infile    != NULL) fclose(modelstate.infile);
    if (modelstate.json_file != NULL) fclose(modelstate.json_file);
    if (modelstate.yaml_file != NULL) fclose(modelstate.yaml_file);
    return 0;

}



