
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
/* make isnan() visible */
/*#define _BSD_SOURCE */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>
#include <stdint.h>

#include "aes128k128d.h"
#include "djenrandommodel.h"
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include "rdrand.h"

#define MODEL_SUMS 0
#define MODEL_PURE 1
#define MODEL_BIASED 2
#define MODEL_CORRELATED 3
#define MODEL_NORMAL 4 
#define MODEL_FILE 5
#define MODEL_LCG 6
#define MODEL_PCG 7
#define MODEL_XORSHIFT 8

#define INFORMAT_01 0
#define INFORMAT_HEX 1
#define INFORMAT_BINARY 2

int aesni_supported;

void display_usage() {
fprintf(stderr,"Usage: djrandom [-bsvh] [-x <bits>] [-y <bits>] [-z <bits>] [-c <generate length>]\n");
fprintf(stderr,"       [-m <|pure(default)|sums|biased|correlated|normal|file>] [-l <left_stepsize>]\n"); 
fprintf(stderr,"       [-r <right_stepsize>] [--stepnoise=<noise on step>] [--bias=<bias>]\n");
fprintf(stderr,"       [--correlation=<correlation>] [--mean=<normal mean>] [--variance=<normal variance>]\n");
fprintf(stderr,"       [--pcg_state_16=<16|32|64>] [--pcg_generator=<LCG|MCG>] [--pcg_of=<XSH_RS|XSH|RR]\n");
fprintf(stderr,"       [-o <output_filename>] [-j <j filename>] [-i <input filename>] [-f <hex|binary|01>]\n");
fprintf(stderr,"       [-k <1K_Blocks>] [-w [1..256]]\n");
fprintf(stderr,"\n");
fprintf(stderr,"Generate random bits with configurable non-uniformities.\n");
fprintf(stderr,"  Author: David Johnston, dj@deadhat.com\n");
fprintf(stderr,"\n");

fprintf(stderr,"  -m, --model=<pure(default)|sums|biased|correlated|lcg|pcg|xorshift|normal|file>   Select random source model\n");
fprintf(stderr,"\nStep Update Metastable Source model (-m sums) Options\n\n");
fprintf(stderr,"  -l, --left=<left_stepsize>     stepsize when moving left as a fraction of sigma_m.\n");
fprintf(stderr,"  -r, --right=<right_stepsize>   stepsize when moving right as a fraction of sigma_m.\n");
fprintf(stderr,"  --stepnoise=<noise on step>    variance of the noise on stepsize. e.g. 0.00001.\n");
fprintf(stderr,"\nBiased model (-m biased) Options\n\n");
fprintf(stderr,"  --bias=<bias>                  bias as a number between 0.0 and 1.0. Only for biased model\n");
fprintf(stderr,"\nCorrelated model (-m correlated) Options\n\n");
fprintf(stderr,"  --correlation=<correlation>    correlation with previous bit as a number between -1.0 and 1.0. Only for correlation model\n");
fprintf(stderr,"\nNormal model (-m normal) Options\n\n");
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
fprintf(stderr,"\nGeneral Options\n\n");
fprintf(stderr,"  -x, --xor=<bits>               XOR 'bits' of entropy together for each output bit\n");
fprintf(stderr,"  -y, --xmin=<bits>              Provides the start of a range of XOR ratios to be chosen at random per sample\n");
fprintf(stderr,"  -z, --xmax=<bits>              Provides the end of a range of XOR ratios to be chosen at random per sample\n");
fprintf(stderr,"  -s, --seed                     seed the internal RNG with /dev/random\n");
fprintf(stderr,"  -c, --cmax=<generate length>   number of PRNG generates before a reseed\n");
fprintf(stderr,"  -v, --verbose                  output the parameters\n");
fprintf(stderr,"\nFile Options\n\n");
fprintf(stderr,"  -o <output_filename>             output file\n");
fprintf(stderr,"  -j, --jfile=<j filename>         filename to push source model internal state to\n");
fprintf(stderr,"  -i, --infile=<input filename>    filename of entropy file for file model\n");
fprintf(stderr,"  -f, --informat=<hex|binary|01>   Format of input file. hex=Ascii hex(default), 4 bit per hex character. binary=raw binary. 01=ascii binary. Non valid characters are ignored\n");
fprintf(stderr,"  -k, --blocks=<1K_Blocks>         Size of output in kilobytes\n");
fprintf(stderr,"\nOutput Format Options\n\n");
fprintf(stderr,"  -b, --binary                output in raw binary format\n");
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
        result = log(x)/log(2);
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
	
		if (isnan(result)) printf(" Got NAN in bias2entropy bias=%f\n",bias);
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

	if (isnan(result)) printf(" Got NAN in xorbias_2bit\n");
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

	if (isnan(result)) printf(" Got NAN in xorbias_3bit\n");
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
	if (model==MODEL_SUMS)
	{
		result = smoothsource(modelstate, rngstate);
		return result;
	}
	else if (model==MODEL_PURE)
 	{
		result = puresource(modelstate, rngstate);
		return result;
	}
	else if (model==MODEL_BIASED)
 	{
		result = biasedsource(modelstate, rngstate);
		return result;
	}
	else if (model==MODEL_CORRELATED)
 	{
		result = correlatedsource(modelstate, rngstate);
		return result;
	}
	else if (model==MODEL_LCG)
 	{
		result = lcgsource(modelstate, rngstate);
		return result;
	}
	else if (model==MODEL_PCG)
 	{
		result = pcgsource(modelstate, rngstate);
		return result;
	}
	else if (model==MODEL_XORSHIFT)
 	{
		result = xorshiftsource(modelstate, rngstate);
		return result;
	}
	else if (model==MODEL_FILE)
 	{
		if (rngstate->input_format==INFORMAT_HEX)
		    result = filesourcehex(modelstate,rngstate);
		else if (rngstate->input_format==INFORMAT_01)
		    result = filesource(modelstate,rngstate);
		else
			result = filesourcebinary(modelstate, rngstate);
		
		return result;
	}
	else return 0;
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
	else if (model==MODEL_BIASED)
 	{
		biasedinit(modelstate, rngstate);
	}
	else if (model==MODEL_CORRELATED)
 	{
		correlatedinit(modelstate, rngstate);
	}
	else if (model==MODEL_NORMAL)
	{
		normalinit(modelstate, rngstate);
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
	int onek;
	int tempindex;
	int xoriter;
	int abort;
	double thevalue;
	int thebit;
	unsigned char thebyte;
	int binary_mode;
	int xormode;
	int xorbits;
	int verbose_mode;
	int kilobytes;
	int no_k=1;
	int ofile;
	char filename[1000];
	char jfilename[1000];
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
	double sums_entropy;
	double postxor_entropy;
	double total_entropy;
	double prob;

	int samplenum;
	int simrun;
	unsigned char thesample[256];
	double floatingpointsamples[256];
	t_rngstate rngstate;
	t_modelstate modelstate;


	/* Defaults */
	binary_mode = 0; /* binary when 1, hex when 0 */
	ofile = 0;       /* use stdout instead of outputfile*/
	model = MODEL_PURE;
	xormode = 0;  /* do xor when 1, else don't do xor */
	xorbits = 3;  /* the number of bits to xor together when xormode=1 */
	verbose_mode = 0;
	kilobytes = 1;
	modelstate.using_jfile = 0;
	modelstate.using_infile = 0;
	input_format = INFORMAT_HEX;
	rngstate.input_format = INFORMAT_HEX;
	rngstate.randseed=0;
	rngstate.rdrand_available=0;
	rngstate.devurandom_available=0;
	using_xor_range=0;
	xmin=0;
	xmax=0;
	aesni_supported = 0;
    
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
	  
	rngstate.c_max=511;
	
	modelstate.left_stepsize = 0.1;
	modelstate.right_stepsize = 0.1;
	modelstate.bias = 0.2;
	modelstate.correlation = 0.1;
	modelstate.mean = 0.0;
	modelstate.variance = 0.1;
	modelstate.using_stepnoise = 0;
	modelstate.stepnoise = 0.0;
	modelstate.sums_bias = 0.0;
	
	modelstate.xorshift_size=32;
	sums_entropy = 0.0;
	postxor_entropy = 0.0;
	total_entropy = 0.0;
    width = 32;
	linewidth = 32;
    
    filename[0] = (char)0;
	jfilename[0] = (char)0;
	infilename[0] = (char)0;
    
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
    
    gotxmin = 0;
    gotxmax = 0;
    char optString[] = "c:m:l:r:B:o:j:i:f:k:w:bxsnvh";
    static const struct option longOpts[] = {
    { "binary", no_argument, NULL, 'b' },
    { "xor", required_argument, NULL, 'x' },
    { "xmin", required_argument, NULL, 0 },
    { "xmax", required_argument, NULL, 0 },
    { "seed", no_argument, NULL, 's' },
    { "noaesni", no_argument, NULL, 'n' },
    { "cmax", required_argument, NULL, 'c' },
    { "model", required_argument, NULL, 'm' },
    { "verbose", required_argument, NULL, 'v' },
    { "left", required_argument, NULL, 'l' },
    { "right", required_argument, NULL, 'r' },
    { "stepnoise", required_argument, NULL, 0 },
    { "bias", required_argument, NULL, 0 },
    { "correlation", required_argument, NULL, 0 },
    { "mean", required_argument, NULL, 0 },
    { "variance", required_argument, NULL, 0 },
    
    { "lcg_a", required_argument, NULL, 0 },
    { "lcg_c", required_argument, NULL, 0 },
    { "lcg_m", required_argument, NULL, 0 },
    { "lcg_truncate", required_argument, NULL, 0 },
    { "lcg_outbits", required_argument, NULL, 0 },
    
    { "pcg_state_size", required_argument, NULL, 0 },
    { "pcg_generator", required_argument, NULL, 0 },
    { "pcg_of", required_argument, NULL, 0 },
    
    { "xorshift_size", required_argument, NULL, 0 },
    
    { "output", required_argument, NULL, 'o' },
    { "jfile", required_argument, NULL, 'j' },
    { "infile", required_argument, NULL, 'i' },
    { "informat", required_argument, NULL, 'f' },
    { "blocks", required_argument, NULL, 'k' },
    { "width", required_argument, NULL, 'w' },
    { "help", no_argument, NULL, 'h' },
    { NULL, no_argument, NULL, 0 }
    };

    opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
    while( opt != -1 ) {
        switch( opt ) {
            case 'b':
                binary_mode = 1;
                break;
                
            case 'x':
                xormode = 1;
                xorbits = atoi(optarg);
                break;
                
            case 's':
                rngstate.randseed = 1;
                break;
                
            case 'n':
                aesni_supported = 0;
                break;
            
            case 'v':
                verbose_mode = 1;
                break;
            
            case 'l':
                modelstate.left_stepsize = atof(optarg);
                break;
                
            case 'r':
                modelstate.right_stepsize = atof(optarg);
                break;
                
            case 'o':
                ofile = 1;
                strcpy(filename,optarg);
                break;
                
            case 'j':
                modelstate.using_jfile = 1;
                strcpy(jfilename,optarg);
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
                    printf("Input file format %s not recognized. Choose from binary, hex or 01.\n",optarg);
                    exit(1);
                }
                rngstate.input_format = input_format;
                break; 

            case 'm':
                if (strcmp(optarg,"sums")==0) model=MODEL_SUMS;
                else if (strcmp(optarg,"pure")==0) model=MODEL_PURE;
                else if (strcmp(optarg,"biased")==0) model=MODEL_BIASED;
                else if (strcmp(optarg,"correlated")==0) model=MODEL_CORRELATED;
                else if (strcmp(optarg,"lcg")==0) model=MODEL_LCG;
                else if (strcmp(optarg,"pcg")==0) model=MODEL_PCG;
                else if (strcmp(optarg,"xorshift")==0) model=MODEL_XORSHIFT;
                else if (strcmp(optarg,"normal")==0) model=MODEL_NORMAL;
                else if (strcmp(optarg,"file")==0) model=MODEL_FILE;
                else
                {
                    printf("model type %s not recognized. Choose from sums, pure, biased, correlated, normal or file.\n",optarg);
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
               
            case 0:     /* long option without a short arg */
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
                if( strcmp( "bias", longOpts[longIndex].name ) == 0 ) {
                    modelstate.bias = atof(optarg);
                }
                if( strcmp( "correlation", longOpts[longIndex].name ) == 0 ) {
                    modelstate.correlation = atof(optarg);
                }
                if( strcmp( "mean", longOpts[longIndex].name ) == 0 ) {
                    modelstate.mean = atof(optarg);
                }
                if( strcmp( "variance", longOpts[longIndex].name ) == 0 ) {
                    modelstate.variance = atof(optarg);
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

	/* Range check the var args */

	abort = 0;
	if (using_xor_range == 1)
	{
		if (gotxmin != gotxmax)
		{
			printf("Error: -y and -z (xor_min and xor_max) must be used together or not at all\n");
			abort = 1;
		}
		if (xmin == xmax)
		{
			printf("Error: -y and -z (xor_min and xor_max) must be different and positive\n");
			abort = 1;
		}
		if (xmin < 1)
		{
			printf("Error: -y (xor min) must be 1 or greater. Provided value = %d\n",xmin);
			abort = 1;
		}
		if (xmax < 1)
		{
			printf("Error: -z (xor max) must be 1 or greater. Provided value = %d\n",xmax);
			abort = 1;
		}

		if (xormode==1)
		{
			printf("Error: -x (xor ratio) cannot be used with the XOR range (-y, -z) options");
			abort = 1;
		}

	}
	if (rngstate.c_max < 1)
	{
        	printf("Error: -c n: cmax must be an integer of 1 or greater. Supplied value = %d\n",rngstate.c_max);
        	abort=1;
	}

	if (xorbits < 1)
	{
        	printf("Error: -x n: XOR ratio out of bounds. n must be an integer of 1 or greater. Supplied value = %d\n",xorbits);
        	abort=1;
        }

	if (modelstate.left_stepsize == 0.0)
	{
        	printf("Error: -l n: left stepsize cannot be 0.0. Supplied value = %f\n",modelstate.left_stepsize);
        	abort=1;
        }
	if (modelstate.right_stepsize == 0.0)
	{
        	printf("Error: -l n: right stepsize cannot be 0.0 . Supplied value = %f\n",modelstate.right_stepsize);
        	abort=1;
        }
	if (modelstate.left_stepsize < 0.0)
	{
        	printf("Error: -l n: left stepsize cannot be negative. Supplied value = %f\n",modelstate.left_stepsize);
        	abort=1;
        }
	if (modelstate.right_stepsize < 0.0)
	{
        	printf("Error: -l n: right stepsize cannot be negative . Supplied value = %f\n",modelstate.right_stepsize);
        	abort=1;
        }
	if (modelstate.stepnoise < 0.0)
	{
        	printf("Error: --stepnoise n: stepnoise cannot be negative . Supplied value = %f\n",modelstate.stepnoise);
        	abort=1;
	}
	if ((modelstate.bias < 0.0) || (modelstate.bias >1.0))
	{
		printf("Error: --bias n: bias must be between 0.0 and 1.0. Supplied value = %f\n",modelstate.bias);
		abort=1;
	}

	if ((modelstate.correlation < -1.0) || (modelstate.correlation >1.0))
	{
		printf("Error: --correlation n: correlation must be between -1.0 and 1.0. Supplied value = %f\n",modelstate.correlation);
		abort=1;
	}

	if (kilobytes < 1)
	{
        	printf("Error: -k n: Output size must be 1 or more kilobytes. n must be an integer of 1 or greater. Supplied value = %d\n",kilobytes);
        	abort=1;
        }
	if ((width >256) || (width < 1))
	{
		printf("Error: Width must be from 1 to 256\n");
		abort = 1;	
	} 
	if (model==MODEL_PCG) {
	    if ((modelstate.pcg_state_size != 16) && (modelstate.pcg_state_size != 32)
	          && (modelstate.pcg_state_size != 64) && (modelstate.pcg_state_size != 128)) {
	        printf("Error: A pcg_size must be one of 16, 32, 64 or 128\n");
	        exit(1);
	    }    
	}
	if (model==MODEL_XORSHIFT) {
	    if ((modelstate.xorshift_size != 32) &&  (modelstate.xorshift_size != 128)) {
	        printf("Error: A xorshift_size must be one of 32 or 128\n");
	        exit(1);
	    }    
	}
	if (model==MODEL_FILE) {
		if (modelstate.using_infile == 0) {
			printf("Error: A file must be provided for the file input model using -i <filename> or --infile=<filename>\n");
		abort = 1;
		} 
	}
	if (abort==1) {
		exit(1);
	}

    /* Test for nondeterministic random source */
    
	if (rngstate.randseed==1) {
        if (rdrand_check_support()==1) {    
            rngstate.rdrand_available=1;	
        }
        else if ((rngstate.devrandom =  fopen("/dev/urandom", "r")) !=NULL) {
            rngstate.devurandom_available = 1;
        }
        else {
            rngstate.rdrand_available=0;
		    rngstate.devurandom_available = 0;
		    fprintf(stderr,"Neither /dev/urandom not RdRand Supported for nondeterministic seeding.");
            exit(1);
		}
	}
	
	/* Print out the job parameters */
	if (verbose_mode==1)
	{
        if (aesni_check_support() == 1)
            printf("AESNI Supported in instruction set\n");
        else
            printf("AESNI Not supported in instruction set\n");
 
		if (binary_mode == 0)
			printf("Format=Hex\n");
		else
			printf("Format=Binary\n");

		if (model == MODEL_SUMS)
		{
			printf("model=sums\n");
			printf("  left stepsize  = %f\n",modelstate.left_stepsize);
			printf("  right stepsize = %f\n",modelstate.right_stepsize);
			if (modelstate.using_stepnoise==0)
			{
				printf("  Not adding noise to stepsize\n");
			}
			else
				printf("  Adding noise of variance %f to stepsize\n",modelstate.stepnoise);
		}

		if (model == MODEL_PURE)
		{
			printf("model=pure\n");
		}

		if (model == MODEL_BIASED)
		{
			printf("model=biased\n");
			printf("  bias  = %f\n",modelstate.bias);
		}

		if (model == MODEL_CORRELATED)
		{
			printf("model=correlated\n");
			printf("  correlation  = %f\n",modelstate.correlation);
		}

		if (model == MODEL_LCG)
		{
			printf("model=linear congruential generator\n");
			printf("  a  = 0x%llx\n",modelstate.lcg_a);
			printf("  c  = 0x%llx\n",modelstate.lcg_c);
			printf("  m  = 0x%llx\n",modelstate.lcg_m);
			printf("  start x = 0x%llx\n",modelstate.lcg_x % modelstate.lcg_m);
			printf("  Output bit field = %d:%d\n",
			        (modelstate.lcg_truncate)+(modelstate.lcg_outbits)-1,modelstate.lcg_truncate);
		}
		if (model == MODEL_PCG)
		{
			printf("model=permuted congruential generator\n");
			printf("  state size = %d\n",modelstate.pcg_state_size);
			if (modelstate.pcg_alg == PCG_MCG)
			    printf("  State update algorithm MCG\n");
			else if (modelstate.pcg_alg == PCG_LCG)
			    printf("  State update algorithm LCG\n");
			else
			    printf("  Unknown State Update Function %d\n", modelstate.pcg_alg);
			        
			if (modelstate.pcg_of == XSH_RS)
			    printf("  Output function XSH_RS\n");
			else if (modelstate.pcg_of == XSH_RR)
			    printf("  Output function XSH_RR\n");
			else
			    printf("  Unknown Output function %d\n", modelstate.pcg_of);
		}
		if (model == MODEL_XORSHIFT)	
		{
			printf("model=XORSHIFT\n");
		}	
		if (model == MODEL_NORMAL)
		{
			printf("model=normal\n");
			printf("  mean  = %f\n",modelstate.mean);
			printf("  variance  = %f\n",modelstate.variance);
		}

		if (model == MODEL_FILE)
		{
			printf("model=file\n");
			printf("  filename  = %s\n",infilename);
			if (rngstate.input_format == INFORMAT_HEX)
			    printf("  File input format Hex\n");
			else if (rngstate.input_format == INFORMAT_BINARY)
			    printf("  File input format Binary\n");
			else if (rngstate.input_format == INFORMAT_01 )
			    printf("  File input format ASCII Binary\n");
			else
			    printf("  Unknown input format :%d\n",rngstate.input_format);
		}

		printf("size = %d kilobytes\n", kilobytes);

		if (xormode == 1)
			printf("XOR mode on, fixed ratio=%d:1\n",xorbits);
		else
			printf("XOR mode off\n");

		if (using_xor_range == 1)
			printf("XOR range mode on, ratios between %d:1 and %d:1, chosen randomly\n",xmin, xmax);
		else
			printf("XOR range mode off\n");

		if (ofile == 0)
			printf("Output to STDOUT\n");
		else
			printf("Output to file %s\n",filename);

		if (rngstate.randseed==1)
		{
			printf("Hardware Random Seeding on. Non deterministic mode\n");
			printf("  Reseed c_max=%d\n",rngstate.c_max);
			if (rngstate.rdrand_available) {
			    printf("  Using RdRand as nondeterministic source\n");
			}
			else if(rngstate.devurandom_available) {
			    printf("  Using /dev/urandom as nondeterministic source");
			}   
		}
		else
		{
			printf("Hardware Random Seeding off. Deterministic mode.\n");
			printf("  Restir c_max=%d\n",rngstate.c_max);
		}
		
		if (modelstate.using_jfile==1)
		{
			printf("Outputting internal per bit bias to file %s\n",jfilename);
		}
	}

	/* open the output file if needed */
	FILE *fp;

	if (ofile==1)
	{
		fp = fopen(filename, "wb");
		if (fp == NULL) {
			perror("failed to open output file for writing");
			exit(1);
		}
	}

	/* open the j file if needed */
	if (modelstate.using_jfile==1)
	{
		modelstate.jfile = fopen(jfilename, "wb");
		if (modelstate.jfile == NULL) {
			perror("failed to open output j file for writing");
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
			perror("failed to open input file for reading");
			exit(1);
		}
		else
		{
			if (verbose_mode ==1)
				printf("opened input file %s for reading\n",infilename);
		}
	}

	/* Initialize the RNG */
	initialize_sim(model, &modelstate, &rngstate);

	/* For each stepsize, perform the simulation over 256 samples */
	/* And do it 4 times */

	/* entropy = 0x00;*/
	/* Pull some bits to let it settle. */
	if (!((model==MODEL_FILE) || (model==MODEL_NORMAL) || (model==MODEL_LCG) || (model==MODEL_PCG)))
	for(i=0; i<128; i++)
	{
		thebit = entropysource(model, &modelstate, &rngstate);
		modelstate.lastbit = thebit;
	}

	/* Start with pulling samples and testing them */

	if (model==MODEL_NORMAL) /* or any other floating point model that is added */
	{
		
		for (simrun =0; simrun < kilobytes; simrun++)
		{
			for (onek=0;onek<4;onek++)
			{
				for(samplenum=0;samplenum<256;samplenum++)
				{
					thevalue = dbl_entropysource(model, &modelstate, &rngstate);
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
						printf("%0.8f\n",floatingpointsamples[j]);
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
	}
	else /* model = sums, pure, biased, correlated, lcg, pcg, xorshift or file*/
	{
        lineindex = 0;
        
		for (simrun =0; ((simrun < kilobytes) || ((model==MODEL_FILE) && (no_k==1))); simrun++)
		{
			for (onek=0;onek<4;onek++)
			{
				for(samplenum=0;samplenum<256;samplenum++)
				{
					thebyte = (unsigned char)0x00;

					/* Pull 8 bits */
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
								prob = modelstate.sums_bias;
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
								prob = modelstate.sums_bias;
								sums_entropy = bias2entropy(prob);
							}
						}
						else /* no xorring */
						{
							thebit = entropysource(model, &modelstate, &rngstate);
							if (rngstate.reached_eof == 1) {
							    if ((samplenum > 0) && (samplenum < 256)) {
							        /*printf("reached EOF with samplenum > 0 = %d\n",samplenum);*/
							        goto eof_with_partial_block;
							    }
							    else {
							        /*printf("Going to EOF with full block. samplenum = %d\n",samplenum);*/
							        goto reached_eof;
							    }
							}
							modelstate.lastbit = thebit;
							prob = modelstate.sums_bias;
							sums_entropy = bias2entropy(prob);
							total_entropy += sums_entropy;
						}

						if ((thebit & 0x01)==1)
						{
							thebyte = (thebyte << 1) | 0x01;
						}
						else
						{
							thebyte = (thebyte << 1) & 0xfe;
						}
					}

					thesample[samplenum]=thebyte;
				}

				/* Output the 256 byte block */
                eof_with_partial_block:

				if (ofile == 1 && binary_mode==1) /* binary to a file */
				{
					fwrite(thesample, samplenum, 1, fp);
				}
				else if (ofile == 0 && binary_mode == 0) /* hex to stdout */
				{
					tempindex = 0;
                    /* printf(" Samplenum %d\n",samplenum);*/
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
                    if (lineindex != 0) printf("\n");
				}
				else /* binary to stdout */
				{
					fwrite(thesample, samplenum, 1, stdout);
				}
			}
		}
		reached_eof:
		
		if (verbose_mode ==1)
		{
			printf("Total Entropy = %F\n",total_entropy);
			printf("Per bit Entropy = %F %% \n",(100.0*(total_entropy/(8.0*kilobytes*1024.0))));
		}
	}
	return 0;

}


