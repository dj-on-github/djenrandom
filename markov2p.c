
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

#include <unistd.h>
#include <string.h>

#include "djenrandommodel.h"
#include "markov2p.h"

#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"


// A library for coverting between points, scc, bias and entropy
// with the 2 parameter markov model.

extern int verbose_mode;

char msymboltext[255];

void print_symbol(uint64_t x, int bitwidth) {
    int i;

    for(i=0;i<bitwidth;i++) {
        if (((x >> (bitwidth-1-i)) & 0x01)==0) msymboltext[i]='0';
        else msymboltext[i]='1';
    }
    msymboltext[bitwidth]=(char)0;
}



// Compute the min entropy per symbol for the
// markov 2 parameter model, given the markov model
// parameters p01 and p10.
double symbol_prob(double p01, double p10, uint64_t x, int bitwidth) {
    double p00;
    double p11;
    double mu;
    double p0;
    double p1;
    double plist0;
    double plist1;
    int bp;
    double p;
    
    int i;
    
    plist0=1.0;
    plist1=1.0;
    
    p00 = 1.0-p01;
    p11 = 1.0-p10;
    mu = p01/(p10+p01);
    p0 = 1.0-mu;
    p1 = mu;
    
    print_symbol(x,bitwidth);
    //fprintf(stderr,"  SYMBOL PROB p01=%f,   p10=%f,  x=%" PRIx64 " = b%s  bitwidth=%d\n",p01,p10,x,symboltext,bitwidth);
    //fprintf(stderr,"              P01 = %f\n", p01);
    //fprintf(stderr,"              P10 = %f\n", p10);
    //fprintf(stderr,"              P00 = %f\n", p00);
    //fprintf(stderr,"              P11 = %f\n", p11);
    //fprintf(stderr,"              mu = %f\n", mu);
    //fprintf(stderr,"              P0 = %f\n", p0);
    //fprintf(stderr,"              P1 = %f\n", p1);
     
    if ((p01==0.5) && (p10==0.5)) return 1.0;
    
    plist0 = 1.0;
    plist1 = 1.0;
    
    if ((x>>(bitwidth-1) & 0x1)==0) {
        plist0 *= p00;
        plist1 *= p10;
    }
    else {
        plist0 *= p01;
        plist1 *= p11;
    }
    
    //fprintf(stderr," plist0=%f  ",plist0);
    //fprintf(stderr," plist1=%f\n",plist1);
    
    for (i=0;i<(bitwidth-2);i++) {
        bp = ((x >> (bitwidth-2-i)) & 0x3);  // Get the bit pair
        //fprintf(stderr,"       bitpair %d = %d ",i,bp);
        if (bp==0) {
            plist0 *= p00;
            plist1 *= p00;
            //fprintf(stderr," plist0=%f * p00(%f)  ",plist0,p00);
            //fprintf(stderr," plist1=%f * p00(%f)\n",plist1,p00);
        } else if (bp==1) {
            plist0 *= p01;
            plist1 *= p01;
            //fprintf(stderr," plist0=%f * p01(%f)  ",plist0,p01);
            //fprintf(stderr," plist1=%f * p01(%f)\n",plist1,p01);
        } else if (bp==2) {
            plist0 *= p10;
            plist1 *= p10;
            //fprintf(stderr," plist0=%f * p10(%f)  ",plist0,p10);
            //fprintf(stderr," plist1=%f * p10(%f)\n",plist1,p10);
        } else if (bp==3) {
            plist0 *= p11;
            plist1 *= p11;
            //fprintf(stderr," plist0=%f * p11(%f)  ",plist0,p11);
            //fprintf(stderr," plist1=%f * p11(%f)\n",plist1,p11);
        }
    }
    
    p = (p0 * plist0) + (p1 * plist1);
    
    //fflush(stdout);
    return p;
    
}

double max(double x, double y) {
    if (x>y) return x;
    if (y>x) return y;
    return x;
}

uint64_t mk_symbol(int prefix, int tbp, int postfix, int bitwidth) {
    int rep;
    int i;
    
    uint64_t pattern;
    
    rep = (bitwidth-2)/2;
    pattern = prefix;
    
    for(i=0;i<rep;i++) {
        pattern = (pattern << 2) + tbp; 
    }
    pattern = (pattern << 1) + postfix;
    
    return pattern;    
}

uint64_t mk_symbol_nopostfix(int prefix, int tbp, int bitwidth) {
    //int rep;
    int i;
    
    uint64_t pattern;
    
    //rep = (bitwidth-2)/2;
    pattern = prefix;

    pattern = prefix;
    for(i=0;i<((bitwidth-1)/2);i++) {
        pattern = (pattern << 2) + tbp; 
    }
    
    return pattern;    
}

int most_probable_transition_pair(double p01, double p10) {
    double p010;
    double p101;
    double p000;
    double p111;
    double p00;
    double p11;
    double p0;
    double p1;
    
    double mu;

    mu = p01/(p10+p01);
    p0 = 1.0-mu;
    p1 = mu;
    
    p00 = 1.0 - p01;
    p11 = 1.0 - p10;
        
    p010 = p0 * p01 * p10;
    p101 = p1 * p10 * p01;
    p000 = p0 * p00 * p00;
    p111 = p1 * p11 * p11;
    
    if      ((p111 >= p000) && (p111 >= p101) && (p111 >= p010)) {
            return P111_MAX;
    }
    else if ((p000 >= p111) && (p000 >= p101) && (p000 >= p010)) {
            return P000_MAX;
    }
    else if ((p101 >= p111) && (p101 >= p000) && (p101 >= p010)) {
            return P101_MAX;
    }
    else if ((p010 >= p111) && (p010 >= p000) && (p010 >= p101)) {
            return P010_MAX;
    }
    
    return EQUIPROBABLE;

}

uint64_t most_probable_symbol_odd(double p01, double p10,int bitwidth) {
    uint64_t mps;
    int i;
        
    if (most_probable_transition_pair(p01, p10) == P000_MAX) {
        mps = 0;
    } else if (most_probable_transition_pair(p01, p10) == P111_MAX) {
        for (i=0; i<((bitwidth-1)>>1); i++) {
            mps = mps << 2;
            mps = mps + 3;
        }
        mps = mps << 1;
        mps = mps + 1;
    } else if (most_probable_transition_pair(p01, p10) == P010_MAX) {
        for (i=0; i<((bitwidth-1)>>1); i++) {
            mps = mps << 2;
            mps = mps + 1;
        }
        mps = mps << 1;
        mps = mps + 0;
    } else if (most_probable_transition_pair(p01, p10) == P101_MAX) {
        for (i=0; i<((bitwidth-1)>>1); i++) {
            mps = mps << 2;
            mps = mps + 2;
        }
        mps = mps << 1;
        mps = mps + 1;
    } else {     // Equiprobable case, any value will do.
        mps = 0;
    }
    return mps;
}

uint64_t most_probable_symbol_even(double p01, double p10,int bitwidth) {
    uint64_t mps;
    int i;
    double p00;
    double p11;
    //double p0;
    //double p1;
    
    //double mu;

    //mu = p01/(p10+p01);
    //p0 = 1.0-mu;
    //p1 = mu;
    
    p00 = 1.0 - p01;
    p11 = 1.0 - p10;
    
    mps = 0;
        
    if (most_probable_transition_pair(p01, p10) == P000_MAX) {
        mps = 0;
    } else if (most_probable_transition_pair(p01, p10) == P111_MAX) {
        for (i=0; i<(bitwidth >> 1); i++) {
            mps = mps << 2;
            mps = mps + 3;
        }
    } else if (most_probable_transition_pair(p01, p10) == P010_MAX) {
        for (i=0; i<((bitwidth-2) >> 1); i++) {
            mps = mps << 2;
            mps = mps + 1;
        }
        mps = mps << 2;
        if (p01 > p00) {
            mps = mps + 1;
        } else {
            mps = mps + 0;
        }
        
    } else if (most_probable_transition_pair(p01, p10) == P101_MAX) {
        for (i=0; i<((bitwidth-2) >> 1); i++) {
            mps = mps << 2;
            mps = mps + 2;
        }
        mps = mps << 2;
        if (p11 > p10) {
            mps = mps + 3;
        } else {
            mps = mps + 2;
        }
    } else {     // Equiprobable case, any value will do.
        mps = 0;
    }
    return mps;
}

uint64_t most_probable_symbol(double p01, double p10,int bitwidth) {
    uint64_t mps;
    
    if ((bitwidth & 0x01)==0x01)
        mps = most_probable_symbol_odd(p01,p10,bitwidth);
    else
        mps = most_probable_symbol_even(p01,p10,bitwidth);
    
    
    if (verbose_mode>1) fprintf(stderr,"   MCV = 0x%" PRIx64 " \n",mps);
    return mps;
    
}

double symbol_max_probability(double p01, double p10,int bitwidth,uint64_t *mcv) {
    double mu;
    double p00;
    double p11;
    double p0;
    double p1;
    uint64_t mps;
    
    double p_0mps;
    double p_1mps;
    double p_mps;
    
    int bits[65];
    int i;
    int j;
    
    for (i=0;i<65;i++) bits[i] = 0;
    
    mu = p01/(p10+p01);
    p0 = 1.0-mu;
    p1 = mu;
    
    p00 = 1.0 - p01;
    p11 = 1.0 - p10;
    
    mps = most_probable_symbol(p01,p10,bitwidth);
    *mcv = mps;
    
    // unpack the symbol bits into an array of bits
    bits[0] = 0;   // first with x[-1]=0
    for (i=0; i<bitwidth; i++) {
        bits[i+1] = (mps >> (bitwidth-1-i)) & 0x01;
    }
    
    if (verbose_mode>1) {
        fprintf(stderr,"   unrolled bits 0 prefix = ");
        for(j=0;j<(bitwidth+1);j++) {
            fprintf(stderr,"%d",bits[j]);
        }
        fprintf(stderr,"\n");
    }
    
    // Compute the symbol probability by going through the
    // bits and multiplying the transition probabilities.
    p_0mps = 1.0;
    if (verbose_mode>1) fprintf(stderr,"   Prob = 1.0");
    for (i=0;i<bitwidth; i++) {
        if      ((bits[i]==0) && (bits[i+1]==0)) {
            p_0mps = p_0mps * p00;
            if (verbose_mode>1) fprintf(stderr, " * P00");
        }
        else if ((bits[i]==0) && (bits[i+1]==1)) {
            p_0mps = p_0mps * p01;
            if (verbose_mode>1) fprintf(stderr, " * P01");
        }                     
        else if ((bits[i]==1) && (bits[i+1]==0)) {
            p_0mps = p_0mps * p10;
            if (verbose_mode>1) fprintf(stderr, " * P10");
        }
        else if ((bits[i]==1) && (bits[i+1]==1)) {
            p_0mps = p_0mps * p11;
            if (verbose_mode>1) fprintf(stderr, " * P11");
        }      
    }    
    if (verbose_mode>1) fprintf(stderr,"\n");

    
    bits[0] = 1;   // then with x[-1]=1
    
    if (verbose_mode>1) {
        fprintf(stderr,"   unrolled bits 1 prefix = ");
        for(j=0;j<(bitwidth+1);j++) {
            fprintf(stderr,"%d",bits[j]);
        }
        fprintf(stderr,"\n");
    }
    
    p_1mps = 1.0;
    if (verbose_mode>1) fprintf(stderr,"   Prob = 1.0");
    for (i=0;i<bitwidth; i++) {
        if      ((bits[i]==0) && (bits[i+1]==0)) {
            p_1mps = p_1mps * p00;
            if (verbose_mode>1) fprintf(stderr, " * P00");
        }
        else if ((bits[i]==0) && (bits[i+1]==1)) {
            p_1mps = p_1mps * p01;
            if (verbose_mode>1) fprintf(stderr, " * P01");
        }                     
        else if ((bits[i]==1) && (bits[i+1]==0)) {
            p_1mps = p_1mps * p10;
            if (verbose_mode>1) fprintf(stderr, " * P10");
        }
        else if ((bits[i]==1) && (bits[i+1]==1)) {
            p_1mps = p_1mps * p11;
            if (verbose_mode>1) fprintf(stderr, " * P11");
        }      
    }    
    if (verbose_mode>1) fprintf(stderr,"\n");
    
    if (verbose_mode>1) {
        fprintf(stderr,"   %sMCV BITS = ",KRED);
        for (i=0; i<bitwidth;i++) {
            fprintf(stderr,"%d",bits[i+1]);
        }
        fprintf(stderr,"%s\n",KWHT);
    }
        
    p_mps = (p0 * p_0mps) + (p1 * p_1mps);
    return p_mps;
}
    
    
double p_to_entropy(double p01, double p10,int bitwidth, double *mcv_prob, uint64_t *mcv) {
    double smp = 0.0;
    double ent;
    uint64_t l_mcv;
    
    smp = symbol_max_probability(p01, p10, bitwidth, &l_mcv);
    *mcv_prob = smp;
    *mcv = l_mcv;
    
    ent = -log2(smp);
    return ent/bitwidth;
}
    
int near(double x,double y, double epsilon) {
    return ((y > x-epsilon) && (y<x+epsilon));
}

void pick_point(double *p01, double *p10, double desired, double epsilon, int bitwidth, t_rngstate* rngstate) {
    int chosen_param;
    int chosen_side;
    int rand1;
    int rand2;
    double startpoint01;
    double startpoint10;
    double endpoint01;
    double endpoint10;
    double choice01;
    double choice10;
    double mcv_prob = -1.0;
    double Hc;
    uint64_t mcv;
    
    double edge_entropy;
    
    do {
        rand1 = getrand16(rngstate);
        rand2 = getrand16(rngstate);
        chosen_param = rand1 & 0x01;
        chosen_side = rand2 & 0x01;
        if (verbose_mode > 1) {
            fprintf(stderr,"      rand1         %04x\n", rand1);
            fprintf(stderr,"      rand2         %04x\n", rand2);
            fprintf(stderr,"      chosen_param  %04x\n", chosen_param);
            fprintf(stderr,"      chosen_side   %04x\n", chosen_side);
        }
        
        if (chosen_param==0) {
            *p01 = (double)chosen_side;
            *p10 = get_rand_double(rngstate);
        }
        else {
            *p10 = (double)chosen_side;
            *p01 = get_rand_double(rngstate);
        }
        edge_entropy=p_to_entropy(*p01, *p10, bitwidth, &mcv_prob, &mcv);
        
    } while (edge_entropy > desired);
    
    startpoint01 = 0.5;
    startpoint10 = 0.5;
    endpoint01 = *p01;
    endpoint10 = *p10;
    
    choice01 = (startpoint01 + endpoint01)/2.0;
    choice10 = (startpoint10 + endpoint10)/2.0;
    Hc = p_to_entropy(choice01, choice10, bitwidth, &mcv_prob, &mcv);
    
    if (verbose_mode > 1) {
    fprintf(stderr,"PICKING for entropy %f\n", desired);
    fprintf(stderr,"                bitwidth  %d\n", bitwidth);
    fprintf(stderr,"      first startpoint01  %f\n", startpoint01);
    fprintf(stderr,"      first startpoint10  %f\n", startpoint10);
    fprintf(stderr,"        first endpoint01  %f\n", endpoint01);
    fprintf(stderr,"        first endpoint10  %f\n", endpoint10);
    fprintf(stderr,"          first mid P01 = %f\n", choice01);
    fprintf(stderr,"          first mid P10 = %f\n", choice10);
    fprintf(stderr,"        start Hc    %f\n", Hc);
    }
    
    fflush(stdout);

    while (!near(Hc, desired, epsilon)) {
        if (verbose_mode>1) fprintf(stderr,"WHILE ...\n");
        if (Hc > desired) {
            startpoint01 = choice01;
            startpoint10 = choice10;
        }
        else {
            endpoint01 = choice01;
            endpoint10 = choice10;
        }
        choice01 = (startpoint01 + endpoint01)/2.0;
        choice10 = (startpoint10 + endpoint10)/2.0;
        
        if (verbose_mode > 1) {
        fprintf(stderr,"          bitwidth  %d\n", bitwidth);       
        fprintf(stderr,"      startpoint01  %f\n", startpoint01);
        fprintf(stderr,"      startpoint10  %f\n", startpoint10);
        fprintf(stderr,"        endpoint01  %f\n", endpoint01);
        fprintf(stderr,"        endpoint10  %f\n", endpoint10);       
        fprintf(stderr,"   mid P01 = %f\n", choice01);
        fprintf(stderr,"   mid P10 = %f\n", choice10);
        }
        Hc = p_to_entropy(choice01,choice10,bitwidth,&mcv_prob, &mcv);
        if (verbose_mode > 1) {
            fprintf(stderr,"   Hc  = %f\n", Hc);
            fprintf(stderr,"   %sMCV Probability = %f%s\n",KCYN,mcv_prob,KWHT);
            fflush(stdout);
        }
    }
    
    if (verbose_mode >1) {
    fprintf(stderr," ** Chose P01 = %f\n", choice01);
    fprintf(stderr," ** Chose P10 = %f\n", choice10);
    }
    *p01 = choice01;
    *p10 = choice10;
    
}


