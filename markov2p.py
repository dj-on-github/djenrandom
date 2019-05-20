#!/usr/bin/env python

from __future__ import print_function
from __future__ import division

import math
import random
import sys
import gmpy2
from gmpy2 import mpfr

verbose_mode = 1 

def biasscc_2_p(bias,scc):
    p01 = bias * (mpfr('1.0') - scc)
    p10 = (mpfr('1.0')-bias)*(mpfr('1.0')-scc)
    return p01,p10

# A library for coverting between points, scc, bias and entropy
# with the 2 parameter markov model.

def print_symbol(x, bitwidth):
    msymboltext = ""
    for i in range(bitwidth):
        if (((x >> (bitwidth-1-i)) & 0x01)==0):
            msymboltext += "0"
        else:
            msymboltext += "1"
    return msymboltext


# Compute the min entropy per symbol for the
# markov 2 parameter model, given the markov model
# parameters p01 and p10.
def symbol_prob(p01, p10, x, bitwidth):
    plist0=mpfr('1.0')
    plist1=mpfr('1.0')
    
    p00 = mpfr('1.0')-p01
    p11 = mpfr('1.0')-p10
    mu = p01/(p10+p01)
    p0 = mpfr('1.0')-mu
    p1 = mu
    
    symboltext = print_symbol(x,bitwidth)
     
    if ((p01==0.5) and (p10==0.5)):
        return mpfr('1.0')
    
    plist0 = mpfr('1.0')
    plist1 = mpfr('1.0')
    
    if ((x>>(bitwidth-1) & 0x1)==0):
        plist0 *= p00
        plist1 *= p10
    else:
        plist0 *= p01
        plist1 *= p11
    
    #printf(stderr," plist0=%f  ",plist0);
    #printf(stderr," plist1=%f\n",plist1);
    
    for i in range(bitwidth-2):
        bp = ((x >> (bitwidth-2-i)) & 0x3) #Get the bit pair
        #printf(stderr,"       bitpair %d = %d ",i,bp);
        if (bp==0):
            plist0 *= p00
            plist1 *= p00
            #printf(stderr," plist0=%f * p00(%f)  ",plist0,p00);
            #printf(stderr," plist1=%f * p00(%f)\n",plist1,p00);
        elif (bp==1):
            plist0 *= p01
            plist1 *= p01
            printf(stderr," plist0=%f * p01(%f)  ",plist0,p01);
            printf(stderr," plist1=%f * p01(%f)\n",plist1,p01);
        elif (bp==2):
            plist0 *= p10
            plist1 *= p10
            #printf(stderr," plist0=%f * p10(%f)  ",plist0,p10);
            #printf(stderr," plist1=%f * p10(%f)\n",plist1,p10);
        elif (bp==3):
            plist0 *= p11
            plist1 *= p11
            #printf(stderr," plist0=%f * p11(%f)  ",plist0,p11);
            #printf(stderr," plist1=%f * p11(%f)\n",plist1,p11);
    
    p = (p0 * plist0) + (p1 * plist1)
    return p

def max(x, y):
    if (x>y):
         return x
    if (y>x):
        return y
    return x

def mk_symbol(prefix, tbp, postfix, bitwidth):
    rep = (bitwidth-2)/2
    pattern = prefix
    
    for i in range(rep):
        pattern = (pattern << 2) + tbp
    pattern = (pattern << 1) + postfix
    
    return pattern

def mk_symbol_nopostfix(prefix, tbp, bitwidth):
    pattern = prefix;
    for i in range((bitwidth-1)//2):
        pattern = (pattern << 2) + tbp
    
    return pattern

def most_probable_transition_pair(p01, p10):
    mu = p01/(p10+p01)
    p0 = mpfr('1.0')-mu
    p1 = mu
    
    p00 = mpfr('1.0') - p01
    p11 = mpfr('1.0') - p10
        
    p010 = p0 * p01 * p10
    p101 = p1 * p10 * p01
    p000 = p0 * p00 * p00
    p111 = p1 * p11 * p11
    
    if      ((p111 >= p000) and (p111 >= p101) and (p111 >= p010)):
        return "P111_MAX"
    elif ((p000 >= p111) and (p000 >= p101) and (p000 >= p010)):
        return "P000_MAX"
    elif ((p101 >= p111) and (p101 >= p000) and (p101 >= p010)):
        return "P101_MAX"
    elif ((p010 >= p111) and (p010 >= p000) and (p010 >= p101)):
            return "P010_MAX"
    
    return "EQUIPROBABLE"

def most_probable_symbol_odd(p01, p10,bitwidth):
    if (most_probable_transition_pair(p01, p10) == "P000_MAX"):
        mps = 0
    elif (most_probable_transition_pair(p01, p10) == "P111_MAX"):
        for i in range ((bitwidth-1)>>1): #(i=0; i<((bitwidth-1)>>1); i++) {
            mps = mps << 2
            mps = mps + 3
        mps = mps << 1;
        mps = mps + 1;
    elif (most_probable_transition_pair(p01, p10) == "P010_MAX"):
        for i in range ((bitwidth-1)>>1): #(i=0; i<((bitwidth-1)>>1); i++) {
            mps = mps << 2
            mps = mps + 1
        mps = mps << 1
        mps = mps + 0
    elif (most_probable_transition_pair(p01, p10) == "P101_MAX"):
        for i in range ((bitwidth-1)>>1): #(i=0; i<((bitwidth-1)>>1); i++) {
            mps = mps << 2
            mps = mps + 2
        mps = mps << 1
        mps = mps + 1
    else: #     // Equiprobable case, any value will do.
        mps = 0
    return mps

def most_probable_symbol_even(p01, p10,bitwidth):
    p00 = mpfr('1.0') - p01
    p11 = mpfr('1.0') - p10
    
    mps = 0
        
    if (most_probable_transition_pair(p01, p10) == "P000_MAX"):
        mps = 0
    elif (most_probable_transition_pair(p01, p10) == "P111_MAX"):
        for i in range (bitwidth>>1):
            mps = mps << 2
            mps = mps + 3
    elif (most_probable_transition_pair(p01, p10) == "P010_MAX"):
        for i in range ((bitwidth-2)>>1):
            mps = mps << 2
            mps = mps + 1
        mps = mps << 2
        if (p01 > p00):
            mps = mps + 1
        else:
            mps = mps + 0
    elif (most_probable_transition_pair(p01, p10) == "P101_MAX"):
        for i in range ((bitwidth-2)>>1):
            mps = mps << 2
            mps = mps + 2
        mps = mps << 2
        if (p11 > p10):
            mps = mps + 3
        else:
            mps = mps + 2
    else: #     // Equiprobable case, any value will do.
        mps = 0
    return mps

def most_probable_symbol(p01, p10,bitwidth):
    if ((bitwidth & 0x01)==0x01):
        mps = most_probable_symbol_odd(p01,p10,bitwidth)
    else:
        mps = most_probable_symbol_even(p01,p10,bitwidth)
    
    
    #if (verbose_mode) print("   MCV = 0x%" PRIx64 " \n" % mps,file=stderr)
    return mps
    
def symbol_max_probability(p01, p10,bitwidth):

    mu = p01/(p10+p01)
    p0 = mpfr('1.0')-mu
    p1 = mu
    
    p00 = mpfr('1.0') - p01
    p11 = mpfr('1.0') - p10
    
    mps = most_probable_symbol(p01,p10,bitwidth)
    mcv = mps
    
    # unpack the symbol bits into an array of bits
    bits = [0,]   # first with x[-1]=0
    for i in range(bitwidth):
        bits.append((mps >> (bitwidth-1-i)) & 0x01)
    
    # Compute the symbol probability by going through the
    # bits and multiplying the transition probabilities.
    p_0mps = mpfr('1.0')
    #if (verbose_mode) fprintf(stderr,"   Prob = mpfr('1.0')");
    for i in range(bitwidth):
        if      ((bits[i]==0) and (bits[i+1]==0)):
            p_0mps = p_0mps * p00
            #if (verbose_mode) fprintf(stderr, " * P00");
        elif ((bits[i]==0) and (bits[i+1]==1)):
            p_0mps = p_0mps * p01
            #if (verbose_mode) fprintf(stderr, " * P01");
        elif ((bits[i]==1) and (bits[i+1]==0)):
            p_0mps = p_0mps * p10
            #if (verbose_mode) fprintf(stderr, " * P10");
        elif ((bits[i]==1) and (bits[i+1]==1)):
            p_0mps = p_0mps * p11
            #if (verbose_mode) fprintf(stderr, " * P11");
    #if (verbose_mode) fprintf(stderr,"\n");

    
    bits[0] = 1 #   // then with x[-1]=1
    
    p_1mps = mpfr('1.0')
    for i in range(bitwidth):
        if ((bits[i]==0) and (bits[i+1]==0)):
            p_1mps = p_1mps * p00
        elif ((bits[i]==0) and (bits[i+1]==1)):
            p_1mps = p_1mps * p01
        elif ((bits[i]==1) and (bits[i+1]==0)):
            p_1mps = p_1mps * p10
        elif ((bits[i]==1) and (bits[i+1]==1)):
            p_1mps = p_1mps * p11;
    p_mps = (p0 * p_0mps) + (p1 * p_1mps)
    return p_mps,mcv
    
    
def p_to_entropy(p01, p10,bitwidth):
    smp = mpfr('0.0')
    
    mcv_prob,mcv = symbol_max_probability(p01, p10, bitwidth)
    
    ent = -gmpy2.log2(mcv_prob)
    return ent/bitwidth, mcv_prob, mcv;
    
def near(x,y, epsilon):
    return ((y > x-epsilon) and (y<x+epsilon))

def pick_point(desired, epsilon, bitwidth):
    while True:
        chosen_param = random.genrandbits(16) & 0x01
        chosen_side = random.genrandbits(16) & 0x01
        
        if (chosen_param==0):
            p01 = chosen_side
            p10 = random.random()
        else:
            p10 = chosen_side
            p01 = random.random()

        edge_entropy,mcv_prob,mcv=p_to_entropy(p01, p10, bitwidth)
        if edge_entropy > desired:
            break
    
    startpoint01 = 0.5
    startpoint10 = 0.5
    endpoint01 = p01
    endpoint10 = p10
    
    choice01 = (startpoint01 + endpoint01)/2.0
    choice10 = (startpoint10 + endpoint10)/2.0
    Hc,mcv_prob,mcv = p_to_entropy(choice01, choice10, bitwidth)
    
    while (not near(Hc, desired, epsilon)):
        if (Hc > desired):
            startpoint01 = choice01
            startpoint10 = choice10
        else:
            endpoint01 = choice01
            endpoint10 = choice10
        choice01 = (startpoint01 + endpoint01)/2.0
        choice10 = (startpoint10 + endpoint10)/2.0
        
    p01 = choice01
    p10 = choice10
    return p01,p10

if __name__ == '__main__':
    gmpy2.get_context().precision=16384
    bias = float(sys.argv[1])
    scc = float(sys.argv[2])
    bits = int(sys.argv[3])
    p01,p10 = biasscc_2_p(bias,scc)

    entropy, mcv_p, mcv =p_to_entropy(p01,p10,bits)
    print("entropy per bit =",entropy)
    print("mcv prob        =",mcv_p)
    print("mcv             = %x" % mcv)

    
