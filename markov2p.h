
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

#define EQUIPROBABLE 0
#define P000_MAX 1
#define P111_MAX 2
#define P101_MAX 3
#define P010_MAX 4

double symbol_prob(double p01, double p10, uint64_t x, int bitwidth) ;
double max(double x, double y) ;
uint64_t mk_symbol(int prefix, int tbp, int postfix, int bitwidth) ;
uint64_t mk_symbol_nopostfix(int prefix, int tbp, int bitwidth) ;
int most_probable_transition_pair(double p01, double p10) ;
uint64_t most_probable_symbol_odd(double p01, double p10,int bitwidth) ;
uint64_t most_probable_symbol_even(double p01, double p10,int bitwidth) ;
uint64_t most_probable_symbol(double p01, double p10,int bitwidth) ;
double symbol_max_probability(double p01, double p10,int bitwidth,uint64_t *mcv) ;
double p_to_entropy(double p01, double p10,int bitwidth, double *mcv_prob, uint64_t *mcv) ;
int near(double x,double y, double epsilon) ;
void pick_point(double *p01, double *p10, double desired, double epsilon, int bitwidth, t_rngstate* rngstate) ;
void make_sample_table(double p01, double p10, int bitwidth, int **sampletable0, int **sampletable1) ;

