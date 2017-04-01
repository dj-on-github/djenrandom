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

extern int rdrand16_step(unsigned short int *);
extern int rdseed16_step(unsigned short int *);

extern int rdrand32_step(unsigned int *);
extern int rdseed32_step(unsigned int *);

extern int rdrand64_step(unsigned long long int *);
extern int rdseed64_step(unsigned long long int *);

extern int rdrand_get_uint_retry(unsigned int retry_limit, unsigned int *dest);
extern int rdseed_get_uint_retry(unsigned int retry_limit, unsigned int *dest);
extern int rdrand_get_n_uints_retry(unsigned int n, unsigned int retry_limit, unsigned int *dest);
extern int rdseed_get_n_uints_retry(unsigned int n, unsigned int retry_limit, unsigned int *dest);

extern int rdrand_get_uint_step(unsigned int *dest);
extern int rdseed_get_uint_step(unsigned int *dest);

extern int rdrand_get_n_uints_step(int n, unsigned int *dest);
extern int rdseed_get_n_uints_step(int n, unsigned int *dest);

extern int rdrand_get_bytes_step(unsigned int n, unsigned char *dest);
extern int rdrand_check_support();
extern int rdseed_check_support();
extern int rdrand_get_n_qints_retry(unsigned int n, unsigned int retry_limit, unsigned long long int *dest);
extern int rdseed_get_n_qints_retry(unsigned int n, unsigned int retry_limit, unsigned long long int *dest);

    
extern int aesni_check_support();

