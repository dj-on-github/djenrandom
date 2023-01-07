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

#include "rdrand.h"

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "aes128k128d.h"

extern int verbose_mode;

typedef struct {
        unsigned int EAX;
        unsigned int EBX;
        unsigned int ECX;
        unsigned int EDX;
} CPUIDinfo;

void get_cpuid_windows(int leaf, CPUIDinfo *info) {
unsigned int a;
unsigned int b;
unsigned int c;
unsigned int d;

asm("\n\
    mov %4, %%eax\n\
    cpuid\n\
    mov %%eax,%0\n\
    mov %%ebx,%1\n\
    mov %%ecx,%2\n\
    mov %%edx,%3":"=r"(a),"=r"(b),"=r"(c),"=r"(d):"r"(leaf):"%eax","%ebx","%ecx","%edx");
    info->EAX = a;
    info->EBX = b;
    info->ECX = c;
    info->EDX = d;
}

/*void get_cpuid_linux(CPUIDinfo *info, const unsigned int func, const unsigned int subfunc)*/
/*{*/
/*asm(".intel_syntax noprefix\n");*/
/*asm("mov r8, rdi\n");*/
/*asm("mov r9, rsi\n");*/
/*asm("mov r10, rdx\n");*/
/*asm("push rax\n");*/
/*asm("push rbx\n");*/
/*asm("push rcx\n");*/
/*asm("push rdx\n");*/
/*asm("mov eax, r9d\n");*/
/*asm("mov ecx, r10d\n");*/
/*asm("cpuid;\n");*/
/*asm("mov DWORD PTR [r8], eax\n");*/
/*asm("mov DWORD PTR [r8+4], ebx\n");*/
/*asm("mov DWORD PTR [r8+8], ecx\n");*/
/*asm("mov DWORD PTR [r8+12], edx\n");*/
/*asm("pop rdx\n");*/
/*asm("pop rcx\n");*/
/*asm("pop rbx\n");*/
/*asm("pop rax\n");*/
/*asm(".att_syntax prefix\n");*/
/*}*/

/* Trying GAS format to make clang happy*/
void get_cpuid_linux(CPUIDinfo *info, const unsigned int func, const unsigned int subfunc)
{
asm(".intel_syntax noprefix;\n\
mov r8, rdi;\n\
mov r9, rsi;\n\
mov r10, rdx;\n\
push rax;\n\
push rbx;\n\
push rcx;\n\
push rdx;\n\
mov eax, r9d;\n\
mov ecx, r10d;\n\
cpuid;\n\
mov DWORD PTR [r8], eax;\n\
mov DWORD PTR [r8+4], ebx;\n\
mov DWORD PTR [r8+8], ecx;\n\
mov DWORD PTR [r8+12], edx;\n\
pop rdx;\n\
pop rcx;\n\
pop rbx;\n\
pop rax;\n\
.att_syntax prefix\n");
}


void get_cpuid(CPUIDinfo *info, const unsigned int func, const unsigned int subfunc) {
    #if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
        get_cpuid_windows(func, info);
    #else
        get_cpuid_linux(info, func, subfunc);
    #endif
}

typedef unsigned int DWORD;

int check_is_intel() {
    CPUIDinfo info;
   
    get_cpuid(&info,0,0);
    if(memcmp((char *)(&info.EBX), "Genu", 4) == 0 &&
        memcmp((char *)(&info.EDX), "ineI", 4) == 0 &&
        memcmp((char *)(&info.ECX), "ntel", 4) == 0) {
            if (verbose_mode) fprintf(stderr,"Is Intel CPU\n");
            return 1;
    }
    if (verbose_mode) fprintf(stderr,"Is not Intel CPU\n");
    return 0;
}

int check_is_amd() {
    CPUIDinfo info;
   
    get_cpuid(&info,0,0);

    if( memcmp((char *)(&info.EBX), "Auth", 4) == 0 &&
        memcmp((char *)(&info.EDX), "enti", 4) == 0 &&
        memcmp((char *)(&info.ECX), "cAMD", 4) == 0) {
            if (verbose_mode) fprintf(stderr,"Is AMD CPU\n");
            return 1;
    }
    if (verbose_mode) fprintf(stderr,"Is not AMD CPU\n");
    return 0;
}

int check_rdrand() {
    CPUIDinfo info;
   
    get_cpuid(&info,1,0);
    if (verbose_mode) {
        if ((info.ECX & 0x40000000)==0x40000000) fprintf(stderr,"RdRand bit is set\n");
    }
    if ((info.ECX & 0x40000000)==0x40000000) return 1;
    //fprintf(stderr,"RdRand bit is not set\n");
    return 0;
}

int check_rdseed() {
    CPUIDinfo info;
   
    get_cpuid(&info,7,0);
   
    if (verbose_mode) {
        if ((info.EBX & 0x00040000)==0x00040000) fprintf(stderr,"RdSeed bit is set\n");
    }
    if ((info.EBX & 0x00040000)==0x00040000) return 1;
    //fprintf(stderr,"RdSeed bit is not set\n");
    return 0;
}

int check_aesni()
{
    CPUIDinfo info;
    get_cpuid(&info,1,0);
    
    if ((info.ECX & 0x02000000)==0x02000000) return 1;
    return 0;   
}

int aesni_check_support()
{
    if ((check_is_intel()==1) || (check_is_amd()==1)){
        if (check_aesni()==1) return 1;
    }
    return 0;
}

int rdrand_check_support()
{
    if ((check_is_intel()==1) || (check_is_amd()==1)){
        if (check_rdrand()==1) return 1;
    }
    return 0;
}

int rdseed_check_support()
{
    if ((check_is_intel()==1) || (check_is_amd()==1)){
        if (check_rdseed()==1) return 1;
    }
    return 0;
}

/***************************************************/
/* Gathers 16 bits of entropy through RDRAND       */
/*   The 16 bit result is zero extended to 32 bits */
/*   Writes that entropy to *therand.              */
/*   Returns 1 on success, or 0 on underflow      */
/***************************************************/

int rdrand16_step(unsigned short int *therand)
{
unsigned short int foo;
int cf_error_status;
asm("\n\
        rdrand %%ax;\n\
        mov $1,%%edx;\n\
        cmovae %%ax,%%dx;\n\
        mov %%edx,%1;\n\
        mov %%ax, %0;":"=r"(foo),"=r"(cf_error_status)::"%ax","%dx");
        *therand = foo;
return cf_error_status;

}

int rdseed16_step(unsigned short int *therand)
{
unsigned short int foo;
int cf_error_status;
asm("\n\
        rdseed %%ax;\n\
        mov $1,%%edx;\n\
        cmovae %%ax,%%dx;\n\
        mov %%edx,%1;\n\
        mov %%ax, %0;":"=r"(foo),"=r"(cf_error_status)::"%ax","%dx");
        *therand = foo;
return cf_error_status;
}

/**********************************************/
/* Gathers 32 bits of entropy through RDRAND  */
/*   Writes that entropy to *therand.         */
/*   Returns 1 on success, or 0 on undeerflow */
/**********************************************/

int rdrand32_step(unsigned int *therand)
{
int foo;
int cf_error_status;
asm("\n\
    rdrand %%eax;\n\
        mov $1,%%edx;\n\
        cmovae %%eax,%%edx;\n\
        mov %%edx,%1;\n\
        mov %%eax,%0;":"=r"(foo),"=r"(cf_error_status)::"%eax","%edx");
        *therand = foo;
return cf_error_status;

}

int rdseed32_step(unsigned int *therand)
{
int foo;
int cf_error_status;
asm("\n\
    rdseed %%eax;\n\
        mov $1,%%edx;\n\
        cmovae %%eax,%%edx;\n\
        mov %%edx,%1;\n\
        mov %%eax,%0;":"=r"(foo),"=r"(cf_error_status)::"%eax","%edx");
        *therand = foo;
return cf_error_status;

}

/**********************************************/
/* Gathers 64 bits of entropy through RDRAND  */
/*   Writes that entropy to *therand.         */
/*   Returns 1 on success, or 0 on underflow  */
/**********************************************/

int rdrand64_step(unsigned long long int *therand)
{
unsigned long long int foo;
int cf_error_status;
asm("\n\
        rdrand %%rax;\n\
        mov $1,%%edx;\n\
        cmovae %%rax,%%rdx;\n\
        mov %%edx,%1;\n\
        mov %%rax, %0;":"=r"(foo),"=r"(cf_error_status)::"%rax","%rdx");
        *therand = foo;
return cf_error_status;
}

int rdseed64_step(unsigned long long int *therand)
{
unsigned long long int foo;
int cf_error_status;
asm("\n\
        rdseed %%rax;\n\
        mov $1,%%edx;\n\
        cmovae %%rax,%%rdx;\n\
        mov %%edx,%1;\n\
        mov %%rax, %0;":"=r"(foo),"=r"(cf_error_status)::"%rax","%rdx");
        *therand = foo;
return cf_error_status;
}

/**************************************************/
/* Uses RdRand to acquire a 32 bit random number  */
/*   Writes that entropy to (unsigned int *)dest. */
/*   Will not attempt retry on underflow          */
/*   Returns 1 on success, or 0 on underflow      */
/**************************************************/

int rdrand_get_uint(unsigned int *dest)
{
    unsigned int therand;
    if (rdrand32_step(&therand))
    {
        *dest = therand;
        return 1;
    }
    else return 0;
}

int rdseed_get_uint(unsigned int *dest)
{
    unsigned int therand;
    if (rdseed32_step(&therand))
    {
        *dest = therand;
        return 1;
    }
    else return 0;
}

int rdrand_get_ulint(unsigned long long int *dest)
{
    unsigned long long int therand;
    if (rdrand64_step(&therand))
    {
        *dest = (unsigned long long int)therand;
        return 1;
    }
    else return 0;
}

int rdseed_get_ulint(unsigned long long int *dest)
{
    unsigned long long int therand;
    if (rdseed64_step(&therand))
    {
        *dest = (unsigned long long int)therand;
        return 1;
    }
    else return 0;
}

/**************************************************/
/* Uses RdRand to acquire a 32 bit random number  */
/*   Writes that entropy to (unsigned int *)dest. */
/*   Will retry up to retry_limit times           */
/*   Returns 1 on success, or 0 on underflow      */
/**************************************************/

int rdrand_get_uint_retry(unsigned int retry_limit, unsigned int *dest)
{
int success;
int count;
unsigned int therand;

  count = 0;

  do
  {
    success=rdrand32_step(&therand);
  } while((success == 0) || (count++ < retry_limit));
  
  if (success == 1)
  {
    *dest = therand;
    return 1;
  }
  else
  {
    return 0;
  }
}

int rdseed_get_uint_retry(unsigned int retry_limit, unsigned int *dest)
{
int success;
int count;
unsigned int therand;

  count = 0;

  do
  {
    success=rdseed32_step(&therand);
  } while((success == 0) || (count++ < retry_limit));
  
  if (success == 1)
  {
    *dest = therand;
    return 1;
  }
  else
  {
    return 0;
  }
}

/****************************************************************/
/* Uses RdRand to acquire a block of n 64 bit random numbers    */
/*   Writes that entropy to (unsigned long long int *)dest[0+]. */
/*   Will retry up to retry_limit times                         */
/*   Returns 1 on success, or 0 on underflow                    */
/****************************************************************/

int rdrand_get_n_qints_retry(unsigned int n, unsigned int retry_limit, unsigned long long int *dest)
{
int success=0;
int count=0;
int i=0;

    for (i=0; i<n; i++)
    {
        count = 0;
        do
        {
                success=rdrand64_step(dest);
        } while((success == 0) && (count++ < retry_limit));
        if (success == 0) return 0;
        dest=&(dest[1]);
    }
    return 1; 
}

int rdseed_get_n_qints_retry(unsigned int n, unsigned int retry_limit, unsigned long long int *dest)
{
int success;
int count;
unsigned int i;

    for (i=0; i<n; i++)
    {
        count = 0;
        do
        {
                success=rdseed64_step(dest);
        } while((success == 0) && (count++ < retry_limit));
        if (success == 0) return 0;
        dest=&(dest[1]);
    }
    return 1; 
}


/****************************************************************/
/* Uses RdRand to acquire a block of n 32 bit random numbers    */
/*   Will uses RdRand64 and RdRand32 to optimize speed          */
/*   Writes that entropy to (unsigned int *)dest[0+].           */
/*   Will retry up to retry_limit times                         */
/*   Returns the number of units successfuly acquired           */
/****************************************************************/

int rdrand_get_n_uints_retry(unsigned int n, unsigned int retry_limit, unsigned int *dest)
{
    int qwords;
    int dwords;
    int i;

    unsigned long long int qrand;
    unsigned int drand;

    int success;
    int count;

    int total_uints;

    unsigned long int *qptr;

    total_uints = 0;
    qptr = (unsigned long int*)dest;

    qwords = n/2;
    dwords = n -(qwords*2);

    for (i=0; i<qwords; i++)
    {
        count = 0;
        do
        {
                success=rdrand64_step(&qrand);
        } while((success == 0) || (count++ < retry_limit));

        if (success == 1) 
        {
            *qptr = qrand;
            qptr++;
            total_uints+=2;
        }
        else (i = qwords);
    }
    if ((qwords > 0) && (success == 0)) return total_uints; 

    dest = (unsigned int*)qptr;
        for (i=0; i<dwords; i++)
        {
        count = 0;
                do
                {
                    success=rdrand32_step(&drand);
                } while((success == 0) || (count++ < retry_limit));

                if (success == 1)
                {
                        *dest = qrand;
            dest++;
                        total_uints++;
                }
        else (i = dwords);
        }
        return total_uints;
}

int rdseed_get_n_uints_retry(unsigned int n, unsigned int retry_limit, unsigned int *dest)
{
    int qwords;
    int dwords;
    int i;

    unsigned long long int qrand;
    unsigned int drand;

    int success;
    int count;

    int total_uints;

    unsigned long int *qptr;

    total_uints = 0;
    qptr = (unsigned long int*)dest;

    qwords = n/2;
    dwords = n -(qwords*2);

    for (i=0; i<qwords; i++)
    {
        count = 0;
        do
        {
                success=rdseed64_step(&qrand);
        } while((success == 0) || (count++ < retry_limit));

        if (success == 1) 
        {
            *qptr = qrand;
            qptr++;
            total_uints+=2;
        }
        else (i = qwords);
    }
    if ((qwords > 0) && (success == 0)) return total_uints; 

    dest = (unsigned int*)qptr;
        for (i=0; i<dwords; i++)
        {
        count = 0;
                do
                {
                    success=rdseed32_step(&drand);
                } while((success == 0) || (count++ < retry_limit));

                if (success == 1)
                {
                        *dest = qrand;
            dest++;
                        total_uints++;
                }
        else (i = dwords);
        }
        return total_uints;
}

/****************************************************************/
/* Uses RdRand to acquire a block of n 32 bit random numbers    */
/*   Will uses RdRand64 and RdRand32 to optimize speed          */
/*   Writes that entropy to (unsigned int *)dest[0+].           */
/*   Will not attempt retry on underflow                        */
/*   Returns the number of units successfuly acquired           */
/****************************************************************/

int rdrand_get_n_uints(int n, unsigned int *dest)
{
        int qwords;
        int dwords;
        int i;

        unsigned long long int qrand;
        unsigned int drand;

        int success=0;

        int total_uints;

        unsigned long int *qptr;

        total_uints = 0;
        qptr = (unsigned long int*)dest;

        qwords = n/2;
        dwords = n -(qwords*2);

        for (i=0; i<qwords; i++)
        {
                if (rdrand64_step(&qrand))
                {
                        *qptr = qrand;
                        qptr++;
                        total_uints+=2;
                }
        else (i = qwords);
        }
        if ((qwords > 0) && (success == 0)) return total_uints;

        dest = (unsigned int*)qptr;
        for (i=0; i<dwords; i++)
        {
                if (rdrand32_step(&drand))
                {
                        *dest = drand;
                        dest++;
                        total_uints++;
                }
        else (i = dwords);
        }
        return total_uints;

}

/****************************************************************/
/* Uses RdRand to acquire a block of random bytes               */
/*   Uses RdRand64 to optimize speed                            */
/*   Writes that entropy to (unsigned int *)dest[0+].           */
/*   Internally will retry up to 10 times un underflow.         */
/*   Returns 1 on success, 0 on failure                         */
/****************************************************************/

int rdrand_get_bytes_step(unsigned int n, unsigned char *dest)
{
unsigned char *start;
unsigned char *residualstart;
unsigned long long int *blockstart;
unsigned int count;
unsigned int residual;
unsigned int startlen;
unsigned long long int i;
unsigned long long int temprand;
unsigned int length;

    /* Compute the address of the first 64 bit aligned block in the destination buffer */
    start = dest;
    if (((unsigned long long int)start % (unsigned long long int)8) == 0)
    {
        blockstart = (unsigned long long int *)start;
        count = n;
        startlen = 0;
    }
    else
    {
        blockstart = (unsigned long long int *)(((unsigned long long int)start & ~(unsigned long long int)7)+(unsigned long long int)8);
        count = n - (8 - (unsigned int)((unsigned long long int)start % 8));
        startlen = (unsigned int)((unsigned long long int)blockstart - (unsigned long long int)start);
    }

    /* Compute the number of 64 bit blocks and the remaining number of bytes */
    residual = count % 8;
    length = count >> 3;
    if (residual != 0)
    {
        residualstart = (unsigned char *)(blockstart + length);
    }

    /* Get a temporary random number for use in the residuals. Failout if retry fails */
    if (startlen > 0)
    {
        if (rdrand_get_n_qints_retry(1, 10, (void *)&temprand) == 0) return 0;
    }

    /* populate the starting misaligned block */
    for (i = 0; i<startlen; i++)
    {
        start[i] = (unsigned char)(temprand & 0xff);
        temprand = temprand >> 8;
    }

    /* populate the central aligned block. Fail out if retry fails */
    if (rdrand_get_n_qints_retry(length, 10, (void *)(blockstart)) == 0) return 0;

    /* populate the final misaligned block */
    if (residual > 0)
    {
        if (rdrand_get_n_qints_retry(1, 10, (void *)&temprand) == 0) return 0;
        for (i = 0; i<residual; i++)
        {
            residualstart[i] = (unsigned char)(temprand & 0xff);
            temprand = temprand >> 8;
        }
    }

        return 1;
}






