/****************************************************************/
/* 802.11 AES Block Cipher                                      */
/* Implements AES128. I.E. 128 bit key with 128 bit data.       */
/*                                                              */
/* Originally (c) 2002, David Johnston                          */
/* This code has been placed in the public domain by the author */
/* Author: David Johnston                                       */
/* Email (home) : dj@deadhat.com                                */
/* Email (general) : david.johnston@ieee.org                    */
/* Email (work) : dj.johnston@intel.com                         */
/****************************************************************/
 
#include <stdlib.h>
#include <stdio.h>

#define DEBUG 0
/*****************************/
/******** SBOX Table *********/
/*****************************/

    unsigned char sbox_table[256] =
    {
        0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5,
        0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
        0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0,
        0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
        0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc,
        0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
        0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a,
        0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
        0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0,
        0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
        0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b,
        0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
        0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85,
        0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
        0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5,
        0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
        0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17,
        0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
        0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88,
        0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
        0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c,
        0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
        0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9,
        0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
        0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6,
        0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
        0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e,
        0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
        0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94,
        0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
        0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68,
        0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16
    };

/*****************************/
/**** Function Prototypes ****/
/*****************************/

void xor_128(unsigned char *a, unsigned char *b, unsigned char *out);
void xor_32(unsigned char *a, unsigned char *b, unsigned char *out);
unsigned char sbox(unsigned char a);
void next_key(unsigned char *key, int round);
void byte_sub(unsigned char *in, unsigned char *out);
void shift_row(unsigned char *in, unsigned char *out);
void mix_column(unsigned char *in, unsigned char *out);
void add_round_key( unsigned char *shiftrow_in,
                    unsigned char *mcol_in,
                    unsigned char *block_in,
                    int round,
                    unsigned char *out);
void aes128k128d(unsigned char *key, unsigned char *data, unsigned char *ciphertext);

/****************************************/
/* aes128k128d()                        */
/* Performs a 128 bit AES encrypt with  */
/* 128 bit data.                        */
/****************************************/
void xor_128(unsigned char *a, unsigned char *b, unsigned char *out)
{
    int i;
    for (i=0;i<16; i++)
    {
        out[i] = a[i] ^ b[i];
    }
}

void xor_32(unsigned char *a, unsigned char *b, unsigned char *out)
{
    int i;
    for (i=0;i<4; i++)
    {
        out[i] = a[i] ^ b[i];
    }
}
unsigned char sbox(unsigned char a)
{
    return sbox_table[(int)a];
}

void next_key(unsigned char *key, int round)
{
    unsigned char rcon;
    unsigned char sbox_key[4];
    unsigned char rcon_table[12] =
    {
        0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80,
        0x1b, 0x36, 0x36, 0x36
    };

    sbox_key[0] = sbox(key[13]);
    sbox_key[1] = sbox(key[14]);
    sbox_key[2] = sbox(key[15]);
    sbox_key[3] = sbox(key[12]);

    rcon = rcon_table[round];
if (DEBUG >1) printf("rcon : %02X\n",rcon);
    xor_32(&key[0], sbox_key, &key[0]);
    key[0] = key[0] ^ rcon;

    xor_32(&key[4], &key[0], &key[4]);
    xor_32(&key[8], &key[4], &key[8]);
    xor_32(&key[12], &key[8], &key[12]);
}

void byte_sub(unsigned char *in, unsigned char *out)
{
    int i;
    for (i=0; i< 16; i++)
    {
        out[i] = sbox(in[i]);
    }
}

void shift_row(unsigned char *in, unsigned char *out)
{
    out[0] =  in[0];
    out[1] =  in[5];
    out[2] =  in[10];
    out[3] =  in[15];
    out[4] =  in[4];
    out[5] =  in[9];
    out[6] =  in[14];
    out[7] =  in[3];
    out[8] =  in[8];
    out[9] =  in[13];
    out[10] = in[2];
    out[11] = in[7];
    out[12] = in[12];
    out[13] = in[1];
    out[14] = in[6];
    out[15] = in[11];
}
void print32(char *astr, unsigned char *bits)
{
unsigned int data;

	data = 0;
	data = (unsigned int)bits[3];
	data = (data << 8) + (unsigned int)bits[2];
	data = (data << 8) + (unsigned int)bits[1];
	data = (data << 8) + (unsigned int)bits[0];
        if (DEBUG > 0)
        {
                printf("%s : %08X\n",astr,data);
        }
}

void mix_column(unsigned char *in, unsigned char *out)
{
    int i;
    unsigned char add1b[4];
    unsigned char add1bf7[4];
    unsigned char rotl[4];
    unsigned char swap_halfs[4];
    unsigned char andf7[4];
    unsigned char rotr[4];
    unsigned char temp[4];
    unsigned char tempb[4];

 if (DEBUG > 1) print32("in",in);
    for (i=0 ; i<4; i++)
    {
        if ((in[i] & 0x80)== 0x80)
            add1b[i] = 0x1b;
        else
            add1b[i] = 0x00;
    }
 if (DEBUG > 1) print32("add1b",add1b);

    swap_halfs[0] = in[2];    /* Swap halfs */
    swap_halfs[1] = in[3];
    swap_halfs[2] = in[0];
    swap_halfs[3] = in[1];
 if (DEBUG > 1) print32("swap_halfs",swap_halfs);

    rotl[0] = in[3];        /* Rotate left 8 bits */
    rotl[1] = in[0];
    rotl[2] = in[1];
    rotl[3] = in[2];
 if (DEBUG > 1) print32("rotl",rotl);

    andf7[0] = in[0] & 0x7f;
    andf7[1] = in[1] & 0x7f;
    andf7[2] = in[2] & 0x7f;
    andf7[3] = in[3] & 0x7f;
if (DEBUG > 1) print32("andf7",andf7);

    for (i = 3; i>0; i--)    /* logical shift left 1 bit */
    {
        andf7[i] = andf7[i] << 1;
        if ((andf7[i-1] & 0x80) == 0x80)
        {
            andf7[i] = (andf7[i] | 0x01);
        }
    }
    andf7[0] = andf7[0] << 1;
    andf7[0] = andf7[0] & 0xfe;
  if (DEBUG > 1) print32("lsl",andf7);

    xor_32(add1b, andf7, add1bf7);
   if (DEBUG > 1) print32("add1bf7", add1bf7);

    xor_32(in, add1bf7, rotr);
 if (DEBUG > 1) print32("pre-rotr", rotr);

    temp[0] = rotr[0];         /* Rotate right 8 bits */
    rotr[0] = rotr[1];
    rotr[1] = rotr[2];
    rotr[2] = rotr[3];
    rotr[3] = temp[0];
  if (DEBUG > 1) print32("rotr", rotr);

    xor_32(add1bf7, rotr, temp);
    xor_32(swap_halfs, rotl,tempb);
    xor_32(temp, tempb, out);
}

void print_round(int round)
{
	if (DEBUG > 0) printf("ROUND %02X\n",round);
}
void print_rk(unsigned char *data)
{
int i;
	if (DEBUG > 0) printf("         RK:");
	if (DEBUG > 0) for(i=0;i<16;i++) printf("%02X",data[i]);
	if (DEBUG > 0) printf("\n");
}

void print_ct(unsigned char *data)
{
int i;
	if (DEBUG > 0) printf("         CT:");
	if (DEBUG > 0) for(i=0;i<16;i++) printf("%02X",data[i]);
	if (DEBUG > 0) printf("\n");
}
void print_bs(unsigned char *data)
{
int i;
	if (DEBUG > 0) printf("         BS:");
	if (DEBUG > 0) for(i=0;i<16;i++) printf("%02X",data[i]);
	if (DEBUG > 0) printf("\n");
}
void print_sr(unsigned char *data)
{
int i;
	if (DEBUG > 0) printf("         SR:");
	if (DEBUG > 0) for(i=0;i<16;i++) printf("%02X",data[i]);
	if (DEBUG > 0) printf("\n");
}
void print_mc(unsigned char *data)
{
int i;
	if (DEBUG > 0) printf("         MC:");
	if (DEBUG > 0) for(i=0;i<16;i++) printf("%02X",data[i]);
	if (DEBUG > 0) printf("\n");
}
void print_ark(unsigned char *data)
{
int i;
	if (DEBUG > 0) printf("        ARK:");
	if (DEBUG > 0) for(i=0;i<16;i++) printf("%02X",data[i]);
	if (DEBUG > 0) printf("\n");
}

void aes128k128d(unsigned char *key, unsigned char *data, unsigned char *ciphertext)
{
    int round;
    int i;
    unsigned char intermediatea[16];
    unsigned char intermediateb[16];
    unsigned char round_key[16];

    for(i=0; i<16; i++) round_key[i] = key[i];

    for (round = 0; round < 11; round++)
    {
        if (round == 0)
        {
	    print_round(round);
            xor_128(round_key, data, ciphertext);
	    print_ark(ciphertext);
            next_key(round_key, round);
	    print_rk(round_key);         
        }
        else if (round == 10)
        {
	    print_round(round);
            byte_sub(ciphertext, intermediatea);
	    print_bs(intermediatea);
            shift_row(intermediatea, intermediateb);
	    print_sr(intermediateb);
            xor_128(intermediateb, round_key, ciphertext);          
	    print_ark(ciphertext);
        }
        else    /* 1 - 9 */
        {
	    print_round(round);
            byte_sub(ciphertext, intermediatea);
	    print_bs(intermediatea);
            shift_row(intermediatea, intermediateb);
	    print_sr(intermediateb);
            mix_column(&intermediateb[0], &intermediatea[0]);
            mix_column(&intermediateb[4], &intermediatea[4]);
            mix_column(&intermediateb[8], &intermediatea[8]);
            mix_column(&intermediateb[12], &intermediatea[12]);
	    print_mc(intermediatea);
            xor_128(intermediatea, round_key, ciphertext);
	    print_ark(ciphertext);
            next_key(round_key, round);
	    print_rk(round_key);         
        }

    }

}


