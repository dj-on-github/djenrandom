
#include <stdio.h>

extern void aesni128k128d(unsigned char *key, unsigned char *data, unsigned char *ciphertext);
extern void aessw128k128d(unsigned char *key, unsigned char *data, unsigned char *ciphertext);

int aesni_supported;

void printsample(unsigned char *thesample)
{
    int tempindex;
    int j;
    int i;
    tempindex = 0;
    for (i=0;i<16;i++) printf("%02X",thesample[i]);
    printf("\n");
}


int main() {

    unsigned char key[16];
    unsigned char data[16];
    unsigned char ciphertext[16];

    int i;

    for(i=0;i<16;i++) {
        key[i]=0x00;
        data[i]=0x00;
        ciphertext[i]=0x00;
    }
    
    key[15]=0x01;
    data[15]=0x10;

    aesni128k128d(key,data,ciphertext);
    printsample(ciphertext);

    aessw128k128d(key,data,ciphertext);
    printsample(ciphertext);

}

