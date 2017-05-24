
#include <wmmintrin.h>

__m128i  aes128_keyexpand(__m128i key) 
{
    key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
    key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
    return _mm_xor_si128(key, _mm_slli_si128(key, 4));
}

void  aesni256k128d(unsigned char *key, const unsigned char *plaintext, const unsigned char *ciphertext) 
{
            __m128i rk[15];
            __m128i m;

            /* 256 bit Key Expansion */
            
            rk[0] = _mm_loadu_si128((const __m128i*) key);
            rk[1] = _mm_loadu_si128((const __m128i*) (key+16));

            rk[2] = _mm_xor_si128(aes128_keyexpand(rk[0]), _mm_shuffle_epi32(_mm_aeskeygenassist_si128(rk[1], 0x01), 0xff));
            rk[3] = _mm_xor_si128(aes128_keyexpand(rk[1]), _mm_shuffle_epi32(_mm_aeskeygenassist_si128(rk[2], 0x00), 0xaa));
            
            rk[4] = _mm_xor_si128(aes128_keyexpand(rk[2]), _mm_shuffle_epi32(_mm_aeskeygenassist_si128(rk[3], 0x02), 0xff));
            rk[5] = _mm_xor_si128(aes128_keyexpand(rk[3]), _mm_shuffle_epi32(_mm_aeskeygenassist_si128(rk[4], 0x00), 0xaa));
            
            rk[6] = _mm_xor_si128(aes128_keyexpand(rk[4]), _mm_shuffle_epi32(_mm_aeskeygenassist_si128(rk[5], 0x04), 0xff));
            rk[7] = _mm_xor_si128(aes128_keyexpand(rk[5]), _mm_shuffle_epi32(_mm_aeskeygenassist_si128(rk[6], 0x00), 0xaa));
            
            rk[8] = _mm_xor_si128(aes128_keyexpand(rk[6]), _mm_shuffle_epi32(_mm_aeskeygenassist_si128(rk[7], 0x08), 0xff));
            rk[9] = _mm_xor_si128(aes128_keyexpand(rk[7]), _mm_shuffle_epi32(_mm_aeskeygenassist_si128(rk[8], 0x00), 0xaa));
            
            rk[10] = _mm_xor_si128(aes128_keyexpand(rk[8]), _mm_shuffle_epi32(_mm_aeskeygenassist_si128(rk[9], 0x10), 0xff));
            rk[11] = _mm_xor_si128(aes128_keyexpand(rk[9]), _mm_shuffle_epi32(_mm_aeskeygenassist_si128(rk[10], 0x00), 0xaa));
            
            rk[12] = _mm_xor_si128(aes128_keyexpand(rk[10]), _mm_shuffle_epi32(_mm_aeskeygenassist_si128(rk[11], 0x20), 0xff));
            rk[13] = _mm_xor_si128(aes128_keyexpand(rk[11]), _mm_shuffle_epi32(_mm_aeskeygenassist_si128(rk[12], 0x00), 0xaa));
            
            rk[14] = _mm_xor_si128(aes128_keyexpand(rk[12]), _mm_shuffle_epi32(_mm_aeskeygenassist_si128(rk[13], 0x40), 0xff));

            // Do the encrypt

            m = _mm_loadu_si128((const __m128i*) plaintext);

            /* first 9 rounds */
            m = _mm_xor_si128(m, rk[0]);
            m = _mm_aesenc_si128(m, rk[1]);
            m = _mm_aesenc_si128(m, rk[2]);
            m = _mm_aesenc_si128(m, rk[3]);
            m = _mm_aesenc_si128(m, rk[4]);
            m = _mm_aesenc_si128(m, rk[5]);
            m = _mm_aesenc_si128(m, rk[6]);
            m = _mm_aesenc_si128(m, rk[7]);
            m = _mm_aesenc_si128(m, rk[8]);
            m = _mm_aesenc_si128(m, rk[9]);
            m = _mm_aesenc_si128(m, rk[10]);
            m = _mm_aesenc_si128(m, rk[11]);
            m = _mm_aesenc_si128(m, rk[12]);
            m = _mm_aesenc_si128(m, rk[13]);
            m = _mm_aesenclast_si128(m, rk[14]);
            _mm_storeu_si128((__m128i*) ciphertext, m);
}



