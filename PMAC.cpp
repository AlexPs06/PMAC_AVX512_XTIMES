#include <iostream>
#include <wmmintrin.h>
#include <immintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <string.h>
#define ALIGN(n) __attribute__ ((aligned(n)))
#define pipeline 1
#define swap_if_le(b) \
      _mm_shuffle_epi8(b,_mm_set_epi8(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))

#define EXPAND_ASSIST(v1,v2,v3,v4,shuff_const,aes_const)                    \
    v2 = _mm_aeskeygenassist_si128(v4,aes_const);                           \
    v3 = _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(v3),              \
                                         _mm_castsi128_ps(v1), 16));        \
    v1 = _mm_xor_si128(v1,v3);                                              \
    v3 = _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(v3),              \
                                         _mm_castsi128_ps(v1), 140));       \
    v1 = _mm_xor_si128(v1,v3);                                              \
    v2 = _mm_shuffle_epi32(v2,shuff_const);                                 \
    v1 = _mm_xor_si128(v1,v2)

#define PRE_COMP_BLOCKS 31     /* Must be between 0 and 31 */
static inline unsigned ntz(unsigned x) {
		static const unsigned char tz_table[32] =
		{ 0,  1, 28,  2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17,  4, 8,
		 31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18,  6, 11,  5, 10, 9};
		return tz_table[((uint32_t)((x & -x) * 0x077CB531u)) >> 27];
	}

typedef ALIGN(16) __m128i block;

using namespace std;
static void AES_128_Key_Expansion(const unsigned char *userkey, void *key);
static inline void AES_encrypt(__m128i tmp, __m128i *out,__m128i *key, unsigned rounds);
static void imprimiArreglo(int tam, unsigned char *in );
static void shift_left(unsigned char *x);
static void shift_right(unsigned char *x);
void pmac_aes_init(unsigned char *K_1,unsigned tag_len);        
static void PMAC(unsigned char *K_1, unsigned char *M, int size, unsigned char *T);
static inline void AES_encrypt_512(__m512i tmp, __m512i *out,__m512i *key, unsigned rounds);
static void AES_cast_128_to_512_key2(__m128i *key,__m512i *key_512);

block L[PRE_COMP_BLOCKS+1];     /* Precomputed L(i) values, L[0] = L */
block L_inv;                    /* Precomputed L/x value             */


// #define test 128

// int main(){
    
//     ALIGN(64) unsigned char plaintext[test];
//     for (int i = 0; i < test; i++)
//     {
//         plaintext[i]=i;
//     }
    
//     ALIGN(16) unsigned char tag[16 ]={ 0x00,0x01,0x02,0x03,
//                                        0x04,0x05,0x06,0x07,
//                                        0x08,0x09,0x0a,0x0b,
//                                        0x0c,0x0d,0x0e,0x0f};
//     ALIGN(16) unsigned char K_1[16 ]={ 0x00,0x01,0x02,0x03,
//                                        0x04,0x05,0x06,0x07,
//                                        0x08,0x09,0x0a,0x0b,
//                                        0x0c,0x0d,0x0e,0x0f};
//     ALIGN(16) unsigned char N[16 ]={ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};


//     pmac_aes_init(K_1,16);        

//     PMAC(K_1, plaintext, test, tag);

//     printf("\n");
//     imprimiArreglo(16, tag);
//     return 0;
// }



 char infoString[]= "PMAC AVX512 XTIMES";  /* Each AE implementation must have a global one */

#ifndef MAX_ITER
#define MAX_ITER 16384
#endif

int main(int argc, char **argv)
{
	/* Allocate locals */
	ALIGN(64) unsigned char pt[MAX_ITER] = {0};
	ALIGN(16) unsigned char key[16]={ 0x00,0x01,0x02,0x03,
                                       0x04,0x05,0x06,0x07,
                                       0x08,0x09,0x0a,0x0b,
                                       0x0c,0x0d,0x0e,0x0f};

     ALIGN(16) unsigned char tag[16 ]={ 0x00,0x01,0x02,0x03,
                                       0x04,0x05,0x06,0x07,
                                       0x08,0x09,0x0a,0x0b,
                                       0x0c,0x0d,0x0e,0x0f};
	char outbuf[MAX_ITER*15+1024];
	int iter_list[2048]; /* Populate w/ test lengths, -1 terminated */
	char *outp = outbuf;
	int iters, i, j, len;
	double Hz,sec;
	double ipi=0, tmpd;
	clock_t c;
	iter_list[0] = 64;
	iter_list[1] = 128;
	iter_list[2] = 256;
	iter_list[3] = 512;
	iter_list[4] = 1024;
	iter_list[5] = 2048;
	iter_list[6] = 4096;
	iter_list[7] = 8192;
	iter_list[8] = 16384;
	iter_list[9] = -1;
    /* Create file for writing data */
	FILE *fp = NULL;
    char str_time[25];
	time_t tmp_time = time(NULL);
	struct tm *tp = localtime(&tmp_time);
	strftime(str_time, sizeof(str_time), "%F %R", tp);
	if ((argc < 2) || (argc > 3)) {
		printf("Usage: %s MHz [output_filename]\n", argv[0]);
		return 0;
	} else {
		Hz = 1e6 * strtol(argv[1], (char **)NULL, 10);
		if (argc == 3)
			fp = fopen(argv[2], "w");
	}

	
    outp += sprintf(outp, "%s ", infoString);

    outp += sprintf(outp, "- Run %s\n\n",str_time);

	// outp += sprintf(outp, "Context: %d bytes\n");
    
	printf("Starting run...\n");fflush(stdout);


	/*
	 * Get time for key setup
	 */
	iters = (int)(Hz/520);
	do {
	
		c = clock();
		for (j = 0; j < iters; j++) {
            pmac_aes_init(key,16);        
		}
		c = clock() - c;
		sec = c/(double)CLOCKS_PER_SEC;
		tmpd = (sec * Hz) / (iters);
		
		if ((sec < 1.2)||(sec > 1.3))
			iters = (int)(iters * 5.0/(4.0 * sec));
		printf("%f\n", sec);
	} while ((sec < 1.2) || (sec > 1.3));
	
	printf("key -- %.2f (%d cycles)\n",sec,(int)tmpd);fflush(stdout);
	outp += sprintf(outp, "Key setup: %d cycles\n\n", (int)tmpd);

	/*
	 * Get times over different lengths
	 */
	iters = (int)(Hz/1000);
	i=0;
	len = iter_list[0];
	while (len >= 0) {
	
		do {
		

            PMAC(key,pt,iter_list[i],tag);

			c = clock();
			for (j = 0; j < iters; j++) {
                PMAC(key,pt,iter_list[i],tag);
			}
			c = clock() - c;
			sec = c/(double)CLOCKS_PER_SEC;
			tmpd = (sec * Hz) / ((double)len * iters);
			
			if ((sec < 1.2)||(sec > 1.3))
				iters = (int)(iters * 5.0/(4.0 * sec));
			
		} while ((sec < 1.2) || (sec > 1.3));
		
		printf("%d -- %.2f  (%6.2f cpb)\n",len,sec,tmpd);fflush(stdout);
		outp += sprintf(outp, "%5d  %6.2f\n", len, tmpd);
		if (len==44) {
			ipi += 0.05 * tmpd;
		} else if (len==552) {
			ipi += 0.15 * tmpd;
		} else if (len==576) {
			ipi += 0.2 * tmpd;
		} else if (len==1500) {
			ipi += 0.6 * tmpd;
		}
		
		++i;
		len = iter_list[i];
	}
	outp += sprintf(outp, "ipi %.2f\n", ipi);
	if (fp) {
        fprintf(fp, "%s", outbuf);
        fclose(fp);
    } else
        fprintf(stdout, "%s", outbuf);

	return ((pt[0]==12) && (pt[10]==34) && (pt[20]==56) && (pt[30]==78));
}


static void PMAC(unsigned char *K_1, unsigned char *M, int size, unsigned char *T){

    int m_blocks = 0;
    if (size%16==0)
        m_blocks=(size/16);
    else
        m_blocks=(size/16) + 1;
    int size_copy=size;
    static __m512i * plain_text = (__m512i*) M;
    unsigned char tmp_2[16] = {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
    __m128i Tag;
    __m128i Offset_128;
    __m512i checksum;
    __m128i checksum_128;
    __m128i keys_128[11];
    __m512i keys_512[11];
    __m128i sum_nonce= _mm_set_epi32(0,0,0,1);

    AES_128_Key_Expansion(K_1, keys_128);
    AES_cast_128_to_512_key2(keys_128, keys_512);
    checksum = _mm512_setzero_si512();
    checksum_128 = _mm_setzero_si128();
    Tag = _mm_setzero_si128();
    


    static union {__m512i bl512; block bl128[4]; } Offset;
    static union {__m512i bl512; block bl128[4]; } tmp;

    Offset.bl512=_mm512_setzero_si512();
    Offset_128 =_mm_setzero_si128();
    
    int j=1;
    int k=1;

    for (int i = 4; i < m_blocks; i+=4){

        Offset.bl128[0] = _mm_xor_si128( L[0] ,Offset.bl128[3]);        
        Offset.bl128[1] = _mm_xor_si128( L[1] ,Offset.bl128[0]);        
        Offset.bl128[2] = _mm_xor_si128( L[0] ,Offset.bl128[1]);        
        Offset.bl128[3] = _mm_xor_si128( L[ntz(i)] ,Offset.bl128[2]);        

        plain_text[j-1] =_mm512_xor_si512(plain_text[j-1],Offset.bl512);
        
        AES_encrypt_512(plain_text[j-1], &plain_text[j-1], keys_512, 10);

        checksum =_mm512_xor_si512(plain_text[j-1],checksum);
        size_copy-=64;
        k+=4;
        j++;

    }
    tmp.bl512=checksum;
    for (int l = 0; l < 4; l++){
        checksum_128=_mm_xor_si128( tmp.bl128[l] ,checksum_128);
    }


    

    if (size_copy==64){
        Offset.bl128[0] = _mm_xor_si128( L[0] ,Offset.bl128[3]);
        k++;        
        Offset.bl128[1] = _mm_xor_si128( L[1] ,Offset.bl128[0]);        
        k++;        
        Offset.bl128[2] = _mm_xor_si128( L[0] ,Offset.bl128[1]);        
        k++;        
        Offset.bl128[3] = _mm_xor_si128( L[ntz(k)] ,Offset.bl128[2]); 

        tmp.bl512 = plain_text[j-1];

        checksum_128=_mm_xor_si128( tmp.bl128[3] ,checksum_128);

        plain_text[j-1] =_mm512_xor_si512(plain_text[j-1],Offset.bl512);

        AES_encrypt_512(plain_text[j-1], &plain_text[j-1], keys_512, 10);
        
        tmp.bl512=plain_text[j-1];
        for (int l = 0; l < 3; l++){
            checksum_128=_mm_xor_si128( tmp.bl128[l] ,checksum_128);
        }
        
        checksum_128=_mm_xor_si128( L_inv ,checksum_128);

    }else{ 
        tmp.bl512 = plain_text[j-1];
        Offset_128=Offset.bl128[3];
        while(size_copy >16){
            Offset_128 = _mm_xor_si128( L[ntz(k)] ,Offset_128);   

            tmp.bl128[(k%4)-1]=_mm_xor_si128(tmp.bl128[(k%4)-1],Offset_128);

            AES_encrypt(tmp.bl128[(k%4)-1], &tmp.bl128[(k%4)-1], keys_128, 10);
        
            checksum_128 =_mm_xor_si128(tmp.bl128[(k%4)-1],checksum_128);
            size_copy-=16;
            k++;
            
        } 


        if (size_copy==16){
            checksum_128 = _mm_xor_si128( tmp.bl128[(k%4)-1] ,checksum_128);        
            checksum_128 = _mm_xor_si128( L_inv ,checksum_128);      

        }else{
            __m128i tmp_3;
            memcpy(tmp_2, M+((k-1)*16) , size_copy);
            int k=0;
            tmp_2[size_copy]=0x80;

            tmp_3 = _mm_setzero_si128();
            tmp_3 = _mm_load_si128((__m128i *)&tmp_2);
            checksum_128 = _mm_xor_si128( tmp_3 ,checksum_128); 
        }

    }
    
    

    AES_encrypt(checksum_128, &Tag, keys_128, 10);
	_mm_store_si128 ((__m128i*)T,Tag);
}


static inline void AES_encrypt(__m128i tmp, __m128i *out,__m128i *key, unsigned rounds){
	int j;
	tmp = _mm_xor_si128 (tmp,key[0]);
	for (j=1; j<rounds; j++)  tmp = _mm_aesenc_si128 (tmp,key[j]);
	tmp = _mm_aesenclast_si128 (tmp,key[j]);
	_mm_store_si128 ((__m128i*)out,tmp);
}

static inline void AES_encrypt_512(__m512i tmp, __m512i *out,__m512i *key, unsigned rounds){
	int j;
	tmp = _mm512_xor_si512 (tmp,key[0]);
	for (j=1; j<rounds; j++)  tmp = _mm512_aesenc_epi128 (tmp,key[j]);
	tmp = _mm512_aesenclast_epi128 (tmp,key[j]);
	_mm512_store_si512 ((__m128i*)out,tmp);
}

static void AES_cast_128_to_512_key2(__m128i *key,__m512i *key_512){
    union {__m128i oa128[4]; __m512i oa512;} oa;
    for(int i = 0; i< 11; i++ ){
        oa.oa128[0] = key[i];
        oa.oa128[1] = key[i];
        oa.oa128[2] = key[i];
        oa.oa128[3] = key[i];
        key_512[i]=oa.oa512;
    }

}

static void AES_128_Key_Expansion(const unsigned char *userkey, void *key)
{
    __m128i x0,x1,x2;
    __m128i *kp = (__m128i *)key;
    kp[0] = x0 = _mm_loadu_si128((__m128i*)userkey);
    x2 = _mm_setzero_si128();
    EXPAND_ASSIST(x0,x1,x2,x0,255,1);   kp[1]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,2);   kp[2]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,4);   kp[3]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,8);   kp[4]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,16);  kp[5]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,32);  kp[6]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,64);  kp[7]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,128); kp[8]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,27);  kp[9]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,54);  kp[10] = x0;
}



void imprimiArreglo(int tam, unsigned char *in )
{

    for (int i = 0; i<tam; i++){
        printf("%02x", in[i] );
    }
    printf("\n" );

}


/************************************************************************* 
 * ocb_aes_init 
 *************************************************************************/
// keystruct *                         /* Init'd keystruct or NULL      */
void pmac_aes_init(unsigned char   *K_1,    /* AES key                       */
             unsigned   tag_len    /* Length of tags to be used     */
             )        /* OCB key structure. NULL means */
                                    /* Allocate/init new, non-NULL   */
                                    /* means init existing structure */
{
    ALIGN(16) unsigned char tmp[16] = {0,};
    ALIGN(16) __m128i tmp_2 = _mm_setzero_si128();
    ALIGN(16) __m128i enc_key[11];
    
    unsigned first_bit, last_bit, i;

        /* Initialize AES keys.   (Note that if one is only going to 
           encrypt, key->rdk can be eliminated */
        AES_128_Key_Expansion(K_1, enc_key);

        AES_encrypt(  tmp_2, (__m128i*)&tmp[0], enc_key, 10);
        /* Precompute L[i]-values. L[0] is synonym of L */

        for (i = 0; i <= PRE_COMP_BLOCKS; i++) {

            memcpy(&L[i], tmp, 16);   /* Copy tmp to L[i] */
            first_bit = tmp[0] & 0x80u;    /* and multiply tmp by x */
            shift_left(tmp);
            if (first_bit) 
                tmp[15] ^= 0x87;
            
        }

        /* Precompute L_inv = L . x^{-1} */
        memcpy(tmp, (unsigned char*)&L[0], 16);
        last_bit = tmp[15] & 0x01;
        shift_right(tmp);

        if (last_bit) {
            tmp[0] ^= 0x80;
            tmp[15] ^= 0x43;
        }
        memcpy(&L_inv, tmp, 16);

        /* Set tag length used for this session */
}
    


/************************************************************************* 
 * shift_left  
 *************************************************************************/
static void
shift_left(unsigned char *x)
/* 128-bit shift-left by 1 bit: *x <<= 1.                                */
{
    int i;
    for (i = 0; i < 15; i++) {
        x[i] = (x[i] << 1) | (x[i+1] & 0x80 ? 1 : 0);
    }
    x[15] = (x[15] << 1);
}

/************************************************************************* 
 * shift_right 
 *************************************************************************/
static void
shift_right(unsigned char *x)
/* 128-bit shift-right by 1 bit:  *x >>= 1                               */
{
    int i;
    for (i = 15; i > 0; i--) {
        x[i] = (x[i] >> 1) | (x[i-1] & 1 ? 0x80u : 0);
    }
    x[0] = (x[0] >> 1);
}



