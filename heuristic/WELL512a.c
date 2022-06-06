#include "WELL512a.h"
// #include <sys/random.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

void WriteEntropyInBuffer(void *buffer, uint64_t size)
{
  char *p = (char *)buffer;
  p[0] =     (char) 123;
  p[1] =     (char) 220;
  p[2] =     (char) 127;
  p[3] =     (char) 216;
  p[4] =     (char) 103;
  p[5] =     (char) 107;
  p[6] =     (char) 129;
  p[7] =     (char)  11;
  p[8] =     (char)  75;
  p[9] =     (char)  86;
  p[10] =    (char)  37;
  p[11] =    (char)  40;
  p[12] =    (char) 156;
  p[13] =    (char) 115;
  p[14] =    (char)  95;
  p[15] =    (char)   3;
  p[16] =    (char) 124;
  p[17] =    (char) 230;
  p[18] =    (char) 134;
  p[19] =    (char)  13;
  p[20] =    (char)  69;
  p[21] =    (char) 162;
  p[22] =    (char)  28;
  p[23] =    (char) 197;
  p[24] =    (char) 120;
  p[25] =    (char)  12;
  p[26] =    (char) 117;
  p[27] =    (char)  67;
  p[28] =    (char) 250;
  p[29] =    (char) 102;
  p[30] =    (char) 229;
  p[31] =    (char)  73;
  p[32] =    (char)  49;
  p[33] =    (char) 132;
  p[34] =    (char) 255;
  p[35] =    (char)  98;
  p[36] =    (char) 111;
  p[37] =    (char)  27;
  p[38] =    (char) 253;
  p[39] =    (char)  77;
  p[40] =    (char) 119;
  p[41] =    (char) 142;
  p[42] =    (char) 187;
  p[43] =    (char) 131;
  p[44] =    (char)   8;
  p[45] =    (char) 140;
  p[46] =    (char) 169;
  p[47] =    (char)  85;
  p[48] =    (char) 151;
  p[49] =    (char)  19;
  p[50] =    (char)  52;
  p[51] =    (char) 247;
  p[52] =    (char) 194;
  p[53] =    (char)  82;
  p[54] =    (char) 183;
  p[55] =    (char) 138;
  p[56] =    (char) 163;
  p[57] =    (char)   7;
  p[58] =    (char)  56;
  p[59] =    (char)   5;
  p[60] =    (char)  93;
  p[61] =    (char)  45;
  p[62] =    (char) 191;
  p[63] =    (char)  82;
}

/* ***************************************************************************** */
/* Copyright:      Francois Panneton and Pierre L'Ecuyer, University of Montreal */
/*                 Makoto Matsumoto, Hiroshima University                        */
/* Notice:         This code can be used freely for personal, academic,          */
/*                 or non-commercial purposes. For commercial purposes,          */
/*                 please contact P. L'Ecuyer at: lecuyer@iro.UMontreal.ca       */
/* ***************************************************************************** */

#define W 32
#define R 16
#define P 0
#define M1 13
#define M2 9
#define M3 5

#define MAT0POS(t,v) (v^(v>>t))
#define MAT0NEG(t,v) (v^(v<<(-(t))))
#define MAT3NEG(t,v) (v<<(-(t)))
#define MAT4NEG(t,b,v) (v ^ ((v<<(-(t))) & b))

#define V0            STATE[state_i                   ]
#define VM1           STATE[(state_i+M1) & 0x0000000fU]
#define VM2           STATE[(state_i+M2) & 0x0000000fU]
#define VM3           STATE[(state_i+M3) & 0x0000000fU]
#define VRm1          STATE[(state_i+15) & 0x0000000fU]
#define VRm2          STATE[(state_i+14) & 0x0000000fU]
#define newV0         STATE[(state_i+15) & 0x0000000fU]
#define newV1         STATE[state_i                 ]
#define newVRm1       STATE[(state_i+14) & 0x0000000fU]

#define FACT 2.32830643653869628906e-10

static uint32_t state_i = 0;
static uint32_t STATE[R];
static uint32_t z0, z1, z2;

void InitWELLRNG512a (uint32_t *init){
   int j;
   state_i = 0;
   for (j = 0; j < R; j++)
     STATE[j] = init[j];
}

double WELLRNG512a (void){

  z0    = VRm1;
  z1    = MAT0NEG (-16,V0)    ^ MAT0NEG (-15, VM1);
  z2    = MAT0POS (11, VM2)  ;
  newV1 = z1                  ^ z2; 
  newV0 = MAT0NEG (-2,z0)     ^ MAT0NEG(-18,z1)    ^ MAT3NEG(-28,z2) ^ MAT4NEG(-5,0xda442d24U,newV1) ;
  state_i = (state_i + 15) & 0x0000000fU;
  return ((double) STATE[state_i]) * FACT;
}

/* ***************************************************************************** */

uint32_t WELLRNG512a_u4 (void)
{

  z0    = VRm1;
  z1    = MAT0NEG (-16,V0)    ^ MAT0NEG (-15, VM1);
  z2    = MAT0POS (11, VM2)  ;
  newV1 = z1                  ^ z2; 
  newV0 = MAT0NEG (-2,z0)     ^ MAT0NEG(-18,z1)    ^ MAT3NEG(-28,z2) ^ MAT4NEG(-5,0xda442d24U,newV1) ;
  state_i = (state_i + 15) & 0x0000000fU;
  return STATE[state_i];
}
