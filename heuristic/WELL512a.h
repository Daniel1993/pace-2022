#ifndef WELL_512_H_GUARD_
#define WELL_512_H_GUARD_

#include <stdint.h>

void WriteEntropyInBuffer(void *buffer, uint64_t size);
uint32_t WELLRNG512a_u4 (void);

/* ***************************************************************************** */
/* Copyright:      Francois Panneton and Pierre L'Ecuyer, University of Montreal */
/*                 Makoto Matsumoto, Hiroshima University                        */
/* Notice:         This code can be used freely for personal, academic,          */
/*                 or non-commercial purposes. For commercial purposes,          */
/*                 please contact P. L'Ecuyer at: lecuyer@iro.UMontreal.ca       */
/* ***************************************************************************** */

void InitWELLRNG512a (uint32_t *init);
double WELLRNG512a (void);

#endif /* WELL_512_H_GUARD_ */


