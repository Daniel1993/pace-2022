#ifndef BITMAP_H_GUARD
#define BITMAP_H_GUARD

#include <stdint.h>
#include <string.h>

typedef uint64_t* bitmap;

#define bitmapAlloc(_bitmap, _nbV) ({ \
  assert(8 == sizeof(uint64_t) && "does not work"); \
  (_bitmap) = (bitmap)calloc((((_nbV) + 63)/64) + 1, sizeof(uint64_t)); \
  *(_bitmap) = _nbV; \
  (_bitmap) += 1; \
})

#define bitmapClear(_bitmap) ({ \
  memset((_bitmap), 0, ((*((_bitmap) - 1) + 63)/64) * sizeof(uint64_t)); \
})

#define bitmapDestroy(_bitmap) ({ \
  free((_bitmap)-1); \
})

#define bitmapCheck(_bitmap, _v) ({ \
  int _loc = (_v) >> 6; \
  uint64_t _bit = 1UL << ((_v)&0x3F); \
  (((uint64_t*)(_bitmap))[_loc]&_bit) == _bit; \
})


#define bitmapPut(_bitmap, _v) ({ \
  int _loc = (_v) >> 6; \
  uint64_t _bit = 1UL << ((_v)&0x3F); \
  ((uint64_t*)(_bitmap))[_loc] |= _bit; \
})

#endif /* BITMAP_H_GUARD */