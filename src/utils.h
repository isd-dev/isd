#ifndef ISD_UTILS_H_
#define ISD_UTILS_H_

#include <cstdlib>

namespace ISD {

  static const int ulong_size = sizeof(unsigned long);
  static const int ulong_bits = ulong_size * 8;

  inline unsigned int size_in_limbs(size_t i) {
    return (i + ulong_bits - 1) / ulong_bits;
  }

  inline unsigned int OFF(size_t i) {
    return i & (ulong_bits - 1);
  }

  inline unsigned int IND(size_t i) {
    return i / ulong_bits;
  }

  inline unsigned int popcount(unsigned long x) {
#if 0 && (__GNUC__ > 3) || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
    return __builtin_popcountl(x);
#else
# ifndef SIZEOF_UNSIGNED_LONG
#  error "You must define SIZEOF_UNSIGNED_LONG"
# endif
# if SIZEOF_UNSIGNED_LONG == 8
      x -=  (x>>1) & 0x5555555555555555UL;
      x  = ((x>>2) & 0x3333333333333333UL) + (x & 0x3333333333333333UL);
      x  = ((x>>4) + x) & 0x0f0f0f0f0f0f0f0fUL;
      x *= 0x0101010101010101UL;
      return  x>>56;
# elif SIZEOF_UNSIGNED_LONG == 4
      x = (0x55555555UL & x) + (0x55555555UL & (x>> 1));
      x = (0x33333333UL & x) + (0x33333333UL & (x>> 2));
      x = (0x0f0f0f0fUL & x) + (0x0f0f0f0fUL & (x>> 4));
      x = (0x00ff00ffUL & x) + (0x00ff00ffUL & (x>> 8));
      x = (0x0000ffffUL & x) + (0x0000ffffUL & (x>>16));
      return x;
# else
#  error "Unsupported sizeof(unsigned long)"
# endif
#endif
  }

  inline unsigned int count_zeros(unsigned long x) {
#if __GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
    return __builtin_ctzl(x);
#else
    unsigned int r = 0;
    x &= -x;  // isolate lowest bit
    if (ulong_bits == 32) {
      if ( x & 0xffff0000UL )  r += 16;
      if ( x & 0xff00ff00UL )  r += 8;
      if ( x & 0xf0f0f0f0UL )  r += 4;
      if ( x & 0xccccccccUL )  r += 2;
      if ( x & 0xaaaaaaaaUL )  r += 1;
    }
    if (ulong_bits == 64) {
      if ( x & 0xffffffff00000000UL )  r += 32;
      if ( x & 0xffff0000ffff0000UL )  r += 16;
      if ( x & 0xff00ff00ff00ff00UL )  r += 8;
      if ( x & 0xf0f0f0f0f0f0f0f0UL )  r += 4;

      if ( x & 0xccccccccccccccccUL )  r += 2;
      if ( x & 0xaaaaaaaaaaaaaaaaUL )  r += 1;
    }
    return r;
#endif
  }

}

#endif  // ISD_UTILS_H_
