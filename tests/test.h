#ifndef __ISD_TEST_H__
#define __ISD_TEST_H__

#ifdef __cplusplus
# include <cstdio>
#else
# include <stdio.h>
#endif  /* __cplusplus */

#include "../src/isd.h"

#ifdef ISD_BENCH
#ifdef __i386
extern __inline__ unsigned long long rdtsc() {
  unsigned long long x;
  __asm__ volatile ("rdtsc" : "=A" (x));
  return x;
}
#elif __amd64
extern __inline__ unsigned long long rdtsc() {
  unsigned long long a, d;
  __asm__ volatile ("rdtsc" : "=a" (a), "=d" (d));
  return (d<<32) | a;
}
#else
# define rdtsc __rdtsc
#endif
#endif  /* ISD_BENCH */

#define DO(test) perform(&test, #test)

int perform(int (* test)(void), const char *name) {
  static int compteur = 0;
  int res = test();
  /* printf("do %s\n", name); */
  if (res)
    printf("[ERROR: test %d - %s]\n", compteur, name);
  ++compteur;
  return res;
}

#endif /* __ISD_TEST_H__ */
