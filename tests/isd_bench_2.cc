/* Copyright 2010 ISD-team */
#include <cmath>

#include "./test.h"

void callback(unsigned long *word) {
  int w=0;
  for(int i=0; i< 1024; i++) {
    int b = isd_get_bit(word, i);
    printf("%d", b);
    w += b;
  }
  printf("\n");
  printf("poids %d \n", w);
}

int test(int i) {
  ISD::Decoder CC(0x5eed);

  int n = 1024;
  int k = 525;
  int w = 50;
  int p = 2;
  int l = 20;
  int m = 5;
  int c = 32;
  int r = 8;

  CC.GuessParameters(n, k);

  CC.SetParameter_max_iteration(i);
  CC.SetParameter_max_words(1);
  CC.SetParameter_px(p);
  CC.SetParameter_py(p);

  CC.SetParameter_m(m);
  CC.SetParameter_w(w);
  CC.SetParameter_l(l);
  CC.SetParameter_c(c);
  CC.SetParameter_r(r);

  CC.SetParameter_callback(&callback);

  CC.MatrixRandomize();

  unsigned long long t = rdtsc();

  isd_stats s;
  s.iteration = 0;
  s.n_found = 0;

  int err_code = CC.Run(&s);

  printf("err_code %d\n", err_code);
  printf("iter %d, found %d\n", s.iteration, s.n_found);

  t = rdtsc() - t;
  double clocks = static_cast<double>(t)/i;
  printf("%g clock/iter %f ", clocks, log2(clocks));
  printf("\n");

  return 0;
}

int run(void) {
  return test(2000);
}

int main(int argc, char *argv[]) {
  int res = 0;
  DO(run);
  return res;
}
