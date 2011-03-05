/* Copyright 2010 ISD-team */
#include <cmath>
#include <string>
#include <vector>

#include "./test.h"

int test(int i) {
  ISD::Decoder CC(0x5eed);

  int n = 512;
  int k = 256;
  int w = 57;
  int p = 2;
  int l = 16;
  int m = 4;
  int c = 6;
  int r = 1;


  CC.GuessParameters(n, k);
  CC.SetParameter_w(w);
  CC.SetParameter_px(p);
  CC.SetParameter_py(p);
  CC.SetParameter_l(l);
  CC.SetParameter_m(m);
  CC.SetParameter_c(c);
  CC.SetParameter_r(r);
  CC.SetParameter_max_iteration(i);
  CC.SetParameter_max_words(100);

  CC.MatrixRandomize();

  unsigned long long t = rdtsc();
  CC.Run(NULL);
  double clocks = static_cast<double>(rdtsc() - t)/i;
  printf("%g clock/iter\n", clocks);

  return 0;
}

int run(void) {
  return test(10000);
}

int main(int argc, char *argv[]) {
  int res = 0;
  DO(run);
  return res;
}
