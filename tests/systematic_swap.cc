/* Copyright 2010 ISD-team */
#include <string>
#include <vector>

#include <gsl/gsl_rng.h>

#include "./test.h"
#include "../src/systematic.h"

int involutive(void) {
  ISD::SystematicGeneratorMatrix M(3417, 153);

  gsl_rng *rng = gsl_rng_alloc(gsl_rng_gfsr4);

  M.Randomize(rng);

  std::string S = M.ToString();
  for (int i = 0; i < 20; ++i)
    if (M.bit(i, i*i)) M.SwapColumns(i, i*i);
  for (int i = 19; i >= 0; --i)
    if (M.bit(i, i*i)) M.SwapColumns(i, i*i);

  if (S != M.ToString())
    return 1;
  else
    return 0;
}

int word(void) {
  ISD::SystematicGeneratorMatrix M(7, 3);
  M.SetFromString("(1111)"
                  "(1100)"
                  "(1010)", '0', '1');
  M.SwapColumns(0, 0);

  std::vector<int> v;
  v.push_back(0);
  v.push_back(1);

  std::string S = M.Code(v);

  if (S != "(0,1,0,1,1,0,0,)")
    return 1;
  else
    return 0;
}

int main(int argc, char *argv[]) {
  int res = 0;
  res |= DO(involutive);
  res |= DO(word);
  return res;
}
