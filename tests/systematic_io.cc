/* Copyright 2010 ISD-team */
#include <string>

#include "./test.h"
#include "../src/systematic.h"

int test0(void) {
  std::string T("000000000"\
                "101010101"\
                "110010011"\
                "111000111"\
                "111101111"\
                "111111111");
  std::string expected = std::string("(\n"\
                                     "(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,),\n"\
                                     "(0,1,0,0,0,0,1,0,1,0,1,0,1,0,1,),\n"\
                                     "(0,0,1,0,0,0,1,1,0,0,1,0,0,1,1,),\n"\
                                     "(0,0,0,1,0,0,1,1,1,0,0,0,1,1,1,),\n"\
                                     "(0,0,0,0,1,0,1,1,1,1,0,1,1,1,1,),\n"\
                                     "(0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,),\n)");
  ISD::SystematicGeneratorMatrix M(15, 6);
  M.SetFromString(T, '0', '1');
  std::string S = M.ToString();
  if (S != expected) {
    printf("expected:\n%s\n obtained:\n%s\n", expected.c_str(), S.c_str());
    return 1;
  }
  return 0;
}

int main(int argc, char *argv[]) {
  int res = 0;
  res |= DO(test0);
  return res;
}
