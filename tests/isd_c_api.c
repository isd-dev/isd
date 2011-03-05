/* Copyright 2010 ISD-team */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include "test.h"

const char *generating_matrix[17] = {
  "1000000000000010001100110100000010100010001110000011111011111010001111",
  "0100000000000000000100000110100010111100111010001101111110101011111001",
  "0010000000000010001111110101110100000100011011100011110101111100100000",
  "0001000000000010001000001000001111101000010110001010000100010011101001",
  "0000100000000000001101010010100010100001110111110101110110011111000000",
  "0000010000000010000000010011010011001110111111101010111011001001100110",
  "0000001000000010001110110101011111100100110011010000011101010000111011",
  "0000000100000000000110111010110011010010010010110100110100010111100111",
  "0000000010000000001101010010011001101001100111111100001101011001000010",
  "0000000001000010001000011111001010111110111101010001001100101010100100",
  "0000000000100000001011101010101011001011001111010011110110000111010011",
  "0000000000010000001011010011001110110010001100101010001001110001100101",
  "0000000000001010000110101101000101101010001000110101111110000010110010",
  "0000000000000110001011101110110110110111001110010100111111111010110100",
  "0000000000000001000010011100011000000101000101011110000000000110010110",
  "0000000000000000101110110011100110101011001000010100101000111111011001",
  "0000000000000000011111110010100111001111111111111101011100010001101001"
};

const char *min_weight_word =
    "0001001000000001000100100001001000001001100000000100011001000101000100";

const char *send_word =
    "1100011000111110010100011100011101011111000010010001000001010111110010";

const char *error_pattern =
    "0001000000010000000000000000000000000001110001000000000000100000010000";

isd_parameters param;
isd_stats stats;

unsigned long *found;

void callback(unsigned long *w) {
  int byte_size = isd_size_in_limbs(param.n) * sizeof(unsigned long);
  if (memcmp(found, w, byte_size) == 0)
    stats.n_found--; /* it was the same */
  else
    memcpy(found, w, byte_size);
}

int min_weight(void) {
  int i, j;
  unsigned long **matrix;
  int err;
  int res;

  param.n = 70;
  param.k = 17;
  param.w = 17;
  param.px = 2;
  param.py = 2;
  param.l = 10;
  param.m = 1;
  param.c = 1;
  param.r = 1;
  param.kx = 8;
  param.ky = 9;
  param.max_iteration = 5000;
  param.max_words = 2;
  param.seed = 0x5eed;
  param.callback = &callback;

  stats.iteration = 0;
  stats.n_found = 0;

  found = (unsigned long *)calloc(isd_size_in_limbs(param.n), sizeof(unsigned long));

  matrix = (unsigned long **)malloc(param.k * sizeof(unsigned long *));
  for (i = 0; i < param.k; ++i)
    matrix[i] = (unsigned long *)calloc(isd_size_in_limbs(param.n), sizeof(unsigned long));;

  for (i = 0; i < param.k; ++i)
    for (j = 0; j < param.n; ++j) {
      if (generating_matrix[i][j] == '1')
        isd_set_bit(matrix[i], j);
    }

  err = isd_low_weight_word((const unsigned long **)matrix, &param, &stats);

  res = 0;

  if (stats.n_found != 1) {
    printf("expected 1 obtained %d\n", stats.n_found);
    res = -1;
  } else {
    for (i = 0; i < param.n; ++i)
      if (isd_get_bit(found, i) != min_weight_word[i] - '0')
        res = -2;
  }

  return res;
}

int decode(void) {
  int i, j;
  /*  isd_parameters param;*/
  unsigned long **matrix;
  unsigned long *received;
  int err;
  int res;

  param.n = 70;
  param.k = 17;
  param.w = 8;
  param.px = 2;
  param.py = 1;
  param.l = 10;
  param.m = 1;
  param.c = 1;
  param.r = 1;
  param.kx = 10;
  param.ky = 10;
  param.max_iteration = 5000;
  param.max_words = 2;
  param.seed = 0x5eed;
  param.callback = &callback;

  stats.iteration = 0;
  stats.n_found = 0;

  matrix = (unsigned long **)malloc(param.k * sizeof(unsigned long *));
  for (i = 0; i < param.k; ++i)
    matrix[i] = (unsigned long *)calloc(isd_size_in_limbs(param.n), sizeof(unsigned long));;

  for (i = 0; i < param.k; ++i)
    for (j = 0; j < param.n; ++j) {
      if (generating_matrix[i][j] == '1')
        isd_set_bit(matrix[i], j);
    }

  received = (unsigned long *)calloc(isd_size_in_limbs(param.n), sizeof(unsigned long));
  for (i = 0; i < param.n; ++i) {
    if (send_word[i] == '1') isd_set_bit(received, i);
    if (error_pattern[i] == '1')
      isd_swap_bit(received, i);
  }

  err = isd_decode((const unsigned long **)matrix, received, &param, &stats);

  res = 0;

  if (stats.n_found != 1) {
    printf("expected 1 obtained %d\n", stats.n_found);
    res = -1;
  } else {
    for (i = 0; i < param.n; ++i)
      if (isd_get_bit(found, i) != error_pattern[i] - '0') {
        printf("bit %d false\n", i);
        res = -1;
      }
  }

  return res;
}

int main(int argc, char *argv[]) {
  int res = 0;
  res |= DO(min_weight);
  res |= DO(decode);
  return res;
}
