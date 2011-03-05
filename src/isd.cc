/*
  Copyright 2010 ISD-team
*/

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <string>
#include <algorithm>
#include <utility>

#include "./utils.h"
#include "./isd.h"
#include "./isd-impl.h"


namespace ISD {

// ****************************************************************************
// Interface, forward to implementation ***************************************
// ****************************************************************************

Decoder::Decoder(unsigned int seed) : impl_(new Impl(seed)) { }

Decoder::~Decoder(void) {
  delete impl_;
}

size_t Decoder::Run(isd_stats *stats) {
  return impl_->Run(stats);
}

void Decoder::SetSizes(int n, int k) {
  impl_->SetSizes(n, k);
}

void Decoder::SetMatrix(int n, int k, const unsigned long **matrix) {
  impl_->SetMatrix(n, k, matrix);
}

void Decoder::SetMatrixFromString
(int n, int k, const std::string &s) {
  impl_->SetMatrixFromString(n, k, s);
}

void Decoder::SetMatrixFromFile
(int n, int k, const FILE *f) {
  impl_->SetMatrixFromFile(n, k, f);
}

void Decoder::SetSystematicMatrix(int n, int k, const unsigned long **matrix) {
  impl_->SetSystematicMatrix(n, k, matrix);
}

void Decoder::SetSystematicMatrixFromString
(int n, int k, const std::string &s) {
  impl_->SetSystematicMatrixFromString(n, k, s);
}

void Decoder::SetSystematicMatrixFromFile
(int n, int k, const FILE *f) {
  impl_->SetSystematicMatrixFromFile(n, k, f);
}

bool Decoder::AppendRow(const unsigned long *row) {
  return impl_->AppendRow(row);
}

void Decoder::MatrixRandomize(void) {
  impl_->MatrixRandomize();
}

void Decoder::SetParameter_w(int w) {
  impl_->SetParameter_w(w);
}

void Decoder::SetParameter_px(int p) {
  impl_->SetParameter_px(p);
}

void Decoder::SetParameter_py(int p) {
  impl_->SetParameter_py(p);
}

void Decoder::SetParameter_m(int m) {
  impl_->SetParameter_m(m);
}

void Decoder::SetParameter_l(int l) {
  impl_->SetParameter_l(l);
}

void Decoder::SetParameter_c(int c) {
  impl_->SetParameter_c(c);
}

void Decoder::SetParameter_r(int r) {
  impl_->SetParameter_r(r);
}

void Decoder::SetParameter_kx(int kx) {
  impl_->SetParameter_kx(kx);
}

void Decoder::SetParameter_ky(int ky) {
  impl_->SetParameter_ky(ky);
}

void Decoder::SetParameter_max_iteration(int maxiter) {
  impl_->SetParameter_max_iteration(maxiter);
}

void Decoder::SetParameter_max_words(int maxiter) {
  impl_->SetParameter_max_words(maxiter);
}

void Decoder::SetParameter_callback(void (*callback)(unsigned long *)) {
  impl_->SetParameter_callback(callback);
}

void Decoder::SetParameters(const isd_parameters &param) {
  impl_->SetParameters(param);
}

void Decoder::GuessParameters(int n, int k) {
  impl_->GuessParameters(n, k);
}

const isd_parameters *Decoder::parameters(void) {
  return impl_->parameters();
}

}  // namespace ISD

/////////////////
// C interface //
/////////////////

int isd_low_weight_word(const unsigned long **matrix,
                        const isd_parameters *param, isd_stats *stats) {
  ISD::Decoder D(param->seed);
  D.SetMatrix(param->n, param->k, matrix);
  D.SetParameters(*param);

  return D.Run(stats);
}

int isd_low_weight_word_syst(const unsigned long **matrix,
                             const isd_parameters *param, isd_stats *stats) {
  ISD::Decoder D(param->seed);
  D.SetSystematicMatrix(param->n, param->k, matrix);
  D.SetParameters(*param);

  return D.Run(stats);
}

/*
int do_decode(ISD::Decoder *D, const isd_parameters *param,
              const unsigned long *word, unsigned long **res) {
  if (D->AppendRow(word)) {
    return do_low_weight(D, param, res);
  } else {
    // word is in the code
    // only one error pattern possible, all zero
    const int limbs_per_word = ISD::IND(param->n + ISD::ulong_bits - 1);
    (*res) = new unsigned long[limbs_per_word];
    for (int i = 0; i < limbs_per_word; ++i)
      (*res)[i] = 0;
    return 1;
  }
}
*/

int isd_decode(const unsigned long **matrix, const unsigned long *word,
               const isd_parameters *param, isd_stats *stats ) {
  isd_parameters p;
  memcpy(&p, param, sizeof(p));
  p.k += 1;

  ISD::Decoder D(p.seed);
  {
    unsigned long **matrix_augmented;
    matrix_augmented = new unsigned long *[param->k + 1];
    for (int i = 0; i < param->k; ++i)
      matrix_augmented[i] = const_cast<unsigned long *>(matrix[i]);
    matrix_augmented[param->k] = const_cast<unsigned long *>(word);

    D.SetMatrix(p.n, p.k, const_cast<const unsigned long **>(matrix_augmented));
    delete[] matrix_augmented;
  }

  D.SetParameters(p);

  return D.Run(stats);
}

// TODO
/*
int isd_decode_syst(const unsigned long **matrix, const unsigned long *word,
                    const isd_parameters *param, isd_stats *stats);
*/

size_t isd_size_in_limbs(size_t nbits) {
  return ISD::size_in_limbs(nbits);
}

void isd_set_bit(unsigned long *word, int i) {
  word[ISD::IND(i)] |= 1UL << ISD::OFF(i);
}

void isd_unset_bit(unsigned long *word, int i) {
  word[ISD::IND(i)] &= ~(1UL << ISD::OFF(i));
}

void isd_swap_bit(unsigned long *word, int i) {
  word[ISD::IND(i)] ^= 1UL << ISD::OFF(i);
}

int isd_get_bit(const unsigned long *word, int i) {
  return ( word[ISD::IND(i)] >> ISD::OFF(i) ) & 1;
}

void isd_fprintf_parameters(FILE *f, const isd_parameters *param) {
  fprintf(f, "length: %d\n", param->n);
  fprintf(f, "dimension: %d\n", param->k);
  fprintf(f, "weight threshold: %d\n", param->w);
  fprintf(f, "seed: %d\n", param->seed);
  fprintf(f, "px: %d\n", param->px);
  fprintf(f, "py: %d\n", param->py);
  fprintf(f, "m: %d\n", param->m);
  fprintf(f, "l: %d\n", param->l);
  fprintf(f, "c: %d\n", param->c);
  fprintf(f, "r: %d\n", param->r);
  fprintf(f, "maximum iteration: %d\n", param->max_iteration);
  fprintf(f, "maximum number of solutions: %d\n", param->max_words);
  fprintf(f, "callback defined: %d\n", param->callback != NULL);
}
