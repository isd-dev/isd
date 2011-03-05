/* Copyright 2010 ISD-team */

#ifndef _ISD_ISD_H_
#define _ISD_ISD_H_

#ifdef __cplusplus
#include <cstddef>  // for size_t
#endif

/*/////////
// C API //
/////////*/

typedef struct {
  int n;
  int k;
  int w;
  int px;
  int py;
  int m;
  int l;
  int c;
  int r;
  int kx;
  int ky;
  int max_iteration;
  int max_words;
  unsigned int seed;
  void (*callback)(unsigned long *);
} isd_parameters;

typedef struct {
  int iteration;
  int n_found;
} isd_stats;

#ifdef __cplusplus
extern "C" {
#endif  /* __cplusplus */

  size_t isd_size_in_limbs(size_t nbits);

  void isd_set_bit(unsigned long *word, int i);
  void isd_unset_bit(unsigned long *word, int i);
  void isd_swap_bit(unsigned long *word, int i);
  int  isd_get_bit(const unsigned long *word, int i);

  int isd_decode(const unsigned long **matrix, const unsigned long *word,
                 const isd_parameters *param, isd_stats *stats);

  /* TODO
  int isd_decode_syst(const unsigned long **matrix, const unsigned long *word,
                      const isd_parameters *param, isd_stats *stats);
  */

  int isd_low_weight_word(const unsigned long **matrix,
                          const isd_parameters *param, isd_stats *stats);

  int isd_low_weight_word_syst(const unsigned long **matrix,
                               const isd_parameters *param, isd_stats *stats);

  void isd_fprintf_parameters(FILE *, const isd_parameters *param);

  int isd_set_default_parameters(int n, int k, int w, isd_parameters *param);

#ifdef __cplusplus
}  /* extern "C" */

/*//////////////////////
// Additional C++ API //
//////////////////////*/

#include <string>
#include <vector>
#include <utility>
#include <memory>
#include <set>

#define DISALLOW_COPY_AND_ASSIGN(TypeName)      \
  TypeName(const TypeName&);                    \
  void operator=(const TypeName&)

namespace ISD {

  class Decoder {
 public:
    explicit Decoder(unsigned int seed = 0x5eed);
    ~Decoder(void);
    void SetSizes(int n, int k);
    void SetMatrix(int n, int k, const unsigned long **matrix);
    void SetMatrixFromString(int n, int k, const std::string &s);
    void SetMatrixFromFile(int n, int k, const FILE *f);
    void SetSystematicMatrix(int n, int k, const unsigned long **matrix);
    void SetSystematicMatrixFromString(int n, int k, const std::string &s);
    void SetSystematicMatrixFromFile(int n, int k, const FILE *f);
    bool AppendRow(const unsigned long *row);
    void MatrixRandomize(void);
    void SetParameter_w(int w);
    void SetParameter_px(int p);
    void SetParameter_py(int p);
    void SetParameter_m(int m);
    void SetParameter_l(int l);
    void SetParameter_c(int c);
    void SetParameter_r(int r);
    void SetParameter_kx(int kx);
    void SetParameter_ky(int ky);
    void SetParameter_max_iteration(int maxiter);
    void SetParameter_max_words(int maxwords);
    void SetParameter_callback(void (*callback)(unsigned long *));
    void SetParameters(const isd_parameters &param);
    void GuessParameters(int n, int k);
    const isd_parameters *parameters(void);
    size_t Run(isd_stats *stats);

 private:
    class Impl;
    Impl *impl_;
    DISALLOW_COPY_AND_ASSIGN(Decoder);
  };


}  /* namespace ISD */

#undef DISALLOW_COPY_AND_ASSIGN

#endif  /* __cplusplus */

#endif  /* _ISD_ISD_H_ */
