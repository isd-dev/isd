/* Copyright 2010 ISD-team */

#ifndef _ISD_DECODE_IMPL_H_
#define _ISD_DECODE_IMPL_H_

#include <vector>
#include <string>
#include <utility>
#include <set>

#include <gsl/gsl_rng.h>

#include "./isd.h"
#include "./radixsorter.h"
#include "./systematic.h"

namespace ISD {

  class Decoder::Impl {
 public:
    // needed for the interface *********
    explicit Impl(unsigned int seed = 0x5eed);
    ~Impl(void);
    size_t Run(isd_stats *stats);
    void SetSizes(int n, int k);
    void MatrixRandomize(void);
    void SetMatrix(int n, int k, const unsigned long ** const matrix);
    void SetMatrixFromString(int n, int k, const std::string &s);
    void SetMatrixFromFile(int n, int k, const FILE *f);
    void SetSystematicMatrix(int n, int k, const unsigned long **matrix);
    void SetSystematicMatrixFromString(int n, int k, const std::string &s);
    void SetSystematicMatrixFromFile(int n, int k, const FILE *f);
    bool AppendRow(const unsigned long *row);
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

 private:
    // do what you want here ************
    // METHODS
    const isd_parameters &parametersImpl(void);
    void CheckParameters(void);
    void RandomSwap(void) throw();
    void FindCollision(void);
    void ChooseXY(void);
    void ChooseZ(void);
    void Collision_gen(void);
    void Collision_p_is_1(void);
    void Collision_p_is_2(void);
    void fill_bHA(const unsigned long *MlimbZ, int offset, unsigned long mask,
                  int level, int maxi, unsigned long partial_sum);
    void fill_collB(const unsigned long *MlimbZ, int offset, unsigned long mask,
                    int level, int maxi, unsigned long partial_sum);
    void fill_collA(const unsigned long *MlimbZ, int offset, unsigned long mask,
                    int level, int maxi, unsigned long partial_sum);
    void Test_weight(int indA, int indB);
    void Test_p_is_1(int indA, int indB);
    void Test_p_is_2(int indA, int indB);

    // DATA
    gsl_rng *myrand_;
    isd_parameters parameters_;
    isd_stats *stats_;
    RadixSorter<8> radixsorter_;
    SystematicGeneratorMatrix M;
    int Z;
    unsigned long *bHA;
    unsigned long *bHB;
    int *v_in_;
    int *v_out_;
    int *X;
    int *Y;
    std::vector<int> vZ;
    std::vector< std::pair<int, int> > collA;
    std::vector< std::pair<int, int> > collB;
    std::vector<int> collAv;
    std::vector<int> collBv;
    std::vector<int> tuple_;
   };

}  // namespace ISD

#endif  // _ISD_DECODE_IMPL_H_
