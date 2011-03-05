/* Copyright 2010 ISD-team */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <utility>

#include <gsl/gsl_rng.h>

#include "./isd-impl.h"

namespace ISD {
/////////////////////
// Implementations //
/////////////////////

Decoder::Impl::Impl(unsigned int seed)
    : myrand_(gsl_rng_alloc(gsl_rng_gfsr4)) {
  gsl_rng_set(myrand_, seed);
  memset(&parameters_, 0, sizeof(isd_parameters));
}

Decoder::Impl::~Impl(void) {
  gsl_rng_free(myrand_);
}

// setters
//////////

void Decoder::Impl::SetParameter_w(int w) {
  parameters_.w = w;
}

void Decoder::Impl::SetParameter_px(int p) {
  parameters_.px = p;
}

void Decoder::Impl::SetParameter_py(int p) {
  parameters_.py = p;
}

void Decoder::Impl::SetParameter_m(int m) {
  parameters_.m = m;
}

void Decoder::Impl::SetParameter_l(int l) {
  parameters_.l = l;
}

void Decoder::Impl::SetParameter_c(int c) {
  parameters_.c = c;
}

void Decoder::Impl::SetParameter_r(int r) {
  parameters_.r = r;
}

void Decoder::Impl::SetParameter_kx(int kx) {
  parameters_.kx = kx;
}

void Decoder::Impl::SetParameter_ky(int ky) {
  parameters_.ky = ky;
}

void Decoder::Impl::SetParameter_max_iteration(int mi) {
  parameters_.max_iteration = mi;
}

void Decoder::Impl::SetParameter_max_words(int mw) {
  parameters_.max_words = mw;
}

void Decoder::Impl::SetParameter_callback(void (*callback)(unsigned long *)) {
  parameters_.callback = callback;
}

void Decoder::Impl::SetSizes(int n, int k) {
  parameters_.n = n;
  parameters_.k = k;
  M.Resize(n, k);
}

void Decoder::Impl::SetParameters(const isd_parameters &param) {
  parameters_ = param;
}

void Decoder::Impl::SetMatrix
(int n, int k, const unsigned long ** const matrix) {
  SetSizes(n, k);

  ISD::FullMatrix F;
  F.Resize(n, k);
  F.Set(matrix);
  M.SetFromFullMatrix(&F);
}

void Decoder::Impl::SetMatrixFromFile
(int n, int k, const FILE *matrix) {
  SetSizes(n, k);

  ISD::FullMatrix F;
  F.Resize(n, k);
  F.SetFromFile(matrix);
  M.SetFromFullMatrix(&F);
}

void Decoder::Impl::SetMatrixFromString
(int n, int k, const std::string &matrix) {
  SetSizes(n, k);
  // TODO()
}

void Decoder::Impl::SetSystematicMatrix
(int n, int k, const unsigned long **matrix) {
  SetSizes(n, k);
  M.Set(matrix);
}

void Decoder::Impl::SetSystematicMatrixFromFile
(int n, int k, const FILE *matrix) {
  SetSizes(n, k);
  // TODO()
  // M.Set(matrix);
}

void Decoder::Impl::SetSystematicMatrixFromString
(int n, int k, const std::string &matrix) {
  SetSizes(n, k);
  M.SetFromString(matrix, '0', '1');
}

// getters
//////////

const isd_parameters *Decoder::Impl::parameters(void) {
  return &parameters_;
}

// parameters stuff TODO()
///////////////////

void Decoder::Impl::GuessParameters(int n, int k) {
  SetSizes(n, k);
  // TODO(): do it
}

void Decoder::Impl::CheckParameters(void) {
  assert(parameters_.n == M.num_rows() + M.num_cols());
  assert(parameters_.k == M.num_rows());

  if ((parameters_.kx | parameters_.ky) == 0) {
    parameters_.kx = parameters_.k / 2;
    parameters_.ky = parameters_.k - parameters_.kx;
  }

  assert(parameters_.kx >= 0);
  assert(parameters_.kx <= parameters_.k);
  assert(parameters_.ky >= 0);
  assert(parameters_.ky <= parameters_.k);

  if (parameters_.max_iteration == 0)
    parameters_.max_iteration = (1UL << 31) - 1;

  // TODO(): avoid stupid values
}

// ISD algorithm stuff
//////////////////////

size_t Decoder::Impl::Run(isd_stats *stats) {
  CheckParameters();
  // allocation and initialization
  {
    if (stats == NULL) {
      stats_ = new isd_stats;
      stats_->iteration = 0;
      stats_->n_found = 0;
    } else {
      stats_ = stats;
    }
    bHA = new(std::nothrow) unsigned long[IND(1UL << parameters_.l)];
    bHB = new(std::nothrow) unsigned long[IND(1UL << parameters_.l)];
    v_in_ = new(std::nothrow) int[parameters_.k];
    for (int i = 0; i < parameters_.k; ++i) v_in_[i] = i;
    v_out_ = new(std::nothrow) int[parameters_.n - parameters_.k];
    for (int i = 0; i < parameters_.n - parameters_.k; ++i) v_out_[i] = i;
    X = new(std::nothrow) int[parameters_.k];
    for (int i = 0; i < parameters_.k; ++i) X[i] = i;
    Y = X + std::min(parameters_.kx + parameters_.ky, parameters_.k) - parameters_.ky;
    vZ = std::vector<int>(M.limbs_per_row());
    for (int i = 0; i < M.limbs_per_row(); ++i) vZ[i] = i;
    M.InitJointPivotTabs(parameters_.r);
    tuple_.resize(std::max(parameters_.px, parameters_.py));
  }
  // iterations
  while ((++stats_->iteration < parameters_.max_iteration) &&
         (stats_->n_found < parameters_.max_words)) {
    RandomSwap();
    FindCollision();
  }
  // clean-up
  {
    delete[] X;
    delete[] v_out_;
    delete[] v_in_;
    delete[] bHB;
    delete[] bHA;
    if (stats == NULL)
      delete stats_;
  }

  return 0;  // everything went ok
}

void Decoder::Impl::RandomSwap(void) throw() {
  do {
    int up_in = parameters_.k;
    int up_out = parameters_.n - parameters_.k;
    for (int i = 0; i < parameters_.c; ++i, up_in--, up_out--) {
      std::swap(v_in_[i], v_in_[i + gsl_rng_uniform_int(myrand_, up_in)]);
      std::swap(v_out_[i], v_out_[i + gsl_rng_uniform_int(myrand_, up_out)]);
    }
  } while (!M.SubmatrixIsInvertible(v_in_, v_out_, parameters_.c));
  // joint pivoting
  {
    int i = parameters_.c % parameters_.r;
    if (i) M.JointPivot(i, v_in_, v_out_);
    while (i < parameters_.c) {
      M.JointPivot(parameters_.r, v_in_ + i, v_out_ + i);
      i += parameters_.r;
    }
  }
}

void Decoder::Impl::FindCollision(void) {
  ChooseXY();
  ChooseZ();
  for (int i = 0; i < parameters_.m; ++i) {
    Z = vZ[i];
    if (parameters_.px == parameters_.py) {
      switch (parameters_.px) {
        case 1: Collision_p_is_1(); break;
        case 2: Collision_p_is_2(); break;
        default: Collision_gen();
      }
    } else {
      Collision_gen();
    }
  }
}

void Decoder::Impl::ChooseXY(void) {
  int up = std::min(parameters_.kx + parameters_.ky, parameters_.k);
  for (int i = 0; i < parameters_.m; ++i, --up)
    std::swap(X[i], X[i + gsl_rng_uniform_int(myrand_, up)]);
}

void Decoder::Impl::ChooseZ(void) {
  int up = M.limbs_per_row();
  for (int i = 0; i < parameters_.m; ++i, --up)
    std::swap(vZ[i], vZ[i + gsl_rng_uniform_int(myrand_, up)]);
}

void Decoder::Impl::fill_bHA(const unsigned long *MlimbZ,
                             const int offset,
                             const unsigned long mask,
                             int level,
                             int maxi,
                             unsigned long partial_sum) {
  // fill the boolean array bHA with the possible value
  // of the \ell bits of the sum set to 1
  if (--level) {
    // sum not finished, recurse
    for (int i = level; i < maxi; ++i) {
      unsigned long p = partial_sum ^ MlimbZ[X[i]];
      fill_bHA(MlimbZ, offset, mask, level, i, p);
    }
  } else {
    // last xor in the sum, do it and fill. (base case)
    for (int i = level; i < maxi; ++i) {
      unsigned long p = partial_sum ^ MlimbZ[X[i]];
      int s = (p >> offset) & mask;
      bHA[IND(s)] |= static_cast<unsigned long>(1) << OFF(s);
    }
  }
}

void Decoder::Impl::fill_collB(const unsigned long *MlimbZ,
                               const int offset,
                               const unsigned long mask,
                               int level,
                               int maxi,
                               unsigned long partial_sum) {
  // fill the boolean array bHB with \ell bits value of colliding pairs
  // set to one;
  // get the colliding tuples for A
  // keep in collBv the corresponding tuples
  // keep in collB the value of the \ell bits and the index in collBv
  if (--level) {
    // sum not finished, recurse
    for (int i = level; i < maxi; ++i) {
      tuple_[level] = i;
      unsigned long p = partial_sum ^ MlimbZ[Y[i]];
      fill_collB(MlimbZ, offset, mask, level, i, p);
    }
  } else {
    // last xor in the sum, do it and fill. (base case)
    for (int i = level; i < maxi; ++i) {
      tuple_[level] = i;
      unsigned long p = partial_sum ^ MlimbZ[Y[i]];
      int s = (p >> offset) & mask;
      const int is = IND(s);
      const int os = OFF(s);
      const int test = bHA[is] & (static_cast<unsigned long>(1) << os);
      if (test) {
        bHB[is] |= test;
        collB.push_back(std::make_pair(s, collBv.size()));
        collBv.insert(collBv.end(), &tuple_[0], &tuple_[parameters_.py]);
      }
    }
  }
}

void Decoder::Impl::fill_collA(const unsigned long *MlimbZ,
                               const int offset,
                               const unsigned long mask,
                               int level,
                               int maxi,
                               unsigned long partial_sum) {
  // get the colliding tuples for A
  // keep in collAv the corresponding tuples
  // keep in collA the value of the \ell bits and the index in collBv
  if (--level) {
    // sum not finished, recurse
    for (int i = level; i < maxi; ++i) {
      tuple_[level] = i;
      unsigned long p = partial_sum ^ MlimbZ[X[i]];
      fill_collA(MlimbZ, offset, mask, level, i, p);
    }
  } else {
    // last xor in the sum, do it and fill. (base case)
    for (int i = level; i < maxi; ++i) {
      tuple_[level] = i;
      unsigned long p = partial_sum ^ MlimbZ[X[i]];
      int s = (p >> offset) & mask;
      if ((bHB[IND(s)] >> OFF(s)) & 1) {
        collA.push_back(std::make_pair(s, collAv.size()));
        collAv.insert(collAv.end(), &tuple_[0], &tuple_[parameters_.px]);
      }
    }
  }
}

void Decoder::Impl::Collision_gen(void) {
  const int bHsize = IND(1UL << parameters_.l);
  for (int i = 0; i < bHsize; ++i) bHA[i] = 0;
  for (int i = 0; i < bHsize; ++i) bHB[i] = 0;
  // setup mask and offset to select \ell bits
  // for efficiency, selected bits fits in one limb
  const unsigned long mask = (1UL << parameters_.l) - 1;
  const int offset = gsl_rng_uniform_int(myrand_, ulong_bits - parameters_.l);
  const unsigned long *MlimbZ = &M.limb(0, Z);

  // get all possible values for A
  fill_bHA(MlimbZ, offset, mask, parameters_.px, parameters_.kx, 0);
  // check colliding values for B
  collB.clear();
  collBv.clear();
  fill_collB(MlimbZ, offset, mask, parameters_.py, parameters_.ky, 0);
  if (!collBv.empty()) {
    // check colliding values for A
    collA.clear();
    collAv.clear();
    fill_collA(MlimbZ, offset, mask, parameters_.px, parameters_.kx, 0);

    // process colliding pairs
    radixsorter_.sort(&collA, parameters_.l);
    radixsorter_.sort(&collB, parameters_.l);
    // std::sort(collA.begin(), collA.end());
    // std::sort(collB.begin(), collB.end());
    // add a sentinel
    collB.push_back(std::make_pair(collB.back().first ^ 1, 0));
    int j_start = 0;
    int j = 0;
    // iterate through colliding pairs
    for (size_t i = 0; i < collA.size(); ++i) {
      const int sum = collA[i].first;
      const int indXi = collA[i].second;
      if (collB[j_start].first != sum)
        j_start = j;
      else
        j = j_start;
      while (collB[j].first == sum) {
        Test_weight(indXi, collB[j].second);
        ++j;
      }
    }
  }
}

void Decoder::Impl::Collision_p_is_1(void) {
  const int bHsize = IND(1UL << parameters_.l);
  for (int i = 0; i < bHsize; ++i) bHA[i] = 0;
  for (int i = 0; i < bHsize; ++i) bHB[i] = 0;
  // setup mask and offset to select \ell bits
  // for efficiency, selected bits fits in one limb
  unsigned long mask = (1UL << parameters_.l) - 1;
  const int offset = gsl_rng_uniform_int(myrand_, ulong_bits - parameters_.l);
  unsigned long *MlimbZ = &M.limb(0, Z);
  // get all possible values for A
  for (int i = 0; i < parameters_.kx; ++i) {
    unsigned long somme = (MlimbZ[X[i]] >> offset) & mask;
    bHA[IND(somme)] |= 1UL << OFF(somme);
  }
  // check colliding values for B
  // and get indices
  collB.clear();
  for (int i = 0; i < parameters_.ky; ++i) {
    unsigned long somme = (MlimbZ[Y[i]] >> offset) & mask;
    const int is = IND(somme);
    const int os = IND(somme);
    const int test = bHA[is] & (static_cast<unsigned long>(1) << os);
    if (test) {
      bHB[is] |= test;
      collB.push_back(std::make_pair(somme, i));
    }
  }
  if (!collB.empty()) {
    // check colliding values for A
    // and get indices
    collA.clear();
    for (int i = 0; i < parameters_.kx; ++i) {
      int somme = (MlimbZ[X[i]] >>offset) & mask;
      if ((bHB[IND(somme)] >> OFF(somme)) & 1)
        collA.push_back(std::make_pair(somme, i));
    }
    // process colliding pairs
    // radixsorter_.sort(&collA, parameters_.l);
    // radixsorter_.sort(&collB, parameters_.l);
    std::sort(collA.begin(), collA.end());
    std::sort(collB.begin(), collB.end());
    // add sentinelle
    collB.push_back(std::make_pair(collB.back().first ^ 1, 0));
    int j_start = 0;
    int j = 0;
    for (size_t i = 0; i < collA.size(); ++i) {
      const int somme = collA[i].first;
      const int Xi = collA[i].second;
      if (collB[j_start].first != somme)
        j_start = j;
      for (j = j_start; collB[j].first == somme; ++j)
        Test_p_is_1(Xi, collB[j].second);
    }
  }
}

void Decoder::Impl::Collision_p_is_2(void) {
  const int bHsize = IND(1 << parameters_.l);
  for (int i = 0; i < bHsize; ++i) bHA[i] = 0;
  for (int i = 0; i < bHsize; ++i) bHB[i] = 0;
  // setup mask and offset to select \ell bits
  // for efficiency, selected bits fits in one limb
  const unsigned int mask = (1 << parameters_.l) - 1;
  const int offset = gsl_rng_uniform_int(myrand_, ulong_bits - parameters_.l);
  const unsigned long *MlimbZ = &M.limb(0, Z);
  // get all possible values for A
  for (int i = 1; i < parameters_.kx; ++i) {
    const unsigned long MlimbXiZ = MlimbZ[X[i]];
    for (int j = 0; j < i; ++j) {
      const int somme = ((MlimbXiZ ^ MlimbZ[X[j]]) >> offset) & mask;
      bHA[IND(somme)] |= static_cast<unsigned long>(1) << OFF(somme);
    }
  }
  // check colliding values for B
  // and get indices
  collB.clear();
  for (int i = 1; i < parameters_.ky; ++i) {
    const unsigned long MlimbYiZ = MlimbZ[Y[i]];
    for (int j = 0; j < i; ++j) {
      const int somme = ((MlimbYiZ ^ MlimbZ[Y[j]]) >> offset) & mask;
      const int is = IND(somme);
      const int os = OFF(somme);
      const int test = bHA[is] & (static_cast<unsigned long>(1) << os);
      if (test) {
        bHB[is] |= test;
        collB.push_back(std::pair<int, int>(somme, i | (j << 16)));
      }
    }
  }
  if (!collB.empty()) {
    // check colliding values for A
    // and get indices
    collA.clear();
    for (int i = 1; i < parameters_.kx; ++i) {
      const unsigned long MlimbXiZ = MlimbZ[X[i]];
      for (int j = 0; j < i; ++j) {
        const int somme = ((MlimbXiZ ^ MlimbZ[X[j]]) >>offset) & mask;
        if ((bHB[IND(somme)] >> OFF(somme)) & 1)
          collA.push_back(std::pair<int, int>(somme, i | (j << 16)));
      }
    }
    // process colliding pairs
    radixsorter_.sort(&collA, parameters_.l);
    radixsorter_.sort(&collB, parameters_.l);
    /*
    std::sort(collA.begin(), collA.end());
    std::sort(collB.begin(), collB.end());
    */
    // add sentinelle
    collB.push_back(std::pair<int, int>(collB.back().first ^ 1, 0));
    int j_start = 0;
    int j = 0;
    for (size_t i = 0; i < collA.size(); ++i) {
      const int somme = collA[i].first;
      const int Xi = collA[i].second;
      if (collB[j_start].first != somme)
        j_start = j;
      else
        j = j_start;
      while (collB[j].first == somme) {
        Test_p_is_2(Xi, collB[j].second);
        ++j;
      }
    }
  }
}

void Decoder::Impl::Test_weight(int indA, int indB) {
  const int *tupleA = &collAv[indA];
  const int *tupleB = &collBv[indB];

  int weight = 0;

  for (int k = 0; k < M.limbs_per_row(); ++k) {
    unsigned long sum = 0;
    for (int i = 0; i < parameters_.px; ++i) sum ^= M.limb(X[tupleA[i]], k);
    for (int i = 0; i < parameters_.py; ++i) sum ^= M.limb(Y[tupleB[i]], k);
    weight += popcount(sum);
    if (weight > parameters_.w)
      break;
  }

  if (weight <= parameters_.w) {
    std::vector<int> pos;
    for (int i = 0; i < parameters_.px; ++i) pos.push_back(X[tupleA[i]]);
    // remove rows added twice
    for (int i = 0; i < parameters_.py; ++i) {
      int r = Y[tupleB[i]];
      std::vector<int>::iterator j;
      j = std::find(pos.begin(), pos.end(), r);
      if (j == pos.end())
        pos.push_back(r);
      else
        pos.erase(j);
    }
    weight += pos.size();
    if ((weight > 0) && (weight <= parameters_.w)) {
      // a non trivial solution is found, add it
      std::vector<unsigned long> found_word = M.Code_v(pos);
      // it's important to increment before
      // to allow modifications in the callback
      stats_->n_found++;
      if (parameters_.callback != NULL)
        parameters_.callback(&(found_word[0]));
    }
  }
}

void Decoder::Impl::Test_p_is_1(int iA, int iB) {
  int weight = 2;
  const int ml = M.limbs_per_row();
  for (int k = 0; k < ml; ++k) {
    weight += popcount(M.limb(X[iA], k) ^ M.limb(Y[iB], k));
    if (weight > parameters_.w)
      break;
  }

  if (weight <= parameters_.w) {
    // a solution is found, add it
    int p0 = X[iA];
    int p1 = Y[iB];
    if (p0 != p1) {
      std::vector<int> pos(2);
      pos[0] = p0;
      pos[1] = p1;
      std::vector<unsigned long> found_word = M.Code_v(pos);
      // it's important to increment before
      // to allow modifications in the callback
      stats_->n_found++;
      if (parameters_.callback != NULL)
        parameters_.callback(&(found_word[0]));
    }
  }
}

void Decoder::Impl::Test_p_is_2(int cA, int cB) {
  int Ai, Aj;
  Ai = X[cA & 0xFFFF];
  Aj = X[cA >> 16];
  int Bi, Bj;
  Bi = Y[cB & 0xFFFF];
  Bj = Y[cB >> 16];
  int weight = 0;
  for (int k = 0; k < M.limbs_per_row(); ++k) {
    weight += popcount(M.limb(Ai, k) ^ M.limb(Aj, k) ^
                       M.limb(Bi, k) ^ M.limb(Bj, k));
    if (weight > parameters_.w) break;
  }

  if (weight <= parameters_.w) {
    // a solution is found, add it
    std::vector<int> pos(2);
    pos[0] = Ai;
    pos[1] = Aj;
    {
      int r = Bi;
      std::vector<int>::iterator j;
      j = std::find(pos.begin(), pos.end(), r);
      if (j == pos.end())
        pos.push_back(r);
      else
        pos.erase(j);
    }
    {
      int r = Bj;
      std::vector<int>::iterator j;
      j = std::find(pos.begin(), pos.end(), r);
      if (j == pos.end())
        pos.push_back(r);
      else
        pos.erase(j);
    }
    weight += pos.size();
    if ((weight > 0) && (weight <= parameters_.w)) {
      std::vector<unsigned long> found_word = M.Code_v(pos);
      // it's important to increment before
      // to allow modifications in the callback
      stats_->n_found++;
      if (parameters_.callback != NULL)
        parameters_.callback(&(found_word[0]));
    }
  }
}

void Decoder::Impl::MatrixRandomize(void) {
  M.Randomize(myrand_);
}

bool Decoder::Impl::AppendRow(const unsigned long *row) {
  ++parameters_.k;
  return M.AppendRow(row);
}

}  // namespace ISD
