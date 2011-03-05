/* Copyright 2010 ISD-team */

#ifndef _ISD_SYSTEMATIC_H_
#define _ISD_SYSTEMATIC_H_

#include <cstdio>
#include <cassert>
#include <vector>
#include <string>
#include <algorithm>

#include "gsl/gsl_rng.h"

#include "./fullmatrix.h"

#define DISALLOW_COPY_AND_ASSIGN(TypeName)     \
  TypeName(const TypeName&);                    \
  void operator=(const TypeName&)

namespace ISD {

  class SystematicGeneratorMatrix {
 public:

    SystematicGeneratorMatrix() : num_rows_(0), num_cols_(0),
        limbs_per_row_(0) {
      limbs_ = new(std::nothrow) unsigned long[0];
      columns_ = new(std::nothrow) int[0];

      JointPivotTab_ = NULL;
      JointPivotTab__ = NULL;
      JointPivotMasks_ = NULL;

      SetFormat("0", "1", ",");
    }

    SystematicGeneratorMatrix(const int len, const int dim)
        : num_rows_(dim),
        num_cols_(len - dim),
        limbs_per_row_(IND(len - dim + ulong_bits - 1)) {
      limbs_ = new(std::nothrow) unsigned long[num_rows_ * limbs_per_row_];
      for (int i = 0; i < num_rows_ * limbs_per_row_; ++i) limbs_[i] = 0;

      columns_ = new(std::nothrow) int[len];
      for (int i = 0; i < len; ++i) columns_[i] = i;

      JointPivotTab_ = NULL;
      JointPivotTab__ = NULL;
      JointPivotMasks_ = NULL;

      SetFormat("0", "1", ",");
    }

    void SetFormat(const std::string z,
                   const std::string o,
                   const std::string s) {
      zero_ = z;
      one_ = o;
      separator_ = s;
    }

    ~SystematicGeneratorMatrix() {
      delete[] limbs_;
      delete[] columns_;
      delete[] JointPivotTab__;
      delete[] JointPivotTab_;
      delete[] JointPivotMasks_;
    }

    void Resize(const int len, const int dim) {
      num_rows_ = dim;
      num_cols_ = len - dim;
      limbs_per_row_ = IND(len - dim + ulong_bits - 1);
      delete[] limbs_;
      limbs_ = new(std::nothrow) unsigned long[num_rows_ * limbs_per_row_];
      for (int i = 0; i < num_rows_ * limbs_per_row_; ++i) limbs_[i] = 0;
      delete[] columns_;
      columns_ = new(std::nothrow) int[len];
      for (int i = 0; i < len; ++i) columns_[i] = i;
    }

    void Set(const unsigned long **matrix) {
      for (int i = 0; i < num_rows_; ++i)
        for (int j = 0; j < limbs_per_row_; ++j)
          limb(i, j) = matrix[i][j];
    }

    //    void Randomize(Random *myrand) {
    void Randomize(gsl_rng *rng) {
      unsigned long mask = (static_cast<unsigned long>(1) << OFF(num_cols_)) - 1;
      for (int row = 0; row < num_rows_; ++row) {
        for (int indice = 0; indice < limbs_per_row_; ++indice) {
          unsigned long res = gsl_rng_get(rng);
          if (SIZEOF_UNSIGNED_LONG == 8)
            res = (res << 32) | gsl_rng_get(rng);
          //          limb(row, indice) = myrand->UInt64_t();
          limb(row, indice) = res;
        }
        // clean last limb
        if (mask)
          limb(row, limbs_per_row_ - 1) &= mask;
      }
    }

    int SetFromFullMatrix(FullMatrix *M) {
      M->DoPivoting();
      Resize(M->num_cols(), M->num_rows());
      for (int i = 0; i < num_rows_; ++i)
        for (int j = 0; j < num_cols_; ++j)
          if (M->bit(i, num_rows_ + j))
            SwapBit(i, j);
      for (int j = 0; j < num_cols_; ++j)
        columns_[j] = M->columns_[j];
      return 0;
    }

    int SetFromString(const std::string T,
                      const char zero,
                      const char one,
                      const bool identity_is_present = false) {
      size_t pos = 0;
      for (int row = 0; row < num_rows_; ++row) {
        if (identity_is_present) {
          for (int column = 0; column < row; ++column)
            while (T[pos] != zero) ++pos;
          while (T[pos] != one) ++pos;
          for (int column = row+1; column < num_rows_; ++column)
            while (T[pos] != zero) ++pos;
        }
        for (int column = 0; column < num_cols_; ++column) {
          int indice = IND(column);
          int offset = OFF(column);
          char c;
          do {
            if (pos == T.size())
              return 1;
            c = T[pos];
            ++pos;
          } while ((c != zero) && (c != one));
          if (c == one)
            limb(row, indice) |= static_cast<unsigned long>(1) << offset;
        }
      }
      return 0;
    }

    inline unsigned long bit(const int row, const int column) const {
      const int indice = IND(column);
      const int offset = OFF(column);
      return (limb(row, indice) >> offset) & 1;
    }

    inline void SwapBit(const int row, const int column) const {
      const int indice = IND(column);
      const int offset = OFF(column);
      limb(row, indice) ^= static_cast<unsigned long>(1) << offset;
    }

    inline void UnsetBit(const int row, const int column) const {
      const int indice = IND(column);
      const int offset = OFF(column);
      limb(row, indice) &= ~(static_cast<unsigned long>(1) << offset);
    }

    inline unsigned long GetAndUnsetBit(const int row, const int column) const {
      const int indice = IND(column);
      const int offset = OFF(column);
      unsigned long t = limb(row, indice) & (static_cast<unsigned long>(1) << offset);
      limb(row, indice) ^= t;
      return t >> offset;
    }

    inline void SetBit(const int row, const int column) const {
      const int indice = IND(column);
      const int offset = OFF(column);
      limb(row, indice) |= static_cast<unsigned long>(1) << offset;
    }

    inline unsigned long &limb(const int row, const int indice) const {
      return limbs_[row * limbs_per_row_ + indice ];
    }

    std::string ToString() const {
      std::string s = "(\n";
      std::string end_line = "" + separator_ + "\n";
      std::vector<int> v(1);
      for (int i = 0; i < num_rows_; ++i) {
        v[0] = i;
        s += Code(v) + end_line;
      }
      s += ")";
      return s;
    }

    void Print(void) const {
      printf("Matrix(GF(2), %d, %d, ", dimension(), length());
      printf("%s", ToString().c_str());
      printf(")\n");
      for (int i = 0; i < length(); ++i)
        printf(" %d", columns_[i]/10);
      printf("\n");
      for (int i = 0; i < length(); ++i)
        printf(" %d", columns_[i]%10);
      printf("\n");
    }

    bool AppendRow(const unsigned long *row_) {
      printf("BEFORE:\n\n");
      Print();
      unsigned long *old_limbs_ = limbs_;
      limbs_ = new unsigned long[(num_rows_ + 1) * limbs_per_row_];
      for (int i = 0; i < num_rows_ * limbs_per_row_; ++i)
        limbs_[i] = old_limbs_[i];
      unsigned long *new_row_ = limbs_ + num_rows_ * limbs_per_row_;
      for (int i = 0; i < limbs_per_row_; ++i)
        new_row_[i] = 0;
      // do pivoting
      for (int i = 0; i < length(); ++i) {
        const int indice = IND(i);
        const int offset = OFF(i);
        if ((row_[indice] >> offset) & 1) {
          if (columns_[i] <= num_rows_) {
            // Add(num_rows_, i);
          } else {
            SwapBit(num_rows_, columns_[i] - num_rows_);
          }
        }
      }
      /*
      // can we find a pivot in the added line?
      {
        int i;
        for (i = 0; i < num_cols_; ++i)
          if (bit(num_rows_, i))
            break;
        if (i == num_cols_) {
          // added row was already in the code
          delete[] limbs_;
          limbs_ = old_limbs_;
          return false;
        }
        std::swap(columns_[num_rows_], columns_[num_rows_ + i]);
      }
      */
      num_rows_++;
      num_cols_--;
      printf("AFTER:\n\n");
      Print();
      return true;
    }

    inline void Add(const int row0, const int row1) {
      for (int i = 0; i < limbs_per_row_; ++i)
        limb(row0, i) ^= limb(row1, i);
    }

    inline void SwapColumns(const int in, const int out) {
      assert(0 <= in  && in  < num_rows_);
      assert(0 <= out && out < num_cols_);
      SwapBit(in, out);
      for (int row = 0; row < num_rows_; ++row)
        if (bit(row, out)) Add(row, in);
      SwapBit(in, out);
      std::swap(columns_[in], columns_[num_rows_ + out]);
    }

    void InitJointPivotTabs(int r) {
      delete[] JointPivotTab__;
      delete[] JointPivotTab_;
      delete[] JointPivotMasks_;
      JointPivotTab_ = new(std::nothrow) unsigned long *[1 << r];
      JointPivotTab__ = new(std::nothrow) unsigned long[limbs_per_row_ << r];
      for (int i = 0; i < (limbs_per_row_ << r); ++i)
        JointPivotTab__[i] = 0;
      for (int i = 0; i < (1 << r); ++i)
        JointPivotTab_[i] = &(JointPivotTab__[limbs_per_row_ * i]);
      JointPivotMasks_ = new(std::nothrow) int[r];
    }

    void JointPivot(const int r, const int *in, const int *out) {
      unsigned long **Tab = JointPivotTab_;
      // first r lines
      for (int i = 0; i < r; ++i) {
        int mask = 0;
        const int row = in[i];
        for (int j = 0; j < r; ++j)
          mask |= bit(row, out[j]) << j;
        JointPivotMasks_[i] = mask;
        for (int j = 0; j < limbs_per_row_; ++j)
          Tab[mask][j] = limb(row, j);
      }
      // setting identity matrix
      for (int i = 0; i < r; ++i) {
        for (int j = 0; j < r; ++j) {
          unsigned long mask = ~(static_cast<unsigned long>(1) << OFF(out[j]));
          Tab[JointPivotMasks_[i]][IND(out[j])] &= mask;
        }
        unsigned long mask = static_cast<unsigned long>(1) << OFF(out[i]);
        Tab[JointPivotMasks_[i]][IND(out[i])] |= mask;
      }
      // precomputation with gray code
      {
        int mask = JointPivotMasks_[0];
        for (int i = 2; i < (1 << r); ++i) {
          const int d = count_zeros(i);
          const int next_mask = mask ^ JointPivotMasks_[d];
          if ((i & (i+1))) {
            for (int j = 0; j < limbs_per_row_; ++j) {
              int added_row = JointPivotMasks_[d];
              Tab[next_mask][j] = Tab[mask][j] ^ Tab[added_row][j];
            }
          }
          mask = next_mask;
        }
      }
      // pivoting
      for (int i = 0; i < r; ++i)
        SwapBit(in[i], out[i]);
      for (int row = 0; row < num_rows_; ++row) {
        int mask = 0;
        for (int k = 0; k < r; ++k)
          mask |= GetAndUnsetBit(row, out[k]) << k;
        // if (mask)
        for (int j = 0; j < limbs_per_row_; ++j)
          limb(row, j) ^= Tab[mask][j];
      }
      for (int i = 0; i < r; ++i) {
        SwapBit(in[i], out[i]);
        std::swap(columns_[in[i]], columns_[num_rows_ + out[i]]);
      }
    }

    std::string Code(const std::vector<int> &pos) const {
      std::vector<char> v(length());
      for (size_t i = 0; i < pos.size(); ++i) {
        int pi = pos[i];
        assert(pi < num_rows_);
        v[columns_[pi]] = 1;
        for (int j = 0; j < num_cols_; ++j)
          v[columns_[num_rows_ + j]] ^= bit(pi, j);
      }
      std::string res = "(";
      for (size_t i = 0; i < v.size(); ++i)
        if (v[i])
          res += one_ + separator_;
        else
          res += zero_ + separator_;
      res += ")";
      return res;
    }

    std::vector<unsigned long> Code_v(const std::vector<int> &pos) const {
      int vsize = IND(length() + ulong_bits - 1);
      std::vector<unsigned long> v(vsize);
      for (int i = 0; i < vsize; ++i)
        v[i] = 0;
      for (size_t i = 0; i < pos.size(); ++i) {
        int pi = pos[i];
        assert(pi < num_rows_);
        int col = columns_[pi];
        v[IND(col)] |= 1UL << OFF(col);
        for (int j = 0; j < num_cols_; ++j) {
          int pj = columns_[num_rows_ + j];
          int indice = IND(pj);
          int offset = OFF(pj);
          v[indice] ^= bit(pi, j) << offset;
        }
      }
      return v;
    }

    bool SubmatrixIsInvertible(int *rows,
                               int *cols,
                               const int n) const {
      assert(n < 64);
      static unsigned long m[64];

      for (int i = 0; i < n; ++i) {
        register unsigned long mi = 0;
        for (int j = 0; j < n; ++j)
          mi |= bit(rows[i], cols[j]) << j;
        m[i] = mi;
      }
      // Gauss pivoting
      for (int r = 0; r < n; ++r) {
        {
          int c;
          for (c = r; c < n; ++c)
            if ((m[c] >> r) & 1) break;
          if (c == n) {  // no pivot found
            return false;
          }
          // swap columns
          std::swap(rows[r], rows[c]);
          std::swap(m[r], m[c]);
        }
        // do pivoting
        for (int k = r+1; k < n; ++k)
          if ((m[k] >> r) & 1) m[k] ^= m[r];
      }
      return true;
    }

    // Accessors
    inline int limbs_per_row(void) const {
      return limbs_per_row_;
    }

    inline int num_rows(void) const {
      return num_rows_;
    }

    inline int num_cols(void) const {
      return num_cols_;
    }

    inline const int *columns(void) const {
      return columns_;
    }

    inline int length(void) const {
      return num_rows_ + num_cols_;
    }

    inline int dimension(void) const {
      return num_rows_;
    }

    // DATA

 private:

    unsigned long *limbs_;
    int num_rows_, num_cols_, limbs_per_row_;
    int *columns_;
    unsigned long **JointPivotTab_;
    unsigned long *JointPivotTab__;
    int *JointPivotMasks_;
    std::string zero_, one_, separator_;

    DISALLOW_COPY_AND_ASSIGN(SystematicGeneratorMatrix);
  };

}

#undef DISALLOW_COPY_AND_ASSIGN

#endif  // _ISD_SYSTEMATIC_H_
