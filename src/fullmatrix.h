/* Copyright 2010 */

#ifndef _ISD_FULLMATRIX_H_
#define _ISD_FULLMATRIX_H_

#include <cstdio>
#include <memory>
#include <string>
#include <algorithm>

#include "./utils.h"

namespace ISD {

  class FullMatrix {
 public:
    FullMatrix(void)
     : columns_(NULL), num_rows_(0), num_cols_(0), limbs_per_row_(0), rank_(0),
        rows_(NULL)
    {}

    ~FullMatrix(void) {
      Free();
    }

    void Free(void) {
      if (rows_ != NULL) {
        for (int i = 0; i < num_rows_ ; ++i)
          delete[] rows_[i];
      }
      delete[] rows_;
      delete[] columns_;
    }

    void Resize(const int len, const int dim) {
      Free();
      num_rows_ = dim;
      num_cols_ = len;
      limbs_per_row_ = IND(len + ulong_bits - 1);
      rows_ = new(std::nothrow) unsigned long *[num_rows_];
      for (int i = 0; i < num_rows_ ; ++i) {
        rows_[i] = new(std::nothrow) unsigned long[limbs_per_row_];
        for (int j = 0; j < limbs_per_row_; ++j) rows_[i][j] = 0;
      }
      columns_ = new(std::nothrow) int[len];
      for (int i = 0; i < len; ++i) columns_[i] = i;
    }

    void Set(const unsigned long ** const matrix) {
      for (int i = 0; i < num_rows_; ++i)
        for (int j = 0; j < limbs_per_row_; ++j)
          limb(i, j) = matrix[i][j];
    }

    void SetFromFile(const FILE *f) {
      // TODO()
    }

    void SetFromString(const std::string &s) {
      // TODO()
    }

    inline unsigned long &limb(const int row, const int indice) const {
      return rows_[ row ][ indice ];
    }

    inline unsigned long bit(const int row, const int column) const {
      const int indice = IND(column);
      const int offset = OFF(column);
      return (limb(row, indice) >> offset) & 1;
    }

    inline void Add(const int row0, const int row1) {
      for (int i = 0; i < limbs_per_row_; ++i)
        limb(row0, i) ^= limb(row1, i);
    }

    void DoPivoting(void) {
      // going down
      {
        int k = num_rows_;
        int i = 0;
        while ( i < k ) {
          {
            int j = i;
            // look for a pivot
            while (j < num_cols_ && !bit(i, columns_[j])) ++j;
            // no pivot found
            if (j == num_cols_) {
              k--;
              std::swap(rows_[i], rows_[k]);
              continue;
            }
            // pivot found
            std::swap(columns_[i], columns_[j]);
          }
          for (int j = i + 1; j < k; ++j)
            if (bit(j, columns_[i]))
              Add(j, i);
          ++i;
        }
        rank_ = k;
      }
      // going up
      for (int i = rank_ - 1; i; --i)
        for (int j = 0; j < i; ++j)
          if (bit(j, columns_[i]))
            Add(j, i);
    }

    void Print(void) {
      printf("Matrix(GF(2), %d, %d,\n", num_rows_, num_cols_);
      for (int i = 0; i < num_rows_; ++i) {
        printf("(");
        for (int j = 0; j < num_cols_; ++j)
          printf("%d,", static_cast<int>(bit(i, j)));
        printf(")\n");
      }
      printf(")\n");
    }

    const int *columns_order(void) const {
      return columns_;
    }

    //  getters

    int num_rows() {
      return num_rows_;
    }

    int num_cols() {
      return num_cols_;
    }

    int *columns_;
 private:
    int num_rows_;
    int num_cols_;
    int limbs_per_row_;
    int rank_;
    unsigned long **rows_;
  };
}  // namespace ISD

#endif  // _ISD_FULLMATRIX_H_
