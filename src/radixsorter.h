/* Copyright 2010 ISD-team */

#ifndef _ISD_RADIX_SORTER_H_
#define _ISD_RADIX_SORTER_H_

#include <vector>
#include <utility>
#include <algorithm>

namespace ISD {

  template<int kRadixBits>
  class RadixSorter {
 public:
    RadixSorter(void) : count(new int[kSize]) {
    }

    ~RadixSorter(void) {
      delete[] count;
    }

    void sort(std::vector< std::pair<int, int> > *v, int nbits) {
      const size_t vsize = v->size();
      if (V.size() < vsize) {
        V.resize(vsize);
        next.resize(vsize);
      }
      const int mask = kSize - 1;
      // first pass
      for (size_t i = 0; i < kSize; ++i)
        count[i] = 0;
      for (size_t i = 0; i < vsize; ++i)
        count[(*v)[i].first & mask]++;
      for (size_t i = 1; i < kSize; ++i)
        count[i] += count[i-1];
      for (size_t i = 0; i < vsize; ++i) {
        int m = (*v)[i].first & mask;
        count[m]--;
        V[count[m]] = (*v)[i];
      }
      // all but last passes
      {
        int d;
        for (d = kRadixBits; d < nbits - kRadixBits; d += kRadixBits) {
          for (size_t i = 0; i < kSize; ++i)
            count[i] = 0;
          for (size_t i = 0; i < vsize; ++i)
            count[(V[i].first >> d) & mask]++;
          for (size_t i = 1; i < kSize; ++i)
            count[i] += count[i-1];
          for (size_t i = 0; i < vsize; ++i) {
            int m = (V[i].first >> d) & mask;
            count[m]--;
            next[count[m]] = V[i];
          }
          std::swap(V, next);
        }
        // last pass
        for (size_t i = 0; i < kSize; ++i)
          count[i] = 0;
        for (size_t i = 0; i < vsize; ++i)
          count[V[i].first >> d]++;
        for (size_t i = 1; i < kSize; ++i)
          count[i] += count[i-1];
        for (size_t i = 0; i < vsize; ++i) {
          int m = V[i].first >> d;
          count[m]--;
          (*v)[count[m]] = V[i];
        }
      }
    }
 private:
    static const size_t kSize = 1 << kRadixBits;
    int *count;
    std::vector< std::pair<int, int> > V;
    std::vector< std::pair<int, int> > next;
  };
}
#endif  // _ISD_RADIX_SORTER_H_
