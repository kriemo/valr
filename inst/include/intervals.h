// intervals.h
//
// Copyright (C) 2016 - 2017 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#ifndef valr__intervals_H
#define valr__intervals_H

#include "valr.h"

// main interval types used in valr
typedef Interval<int>      ivl_t ;
typedef std::vector<ivl_t> ivl_vector_t ;
typedef IntervalTree<int>  ivl_tree_t ;

// the value field of intervals in the returned vector correspond to the index
// of the interval in the original dataframe (i.e., the values of the
// SlicingIndex)

inline ivl_vector_t makeIntervalVector(DataFrame df, SlicingIndex si,
                                       std::string col_start = "start",
                                       std::string col_end = "end") {

  ivl_vector_t ivls ;

  IntegerVector starts = df[col_start] ;
  IntegerVector ends   = df[col_end] ;

  int size = si.size() ;

  for (int i = 0; i < size; ++i) {
    int j = si[i] ;
    ivls.push_back(ivl_t(starts[j], ends[j], j)) ;
  }
  return ivls ;
}

// Interval class that keeps chromosome and row_number as value
template <class T, typename K = int>
class ChromInterval {
public:
  std::string chrom;
  K start;
  K stop;
  T value;
  ChromInterval(const std::string& c, K s, K e, const T& v)
    : chrom(c)
    , start(s)
    , stop(e)
    , value(v)
  { }
};

// secondary interval types used for chromsweep algorithm
typedef ChromInterval<int> chr_ivl_t ;
typedef std::vector<chr_ivl_t> chr_ivl_vector_t ;

inline chr_ivl_vector_t makeChromIntervalVector(DataFrame df,
                                              std::string col_chrom = "chrom",
                                              std::string col_start = "start",
                                              std::string col_end = "end") {

  chr_ivl_vector_t ivls ;

  std::vector<std::string> chroms = df[col_chrom] ;
  IntegerVector starts = df[col_start] ;
  IntegerVector ends   = df[col_end] ;

  size_t size = df.nrows() ;

  for (int i = 0; i < size; ++i) {
    ivls.push_back(chr_ivl_t(chroms[i], starts[i], ends[i], i)) ;
  }
  return ivls ;
}

#endif
