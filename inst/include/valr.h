// valr.h
//
// Copyright (C) 2016 - 2017 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#ifndef valr__valr_H
#define valr__valr_H

// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
using namespace Rcpp ;

#include <dplyr.h>
using namespace dplyr ;

#include "IntervalTree.h"
#include "intervals.h"
#include "group_apply.h"
#include "genome.h"
#include "random.h"
#include "DataFrameBuilder.h"

#include "hts.h"
#include "sam.h"
#include "bgzf.h"

// class for handling bam file opening and closing
class BamReader {
public:
  samFile* in;
  hts_idx_t* idx;
  BGZF* bz;
  bam_hdr_t* header ;
  BamReader(const std::string& bampath) ;

  ~BamReader(){
    hts_idx_destroy(idx);
    sam_close(in);
  }
};


#endif
