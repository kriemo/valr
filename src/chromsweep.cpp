// chromsweep.cpp
//
// Copyright (C) 2016 - 2017 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.
// Chromsweep  algorithm rewritten in Rcpp
// https://github.com/arq5x/chrom_sweep/blob/master/chrom_sweep.py

#include "valr.h"

int overlaps(chr_ivl_t a, chr_ivl_t b){
  return std::min(a.stop, b.stop) - std::max(a.start, b.start) ;
}

bool after(chr_ivl_t a, chr_ivl_t b){
  return a.start > b.stop ;
}

void scan_cache(chr_ivl_vector_t::iterator& curr_qy,
                const chr_ivl_vector_t& query,
                chr_ivl_vector_t& db_cache,
                chr_ivl_vector_t& hits){

  if(curr_qy == query.end()){
    return  ;
  }

  chr_ivl_vector_t temp_cache ;

  for (chr_ivl_vector_t::iterator curr_db = db_cache.begin() ;
       curr_db != db_cache.end(); ++curr_db){
    auto past = after(*curr_qy, *curr_db) ;
    if((curr_qy->chrom == curr_db->chrom) && !past){
      temp_cache.push_back(*curr_db) ;
      auto overlap = overlaps(*curr_qy, *curr_db) ;
      if (overlap >= 0) {
        curr_db->value = overlap ;
        hits.push_back(*curr_db) ;
      }
    }
  }
  db_cache = temp_cache ;
}

void report_hits(chr_ivl_t current_qy,
                 chr_ivl_vector_t hits,
                 std::vector<int>& indices_x,
                 chr_ivl_vector_t& y_ivls){

  for(auto it:hits){
    indices_x.push_back(current_qy.value) ;
    y_ivls.push_back(it);
  }
}

chr_ivl_vector_t::iterator get_next(chr_ivl_vector_t::iterator& ivls,
                                    chr_ivl_vector_t container){
  if(ivls != container.end()){
    auto out_it = std::next(ivls, 1) ;
    return out_it ;
  } else {
    return container.end() ;
  }
}

bool check_chrom_order(std::string qy_chrom,
                       std::string db_chrom,
                       const std::vector<std::string>& chroms){

  auto qy_it = std::find(chroms.begin(), chroms.end(), qy_chrom);
  auto db_it = std::find(chroms.begin(), chroms.end(), db_chrom);
  if (qy_it != chroms.end() && db_it != chroms.end())
  {
    auto index = std::distance(qy_it, db_it);
    if(index > 0) {
      // query chrom is less than db chrom
      return true ;
    } else {
      return false ;
    }
  } else
  {
    Rcpp::stop("query or database does not contain a chromosome found in genome file") ;
  }

}

void chrom_check(chr_ivl_vector_t::iterator& curr_qy,
                 chr_ivl_vector_t::iterator& curr_db,
                 const chr_ivl_vector_t& query,
                 const chr_ivl_vector_t& database,
                 chr_ivl_vector_t& db_cache,
                 chr_ivl_vector_t& hits,
                 std::vector<int>& indices_x,
                 chr_ivl_vector_t& y_ivls,
                 const std::vector<std::string>& chroms){

  if(curr_db == database.end() || curr_qy->chrom == curr_db->chrom){
    return ;
  }

  if(!check_chrom_order(curr_qy->chrom, curr_db->chrom, chroms)){
    auto tmp_curr_db = curr_db ;
    while(tmp_curr_db != database.end() &&
          check_chrom_order(tmp_curr_db->chrom, curr_qy->chrom, chroms)){
      ++tmp_curr_db ;
    }
    curr_db = tmp_curr_db ;
    db_cache.clear() ;

  } else if(check_chrom_order(curr_qy->chrom, curr_db->chrom, chroms)){
    auto tmp_curr_qy = curr_qy ;
    while(tmp_curr_qy != query.end() && (tmp_curr_qy->chrom == curr_qy->chrom)){
      scan_cache(tmp_curr_qy, query, db_cache, hits) ;
      report_hits(*tmp_curr_qy, hits, indices_x, y_ivls) ;
      ++tmp_curr_qy ;
      hits.clear() ;
    }
    while(tmp_curr_qy != query.end() &&
          check_chrom_order(tmp_curr_qy->chrom, curr_db->chrom, chroms)){
      report_hits(*tmp_curr_qy, hits, indices_x, y_ivls) ;
      ++tmp_curr_qy ;
    }
    curr_qy = tmp_curr_qy ;
    db_cache.clear() ;
  }
}

void sweep(chr_ivl_vector_t query,
           chr_ivl_vector_t database,
           chr_ivl_vector_t& ivls_y,
           std::vector<int>& indices_x,
           const std::vector<std::string>& chroms) {

  chr_ivl_vector_t db_cache ;
  chr_ivl_vector_t hits ;

  // initialize iterators for query and db
  auto curr_qy = query.begin() ;
  auto curr_db = database.begin() ;

  while(curr_qy != query.end()){
    //chrom check modifies curr_qy, could be advanced to end of vector
    chrom_check(curr_qy, curr_db, query, database,
                db_cache, hits, indices_x, ivls_y, chroms) ;
    //scan_cache checks for end of vector
    scan_cache(curr_qy, query, db_cache, hits) ;
    while (curr_db != database.end() && curr_qy != query.end() &&
           curr_qy->chrom == curr_db->chrom && !after(*curr_db, *curr_qy)){
      auto overlap = overlaps(*curr_qy, *curr_db) ;
      if(overlap >= 0){
        curr_db->value = overlap ;
        hits.push_back(*curr_db) ;
      }
      db_cache.push_back(*curr_db) ;
      ++curr_db ;
    }
    if(curr_qy != query.end()){
      report_hits(*curr_qy, hits, indices_x, ivls_y) ;
      hits.clear() ;
      ++curr_qy ;
    }
  }
}


// [[Rcpp::export]]
DataFrame chromsweep_impl(DataFrame query_df,
                          DataFrame database_df,
                          DataFrame chroms){

  auto query_ivls = makeChromIntervalVector(query_df) ;
  auto database_ivls = makeChromIntervalVector(database_df) ;

  std::vector<int> indices_x ;
  chr_ivl_vector_t ivls_y ;
  std::vector<std::string > chrom_order = chroms["chrom"] ;

  sweep(query_ivls, database_ivls, ivls_y, indices_x, chrom_order) ;

  DataFrame subset_x = DataFrameSubsetVisitors(query_df, names(query_df)).subset(indices_x, "data.frame");
  auto ncol_x = subset_x.size() ;
  // update the y subsets once bam support is setup
  auto ncol_y = subset_x.size() ;

  CharacterVector names(ncol_x + ncol_y) ;
  CharacterVector names_x = subset_x.attr("names") ;
  CharacterVector names_y = names_x ;

  // replacing y chrom with overlap, same number of cols
  List out(ncol_x + ncol_y) ;

  // x names, data
  for (int i = 0; i < ncol_x; i++) {
    auto name_x = as<std::string>(names_x[i]) ;
    if (name_x != "chrom") {
      name_x += ".x" ;
    }
    names[i] = name_x ;
    out[i] = subset_x[i] ;
  }

  // y names, data
  std::vector<int> y_overlaps ;
  std::vector<int> y_starts ;
  std::vector<int> y_stops ;

  for(auto it:ivls_y){
    y_overlaps.push_back(it.value) ;
    y_starts.push_back(it.start) ;
    y_stops.push_back(it.stop) ;
  }

  for (int i = 0; i < ncol_y; i++) {
    auto name_y = as<std::string>(names_y[i]) ;
    if (name_y == "chrom") continue ;
    name_y += ".y" ;
    names[i + ncol_x - 1] = name_y ;
  }

  out[ncol_x + ncol_y - 3] = y_starts ;
  out[ncol_x + ncol_y - 2] = y_stops ;

  // overlaps
  out[ncol_x + ncol_y - 1] = y_overlaps ;
  names[ncol_x + ncol_y - 1] = ".overlap" ;

  out.attr("names") = names ;
  out.attr("class") = classes_not_grouped() ;
  auto nrows = subset_x.nrows() ;
  set_rownames(out, nrows) ;

  return out ;
}

/***R
library(valr)
library(dplyr)

genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 1e6,
    "chr2", 1e7
)

n <- 1e5
x <- bed_random(genome, n = n)
y <- bed_random(genome, n = n)

library(microbenchmark)
microbenchmark(
  valr::chromsweep_impl(x, y, genome),
  bed_intersect(x, y),
  times = 3,
  unit = 's'
)
*/
