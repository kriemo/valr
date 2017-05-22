#include "rcpphtslib.h"

// bamreader class definition
BamReader::BamReader(const std::string& bampath){
  const char* cbampath = bampath.c_str(); //c-style string pointing to path (must not have ~)
  in = sam_open(cbampath, "rb"); // file object
  if (in == NULL) {
    stop("Failed to open BAM file, check filepath " + bampath);
  }

  bz = in->fp.bgzf ; // bgzf file pointer
  idx = bam_index_load(cbampath); // load BAM index
  header = sam_hdr_read(in) ; // initialize header object

  if (idx == 0) {
    stop("BAM index file is not available for " + bampath);
  }
}

