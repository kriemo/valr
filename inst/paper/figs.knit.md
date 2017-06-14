---
title: "F1000 paper figures"
author: "Kent Riemondy RBI"
date: "June 6, 2017"
---




```r
library(valr)
library(tidyverse)
library(cowplot)
```


## Figure 1


```r

# theme to make pretty for publication

pub_theme <- theme(axis.text = element_text(size = 20),
                   axis.title = element_text(size = 20),
                   plot.title = element_text(size = 20),
                   strip.text = element_text(size = 20))

x <- tibble::tribble(
  ~chrom, ~start, ~end,
  'chr1', 25,     50,
  'chr1', 100,    125
)

y <- tibble::tribble(
  ~chrom, ~start, ~end,
  'chr1', 30,     75
)

a <- bed_glyph(bed_intersect(x, y)) + pub_theme

x <- tibble::tribble(
  ~chrom, ~start, ~end,
  'chr1',      1,      50,
  'chr1',      10,     75,
  'chr1',      100,    120
)

b <- bed_glyph(bed_merge(x)) + pub_theme
plot_grid(a, b, align = 'h', nrow = 1, labels="AUTO", scale = c(.9, .9))
#> Warning in align_plots(plotlist = plots, align = align): Graphs cannot be
#> horizontally aligned. Placing graphs unaligned.
```

<img src="figs_files/figure-html/fig1-1.png" width="672" />

```r
ggsave("figure1.pdf", width = 11.2, height = 5.6)
```

## Figure 2

```r
# `valr_example()` identifies the path of example files
bedfile <- valr_example('genes.hg19.chr22.bed.gz')
genomefile <- valr_example('hg19.chrom.sizes.gz')
bgfile  <- valr_example('hela.h3k4.chip.bg.gz')

genes <- read_bed(bedfile, n_fields = 6)
genome <- read_genome(genomefile)
y <- read_bedgraph(bgfile)

# generate 1 bp TSS intervals, `+` strand only
tss <- genes %>%
  filter(strand == '+') %>%
  mutate(end = start + 1)

# 1000 bp up and downstream
region_size <- 1000
# 50 bp windows
win_size <- 50

# add slop to the TSS, break into windows and add a group
x <- tss %>%
  bed_slop(genome, both = region_size) %>%
  bed_makewindows(genome, win_size)

x
#> # A tibble: 13,530 x 7
#>    chrom    start      end      name score strand .win_id
#>    <chr>    <int>    <int>     <chr> <chr>  <chr>   <int>
#>  1 chr22 16161065 16161115 LINC00516     3      +       1
#>  2 chr22 16161115 16161165 LINC00516     3      +       2
#>  3 chr22 16161165 16161215 LINC00516     3      +       3
#>  4 chr22 16161215 16161265 LINC00516     3      +       4
#>  5 chr22 16161265 16161315 LINC00516     3      +       5
#>  6 chr22 16161315 16161365 LINC00516     3      +       6
#>  7 chr22 16161365 16161415 LINC00516     3      +       7
#>  8 chr22 16161415 16161465 LINC00516     3      +       8
#>  9 chr22 16161465 16161515 LINC00516     3      +       9
#> 10 chr22 16161515 16161565 LINC00516     3      +      10
#> # ... with 13,520 more rows

# map signals to TSS regions and calculate summary statistics.
res <- bed_map(x, y, win_sum = sum(value, na.rm = TRUE)) %>%
  group_by(.win_id) %>%
  summarize(win_mean = mean(win_sum, na.rm = TRUE),
            win_sd = sd(win_sum, na.rm = TRUE))

res
#> # A tibble: 41 x 3
#>    .win_id win_mean    win_sd
#>      <int>    <dbl>     <dbl>
#>  1       1 100.8974  85.83423
#>  2       2 110.6829  81.13521
#>  3       3 122.9070  99.09635
#>  4       4 116.2800  96.30098
#>  5       5 116.3500 102.33773
#>  6       6 124.9048  95.08887
#>  7       7 122.9437  94.39792
#>  8       8 127.5946  91.47407
#>  9       9 130.2051  95.71309
#> 10      10 130.1220  88.82809
#> # ... with 31 more rows

x_labels <- seq(-region_size, region_size, by = win_size * 5)
x_breaks <- seq(1, 41, by = 5)

sd_limits <- aes(ymax = win_mean + win_sd, ymin = win_mean - win_sd)

p <- ggplot(res, aes(x = .win_id, y = win_mean)) +
  geom_point(size = 0.25) + geom_pointrange(sd_limits, size = 0.1) + 
  scale_x_continuous(labels = x_labels, breaks = x_breaks) + 
  xlab('Position (bp from TSS)') + ylab('Signal') + 
  theme_classic() +
  pub_theme

plot_grid(p, align = 'h', nrow = 1, labels="AUTO")
```

<img src="figs_files/figure-html/fig2-1.png" width="672" />

```r
ggsave("figure2.pdf", width = 5.6, height = 4)
```

## Figure 3


```r
# Load required packages
library(valr)
library(tidyverse)
library(scales)
library(microbenchmark)
library(stringr)
library(cowplot)

# Set-up test genome
genome <- read_genome(valr_example('hg19.chrom.sizes.gz'))
write_tsv(genome, "genome.txt", col_names=FALSE)

# number of intervals
n <- 1e6
# number of timing reps
nrep <- 10

seed_x <- 1010486
x <- bed_random(genome, n = n, seed = seed_x)
write_tsv(x, "x.bed", col_names=FALSE)
seed_y <- 9283019
y <- bed_random(genome, n = n, seed = seed_y)
write_tsv(y, "y.bed", col_names=FALSE)

res <- microbenchmark(
  # randomizing functions
  bed_random(genome, n = n, seed = seed_x),
  bed_shuffle(x, genome, seed = seed_x),
  # # single tbl functions
  bed_slop(x, genome, both = 1000),
  bed_flank(x, genome, both = 1000),
  bed_merge(x),
  bed_cluster(x),
  bed_complement(x, genome),
  # multi tbl functions
  bed_closest(x, y),
  bed_intersect(x, y),
  bed_map(x, y, .n = n()),
  bed_subtract(x, y),
  # stats
  bed_reldist(x, y),
  bed_jaccard(x, y),
  bed_fisher(x, y, genome),
  #utils 
  bed_makewindows(x, genome, win_size = 100),
  times = nrep,
  unit = 's')

load_file <- function(filename) {
  if(endsWith(filename, ".txt")){
    read_genome(filename)
  } else if (endsWith(filename, ".bed")) {
    read_bed(filename)
  }
}

genome_file <- "genome.txt"
x_bed <- "x.bed"
y_bed  <- "y.bed"

res_hard <- microbenchmark(
  # randomizing functions
  bed_random(load_file(genome_file), n = n, seed = seed_x),
  bed_shuffle(load_file(x_bed), load_file(genome_file), seed = seed_x),
  # # single tbl functions
  bed_slop(load_file(x_bed), load_file(genome_file), both = 1000),
  bed_flank(load_file(x_bed), load_file(genome_file), both = 1000),
  bed_merge(load_file(x_bed)),
  bed_cluster(load_file(x_bed)),
  bed_complement(load_file(x_bed), load_file(genome_file)),
  # multi tbl functions
  bed_closest(load_file(x_bed), load_file(y_bed)),
  bed_intersect(load_file(x_bed), load_file(y_bed)),
  bed_map(load_file(x_bed), load_file(y_bed), .n = n()),
  bed_subtract(load_file(x_bed), load_file(y_bed)),
  # stats
  bed_reldist(load_file(x_bed), load_file(y_bed)),
  bed_jaccard(load_file(x_bed), load_file(y_bed)),
  bed_fisher(load_file(x_bed), load_file(y_bed), load_file(genome_file)),
  #utils 
  bed_makewindows(load_file(x_bed), load_file(genome_file), win_size = 100),
  times = nrep,
  unit = 's')

# covert nanoseconds to seconds
res <- res %>%
  as_tibble() %>%
  mutate(time = time / 1e9)

res_hard <- res_hard %>%
  as_tibble() %>%
  mutate(time = time / 1e9)

# futz with the x-axis
sts <- boxplot.stats(res$time)$stats
# filter out outliers (currently filtering fisher and map)
# res <- filter(res, time <= max(sts) * 1.05)
```


```bash
# Function to benchmark any bash command
repeats=10

bench_tests() {
    # --------------------------------------------------------------------------
    # Benchmark loop
    # --------------------------------------------------------------------------
    echo 'Benchmarking ' $1'...';

    # Run the given command [repeats] times
    for (( i = 1; i <= $repeats ; i++ ))
    do
        # Indicate the command we just ran in the csv file
        echo -ne $3"\t" >> $2
        # runs time function for the called script, output in a comma seperated
        # format output file specified with -o command and -a specifies append
        # requires gnu-time not BSD
        gtime -f %e -o ${2} -a ${1} > /dev/null 2>&1
    done;

}

# bedtools time
# set-up files for bedtools function
bedtools sort -i y.bed > y1.bed
bedtools sort -i x.bed > x1.bed
sort -k1,1 genome.txt > genome1.txt

## Run timing test
rm -f time.txt

bench_tests "bedtools random -g genome1.txt -n 1000000 -seed 1010486" "time.txt" "random"
bench_tests "bedtools shuffle -i x.bed -g genome1.txt -seed 1010486" "time.txt" "shuffle"
# single tbl function
bench_tests "bedtools slop -i x1.bed -g genome1.txt -b 1000" "time.txt" "slop"
bench_tests "bedtools flank -i x1.bed -g genome1.txt -b 1000" "time.txt" "flank"
bench_tests "bedtools merge -i x1.bed" "time.txt" "merge"
bench_tests "bedtools cluster -i x1.bed" "time.txt" "cluster"
bench_tests "bedtools complement -i x1.bed -g genome1.txt" "time.txt" "complement"
# multi tbl functions
bench_tests "bedtools closest -a x1.bed -b y1.bed" "time.txt" "closest"
bench_tests "bedtools intersect -a x1.bed -b y1.bed" "time.txt" "intersect"
bench_tests "bedtools map -a x1.bed -b y1.bed -c 1 -o count" "time.txt" "map"
bench_tests "bedtools subtract -a x1.bed -b y1.bed" "time.txt" "subtract"
# stats
bench_tests "bedtools reldist -a x1.bed -b y1.bed" "time.txt" "reldist"
bench_tests "bedtools jaccard -a x1.bed -b y1.bed" "time.txt" "jaccard"
bench_tests "bedtools fisher -a x1.bed -b y1.bed -g genome1.txt" "time.txt" "fisher"
#utils 
bench_tests "bedtools makewindows -b x1.bed -w 100" "time.txt" "makewindows"

# cleanup data
rm -f x.bed y.bed x1.bed y1.bed genome1.txt genome.txt
```


```r
# Code to process data from previous chunks (above)

# This is to import the data from bedtools (bash) and graph it
times <- read_tsv("time.txt", col_names = c("expr", "time"))
#> Parsed with column specification:
#> cols(
#>   expr = col_character(),
#>   time = col_double()
#> )
dat <- list(times, res_hard)
names(dat) <- c("BEDTools", "valr")
dat <- bind_rows(dat, .id = "package")
#> Warning in bind_rows_(x, .id): binding character and factor vector,
#> coercing into character vector

#mean_dat <- group_by(dat, expr, package) %>% 
#             summarize(mean = mean(time), se = sd(time))

dat$func <- ifelse(dat$package == "BEDTools", 
                   dat$expr, 
                   str_split(dat$expr, "\\(",  simplify = T)[, 1] %>% 
  str_split(., "_", simplify = T) %>% 
  .[,2])

# remove makewindows, really throws off proportion of graph
dat <- filter(dat, func != "makewindows")
valr_data <- filter(dat, package == "valr")

plot_a <- ggplot(res, aes(x=reorder(expr, time), y = time)) +
  geom_boxplot(width = 0.75, size = 0.33, outlier.shape = NA, fill = RColorBrewer::brewer.pal(3, "Set1")[2]) +
  coord_flip() +
  theme_bw() +
  labs(
    y='execution time (seconds)',
    x='',
    subtitle=paste0(comma(n), " random x/y intervals,\n", comma(nrep), " repetitions")) +
  theme_classic() +
  pub_theme +
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(size = 18),
    axis.text.x = element_text(angle = 90)
  )


plot_b <- ggplot(dat, aes(x=reorder(func, time), y = time, fill = package)) +
  geom_boxplot(width = 0.75, size = 0.33, outlier.shape = NA) +
  scale_fill_brewer(palette = "Set1") + 
  coord_flip() +
  theme_bw() +
  labs(
    y='execution time (seconds)',
    x='') +
  theme_classic() +
  pub_theme +
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 90)
  )

plot_grid(plot_a, plot_b, align = 'h', nrow = 1, labels="AUTO", label_size = 28)
#> Warning in align_plots(plotlist = plots, align = align): Graphs cannot be
#> horizontally aligned. Placing graphs unaligned.
```

<img src="figs_files/figure-html/unnamed-chunk-3-1.png" width="672" />

```r
ggsave("figure3.pdf", width = 20, height = 10)
```



## Code demonstrations


```r
library(valr)
# function to retrieve path to example data
bed_filepath <- valr_example("3fields.bed.gz") 
read_bed(bed_filepath) 
#> # A tibble: 10 x 3
#>    chrom  start    end
#>    <chr>  <int>  <int>
#>  1  chr1  11873  14409
#>  2  chr1  14361  19759
#>  3  chr1  14406  29370
#>  4  chr1  34610  36081
#>  5  chr1  69090  70008
#>  6  chr1 134772 140566
#>  7  chr1 321083 321115
#>  8  chr1 321145 321207
#>  9  chr1 322036 326938
#> 10  chr1 327545 328439

#using URL
read_bed("https://github.com/rnabioco/valr/raw/master/inst/extdata/3fields.bed.gz")
#> # A tibble: 10 x 3
#>    chrom  start    end
#>    <chr>  <int>  <int>
#>  1  chr1  11873  14409
#>  2  chr1  14361  19759
#>  3  chr1  14406  29370
#>  4  chr1  34610  36081
#>  5  chr1  69090  70008
#>  6  chr1 134772 140566
#>  7  chr1 321083 321115
#>  8  chr1 321145 321207
#>  9  chr1 322036 326938
#> 10  chr1 327545 328439
```

Using mysql

```r
library(dplyr, warn.conflicts = F)
hg19 <- db_ucsc("hg19")
tbl(hg19, "refGene")
#> Warning in .local(conn, statement, ...): Unsigned INTEGER in col 0 imported
#> as numeric
#> Warning in .local(conn, statement, ...): Unsigned INTEGER in col 4 imported
#> as numeric
#> Warning in .local(conn, statement, ...): Unsigned INTEGER in col 5 imported
#> as numeric
#> Warning in .local(conn, statement, ...): Unsigned INTEGER in col 6 imported
#> as numeric
#> Warning in .local(conn, statement, ...): Unsigned INTEGER in col 7 imported
#> as numeric
#> Warning in .local(conn, statement, ...): Unsigned INTEGER in col 8 imported
#> as numeric
#> Warning in .local(conn, statement, ...): Unsigned INTEGER in col 0 imported
#> as numeric
#> Warning in .local(conn, statement, ...): Unsigned INTEGER in col 4 imported
#> as numeric
#> Warning in .local(conn, statement, ...): Unsigned INTEGER in col 5 imported
#> as numeric
#> Warning in .local(conn, statement, ...): Unsigned INTEGER in col 6 imported
#> as numeric
#> Warning in .local(conn, statement, ...): Unsigned INTEGER in col 7 imported
#> as numeric
#> Warning in .local(conn, statement, ...): Unsigned INTEGER in col 8 imported
#> as numeric
#> # Source:   table<refGene> [?? x 16]
#> # Database: mysql 5.6.26-log [genomep@genome-mysql.cse.ucsc.edu:/hg19]
#>      bin         name chrom strand   txStart     txEnd  cdsStart    cdsEnd
#>    <dbl>        <chr> <chr>  <chr>     <dbl>     <dbl>     <dbl>     <dbl>
#>  1   585    NR_024540  chr1      -     14361     29370     29370     29370
#>  2  1759    NR_110022  chrX      + 153991016 154005964 154005964 154005964
#>  3   604    NR_033380  chrX      +   2527305   2575270   2575270   2575270
#>  4    75    NR_033381  chrY      +   2477305   2506370   2506370   2506370
#>  5   915    NM_182508 chr13      +  43355685  43365685  43358203  43362926
#>  6   923    NM_181845 chr19      +  44331472  44353050  44337656  44352793
#>  7   586    NR_039983  chr1      -    134772    140566    140566    140566
#>  8  1964    NR_028327  chr5      + 180750506 180755196 180755196 180755196
#>  9   126 NM_001271875 chr17      -  55938603  56032684  55943837  55962925
#> 10  2278    NR_125989  chr1      - 222001007 222014008 222014008 222014008
#> # ... with more rows, and 8 more variables: exonCount <dbl>,
#> #   exonStarts <chr>, exonEnds <chr>, score <int>, name2 <chr>,
#> #   cdsStartStat <chr>, cdsEndStat <chr>, exonFrames <chr>
```


### Interval statistics

Estimates of significance for interval overlaps can be obtained by combining `bed_shuffle()`, `bed_random()` and the `sample_` functions from `dplyr` with interval statistics in `valr`.

Here we examine the extent of overlap of repeat classes with exons in the human genome (on `chr22` only, for simplicity) using the jaccard similarity index. `bed_jaccard()` implements the jaccard test to examine the similarity between two sets of genomic intervals. Using `bed_shuffle()` and `replicate()` we generate a `data_frame` containing 100 sets of randomly placed intervals then calculate the jaccard index for each set to generate a null-distribution of jaccard scores. Finally an empirical p-value is then calculated from the null-distribution.


```r
library(tidyverse, warn.conflicts = F)

repeats <- read_bed(valr_example('hg19.rmsk.chr22.bed.gz'), n_fields = 6) 
genome <- read_genome(valr_example('hg19.chrom.sizes.gz'))
genes <- read_bed12(valr_example('hg19.refGene.chr22.bed.gz'))
# convert bed12 to bed with exons
exons <- bed12_to_exons(genes)

# function to repeat interval shuffling
shuffle_intervals <- function(n, .data, genome) {
  replicate(n, bed_shuffle(.data, genome), simplify = FALSE) %>%
    bind_rows(.id = 'rep') %>%
    group_by(rep) %>% nest()
}
nreps <- 100
shuffled <- shuffle_intervals(n = nreps, repeats, genome) %>%
  mutate(jaccard = data %>%
           map(bed_jaccard, repeats) %>%
           map_dbl("jaccard"))
shuffled
#> # A tibble: 100 x 3
#>      rep                  data      jaccard
#>    <chr>                <list>        <dbl>
#>  1     1 <tibble [10,000 x 6]> 0.0004337515
#>  2     2 <tibble [10,000 x 6]> 0.0003501136
#>  3     3 <tibble [10,000 x 6]> 0.0006080318
#>  4     4 <tibble [10,000 x 6]> 0.0005768465
#>  5     5 <tibble [10,000 x 6]> 0.0002537739
#>  6     6 <tibble [10,000 x 6]> 0.0004480599
#>  7     7 <tibble [10,000 x 6]> 0.0001777942
#>  8     8 <tibble [10,000 x 6]> 0.0003448734
#>  9     9 <tibble [10,000 x 6]> 0.0003303327
#> 10    10 <tibble [10,000 x 6]> 0.0003312376
#> # ... with 90 more rows

obs <- bed_jaccard(repeats, exons)
obs
#> # A tibble: 1 x 4
#>    len_i   len_u    jaccard     n
#>    <dbl>   <dbl>      <dbl> <dbl>
#> 1 112123 4132109 0.02789139   805

pvalue <- sum(shuffled$jaccard >= obs$jaccard) + 1 /(nreps + 1)
pvalue
#> [1] 0.00990099
```


## Correlations among DNase I hypersensitive sites 

Here we use `bed_jaccard` for a large-scale comparison of related datasets. As shown in the [`BEDtools` tutorial](http://quinlanlab.org/tutorials/bedtools/bedtools.html), we can measure the similarity of DNaseI hypersensitivity sites for 20 fetal tissue samples. 

This data was taken from [Maurano *et al.* Systematic Localization of Common Disease-Associated Variation in Regulatory DNA. (2012) *Science*](www.sciencemag.org/content/337/6099/1190.short).


```r
library(valr)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(stringr)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
```

First read all 20 BED files containing DNase I hypersensitivity sites from 20 fetal tissues.


```r
dnase_files <- list.files('data/dnasei',pattern = 'merge.bed.gz', full.names = TRUE)
data <- dnase_files %>% map(read_bed, n_fields = 4)
```

Then generate a 20x20 table containing a Jaccard statistic for each of the 400 pairwise comparisons.


```r
res <- data %>%
  cross2(.,.) %>%
  map(lift(bed_jaccard)) %>% 
  map("jaccard") %>%
  flatten_dbl() %>%
  matrix(nrow = 20, ncol = 20)
```

We also need to generate labels for the table from the file names.


```r
# names are tissue + sample_num
col_names <- dnase_files %>%
  str_split('/') %>% map(`[[`, 3) %>%
  str_split('\\.') %>% map(`[[`, 1) %>%
  str_split('-') %>% map(`[[`, 1) %>%
  flatten_chr() %>%
  str_replace('^f', '') %>%
  str_c(str_c('-', seq(length(.))))

colnames(res) <- col_names
rownames(res) <- col_names
```

Now the Jaccard coefficients can be visualized in heatmap form.


```r
Heatmap(res, color_space = 'Blues')
```

Finally we can do some PCA analysis on the Jaccard coefficients to identify clusters of related cell types.


```r
pca <- broom::tidy(prcomp(res)) %>% as_data_frame()

pca_comps <- filter(pca, PC <= 2) %>%
  tidyr::spread(PC, value) %>%
  setNames(c('label','PC1','PC2'))

ggplot(pca_comps, aes(PC1, PC2)) +
  geom_point(size = 3, color = 'red') +
  geom_text_repel(aes(label = label))
```
