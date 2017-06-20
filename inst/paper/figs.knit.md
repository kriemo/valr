---
title: "F1000 paper figures"
author: "Kent Riemondy RBI"
date: "June 6, 2017"
---




```r
library(valr)
library(tidyverse)
#> Loading tidyverse: ggplot2
#> Loading tidyverse: tibble
#> Loading tidyverse: tidyr
#> Loading tidyverse: readr
#> Loading tidyverse: purrr
#> Loading tidyverse: dplyr
#> Conflicts with tidy packages ----------------------------------------------
#> filter(): dplyr, stats
#> lag():    dplyr, stats
library(cowplot)
#> 
#> Attaching package: 'cowplot'
#> The following object is masked from 'package:ggplot2':
#> 
#>     ggsave
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
bedfile <- valr_example("genes.hg19.chr22.bed.gz")
genomefile <- valr_example("hg19.chrom.sizes.gz")
bgfile  <- valr_example("hela.h3k4.chip.bg.gz")

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
  xlab("Position (bp from TSS)") + ylab("Signal") + 
  theme_classic() +
  pub_theme +
  theme(axis.text.x = element_text(angle = 90))

plot_grid(p, align = "h", nrow = 1, labels="AUTO")
```

<img src="figs_files/figure-html/fig2-1.png" width="672" />

```r
ggsave("figure2.pdf", width = 5.6, height = 5.6)
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
genome <- read_genome(valr_example("hg19.chrom.sizes.gz"))
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
  geom_boxplot(width = 0.75, size = 0.33, outlier.shape = NA, fill = grDevices::grey.colors(2)[2]) +
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
  scale_fill_grey() + 
  coord_flip() +
  geom_vline(xintercept = seq_along(unique(dat$func)) + 0.5,
             color = "grey") +
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

### file I/O

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

### syntax demo


```r
library(valr)
library(dplyr)

snps <- read_bed(valr_example("hg19.snps147.chr22.bed.gz"), n_fields = 6)
genes <- read_bed(valr_example("genes.hg19.chr22.bed.gz"), n_fields = 6)

# find snps in intergenic regions
intergenic <- bed_subtract(snps, genes)
# distance from intergenic snps to nearest gene
nearby <- bed_closest(intergenic, genes)

nearby %>%
  select(starts_with("name"), .overlap, .dist) %>%
  filter(abs(.dist) < 1000)
#> # A tibble: 285 x 4
#>         name.x            name.y .overlap .dist
#>          <chr>             <chr>    <int> <int>
#>  1   rs2261631             P704P        0  -267
#>  2 rs570770556             POTEH        0  -912
#>  3 rs538163832             POTEH        0  -952
#>  4   rs9606135            TPTEP1        0  -421
#>  5  rs11912392 ANKRD62P1-PARP4P3        0   104
#>  6   rs8136454          BC038197        0   355
#>  7   rs5992556              XKR3        0  -455
#>  8 rs114101676              GAB4        0   473
#>  9  rs62236167             CECR7        0   261
#> 10   rs5747023             CECR1        0  -386
#> # ... with 275 more rows
```

### Figure 1 code demo


```r
x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 25,     50,
  "chr1", 100,    125
)

y <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 30,     75
)

bed_glyph(bed_intersect(x, y))
```

<img src="figs_files/figure-html/unnamed-chunk-5-1.png" width="672" />


```r
x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1",      1,      50,
  "chr1",      10,     75,
  "chr1",      100,    120
)

bed_glyph(bed_merge(x))
```

<img src="figs_files/figure-html/unnamed-chunk-6-1.png" width="672" />

### grouping data demo


```r
x <- tibble::tribble(
  ~chrom, ~start, ~end, ~strand,
  "chr1", 1,      100,  "+",
  "chr1", 50,     150,  "+",
  "chr2", 100,    200,  "-"
)

y <- tibble::tribble(
  ~chrom, ~start, ~end, ~strand,
  "chr1", 50,     125,  "+",
  "chr1", 50,     150,  "-",
  "chr2", 50,     150,  "+"
)

# intersect tbls by strand
x <- group_by(x, strand)
y <- group_by(y, strand)

bed_intersect(x, y)
#> # A tibble: 2 x 8
#>   chrom start.x end.x strand.x start.y end.y strand.y .overlap
#>   <chr>   <dbl> <dbl>    <chr>   <dbl> <dbl>    <chr>    <int>
#> 1  chr1       1   100        +      50   125        +       50
#> 2  chr1      50   150        +      50   125        +       75
```

Comparisons between intervals on opposite strands are done using the `flip_strands()` function:


```r
x <- group_by(x, strand)

y <- flip_strands(y)
y <- group_by(y, strand)

bed_intersect(x, y)
#> # A tibble: 3 x 8
#>   chrom start.x end.x strand.x start.y end.y strand.y .overlap
#>   <chr>   <dbl> <dbl>    <chr>   <dbl> <dbl>    <chr>    <int>
#> 1  chr2     100   200        -      50   150        -       50
#> 2  chr1       1   100        +      50   150        +       50
#> 3  chr1      50   150        +      50   150        +      100
```

### Figure 2 code demo

```r
bedfile <- valr_example("genes.hg19.chr22.bed.gz")
genomefile <- valr_example("hg19.chrom.sizes.gz")
bgfile  <- valr_example("hela.h3k4.chip.bg.gz")

genes <- read_bed(bedfile, n_fields = 6)
genome <- read_genome(genomefile)
y <- read_bedgraph(bgfile)
```

Then we generate 1 bp intervals to represent transcription start sites (TSSs). We focus on `+` strand genes, but `-` genes are easily accomodated by filtering them and using `bed_makewindows()` with `reversed` window numbers.


```r
# generate 1 bp TSS intervals, "+" strand only
tss <- genes %>%
  filter(strand == "+") %>%
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
```

Now we use the `.win_id` group with `bed_map()` to caluclate a sum by mapping `y` signals onto the intervals in `x`. These data are regrouped by `.win_id` and a summary with `mean` and `sd` values is calculated.


```r
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
```

Finally, these summary statistics are used to construct a plot that illustrates histone density surrounding TSSs.


```r
library(ggplot2)

x_labels <- seq(-region_size, region_size, by = win_size * 5)
x_breaks <- seq(1, 41, by = 5)

sd_limits <- aes(ymax = win_mean + win_sd, ymin = win_mean - win_sd)

p <- ggplot(res, aes(x = .win_id, y = win_mean)) +
  geom_point(size = 0.25) + geom_pointrange(sd_limits, size = 0.1) + 
  scale_x_continuous(labels = x_labels, breaks = x_breaks) + 
  xlab("Position (bp from TSS)") + ylab("Signal") + 
  theme_classic()
```


### Interval statistics

Estimates of significance for interval overlaps can be obtained by combining `bed_shuffle()`, `bed_random()` and the `sample_` functions from `dplyr` with interval statistics in `valr`.

Here we examine the extent of overlap of repeat classes with exons in the human genome (on `chr22` only, for simplicity) using the jaccard similarity index. `bed_jaccard()` implements the jaccard test to examine the similarity between two sets of genomic intervals. Using `bed_shuffle()` and `replicate()` we generate a `data_frame` containing 100 sets of randomly placed intervals then calculate the jaccard index for each set to generate a null-distribution of jaccard scores. Finally an empirical p-value is then calculated from the null-distribution.


```r
library(tidyverse, warn.conflicts = F)

repeats <- read_bed(valr_example("hg19.rmsk.chr22.bed.gz"), n_fields = 6) 
genome <- read_genome(valr_example("hg19.chrom.sizes.gz"))
genes <- read_bed12(valr_example("hg19.refGene.chr22.bed.gz"))
# convert bed12 to bed with exons
exons <- bed12_to_exons(genes)

# function to repeat interval shuffling
shuffle_intervals <- function(n, .data, genome) {
  replicate(n, bed_shuffle(.data, genome), simplify = FALSE) %>%
    bind_rows(.id = "rep") %>%
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
#>  1     1 <tibble [10,000 x 6]> 0.0003388967
#>  2     2 <tibble [10,000 x 6]> 0.0004965988
#>  3     3 <tibble [10,000 x 6]> 0.0002974843
#>  4     4 <tibble [10,000 x 6]> 0.0006899870
#>  5     5 <tibble [10,000 x 6]> 0.0004678412
#>  6     6 <tibble [10,000 x 6]> 0.0001726937
#>  7     7 <tibble [10,000 x 6]> 0.0004694941
#>  8     8 <tibble [10,000 x 6]> 0.0004660410
#>  9     9 <tibble [10,000 x 6]> 0.0006846643
#> 10    10 <tibble [10,000 x 6]> 0.0002143829
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
