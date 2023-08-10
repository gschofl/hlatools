
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hlatools

## Overview

`hlatools` provides a collection of tools to work with
[IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/) data

## Installation

Install from GitHub using `remotes`:

``` r
remotes::install_github("gschofl/hlatools")
```

## Usage

### Get HLA data

  - `clone_IMGTHLA()`: Clone the *ANHIG/IMGTHLA* repo on
    [GitHub](https://github.com/ANHIG/IMGTHLA) to a location specified
    by the option `hlatools.local_repos`. Defaults to *\~/local/db/*.

  - `pull_IMGTHLA()`: Send a pull request to *ANHIG/IMGTHLA*.

  - `check_IMGTHLA()`: Check if your local repository is up-to-date.

### Read HLA data

A HLA gene can be read into a `HLAGene` object:

``` r
x <- hlatools::HLAGene("DPB1")
x
#> IPD-IMGT/HLA database <3.41.0>; Locus <HLA-DPB1>
#> An object of class 'HLAAllele' for 'HLA-DPB1'
#> DataFrame with 1584 rows and 6 columns
#>               allele_name        g_group     p_group      cwd_status    SeqLen
#>               <character>    <character> <character>     <character> <integer>
#> 1    HLA-DPB1*01:01:01:01 DPB1*01:01:01G DPB1*01:01P          Common     11468
#> 2    HLA-DPB1*01:01:01:02 DPB1*01:01:01G DPB1*01:01P Not CWD defined     11468
#> 3    HLA-DPB1*01:01:01:03 DPB1*01:01:01G DPB1*01:01P Not CWD defined     11468
#> 4    HLA-DPB1*01:01:01:04 DPB1*01:01:01G DPB1*01:01P Not CWD defined     11468
#> 5    HLA-DPB1*01:01:01:05 DPB1*01:01:01G DPB1*01:01P Not CWD defined     11030
#> ...                   ...            ...         ...             ...       ...
#> 1580     HLA-DPB1*1116:01           None        None Not CWD defined       677
#> 1581     HLA-DPB1*1117:01 DPB1*05:01:01G DPB1*05:01P Not CWD defined       677
#> 1582     HLA-DPB1*1118:01 DPB1*05:01:01G DPB1*05:01P Not CWD defined       657
#> 1583     HLA-DPB1*1119:01 DPB1*05:01:01G DPB1*05:01P Not CWD defined       657
#> 1584     HLA-DPB1*1120:01 DPB1*05:01:01G DPB1*05:01P Not CWD defined       657
#>                                                       FeatureType
#>                                                       <character>
#> 1    UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 2    UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 3    UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 4    UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 5    UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> ...                                                           ...
#> 1580                                          Exon:Exon:Exon:Exon
#> 1581                                          Exon:Exon:Exon:Exon
#> 1582                                               Exon:Exon:Exon
#> 1583                                               Exon:Exon:Exon
#> 1584                                               Exon:Exon:Exon
```

A number of accessor functions can be used to work with these objects:

``` r
x <- hlatools::HLAGene("DPA1")
hlatools::db_version(x)
#> [1] '3.41.0'
hlatools::locusname(x)
#> [1] "HLA-DPA1"
x1 <- x[hlatools::is_complete(x)][1:10]
hlatools::allele_name(x1)
#>  [1] "HLA-DPA1*01:03:01:01" "HLA-DPA1*01:03:01:02" "HLA-DPA1*01:03:01:03"
#>  [4] "HLA-DPA1*01:03:01:04" "HLA-DPA1*01:03:01:05" "HLA-DPA1*01:03:01:06"
#>  [7] "HLA-DPA1*01:03:01:07" "HLA-DPA1*01:03:01:08" "HLA-DPA1*01:03:01:09"
#> [10] "HLA-DPA1*01:03:01:10"
hlatools::cwd_status(x1)
#>  [1] "Common"          "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [5] "Not CWD defined" "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [9] "Not CWD defined" "Not CWD defined"
hlatools::ancestry(x1)
#>  [1] "Caucasoid:Oriental" "Caucasoid:Black"    "Caucasoid"         
#>  [4] "Caucasoid"          "Caucasoid"          "Oriental"          
#>  [7] "Oriental"           "Oriental"           "Caucasoid"         
#> [10] "Caucasoid"
hlatools::g_group(x1)
#>  [1] "DPA1*01:03:01G" "DPA1*01:03:01G" "DPA1*01:03:01G" "DPA1*01:03:01G"
#>  [5] "DPA1*01:03:01G" "DPA1*01:03:01G" "DPA1*01:03:01G" "DPA1*01:03:01G"
#>  [9] "DPA1*01:03:01G" "DPA1*01:03:01G"
hlatools::p_group(x1)
#>  [1] "DPA1*01:03P" "DPA1*01:03P" "DPA1*01:03P" "DPA1*01:03P" "DPA1*01:03P"
#>  [6] "DPA1*01:03P" "DPA1*01:03P" "DPA1*01:03P" "DPA1*01:03P" "DPA1*01:03P"
hlatools::sequences(x1)
#> DNAStringSet object of length 10:
#>      width seq                                              names               
#>  [1]  9775 GGTGGACCTGAAAGAAAGATTAA...AAAAAAGGAATTGTTTAAAGTA HLA-DPA1*01:03:01:01
#>  [2]  9775 GGTGGACCTGAAAGAAAGATTAA...AAAAAAGGAATTGTTTAAAGTA HLA-DPA1*01:03:01:02
#>  [3]  9775 GGTGGACCTGAAAGAAAGATTAA...AAAAAAGGAATTGTTTAAAGTA HLA-DPA1*01:03:01:03
#>  [4]  9775 GGTGGACCTGAAAGAAAGATTAA...AAAAAAGGAATTGTTTAAAGTA HLA-DPA1*01:03:01:04
#>  [5]  9757 GGTGGACCTGAAAGAAAGATTAA...AAAAAAGGAATTGTTTAAAGTA HLA-DPA1*01:03:01:05
#>  [6]  9479 AAATTCTCCCATCTCTTCCCCAG...AAAAAAGGAATTGTTTAAAGTA HLA-DPA1*01:03:01:06
#>  [7]  9477 AAATTCTCCCATCTCTTCCCCAG...AAAAAAGGAATTGTTTAAAGTA HLA-DPA1*01:03:01:07
#>  [8]  9478 AAATTCTCCCATCTCTTCCCCAG...AAAAAAGGAATTGTTTAAAGTA HLA-DPA1*01:03:01:08
#>  [9]  5265 TTACCCAGCAACAGAGAATGTCA...CCGATAACTAACTGAGTAGTTA HLA-DPA1*01:03:01:09
#> [10]  5735 TTTTTACATCTCTTTCTCTAACT...GGTGGGTGCCTGTAACTACTTA HLA-DPA1*01:03:01:10
hlatools::exon(x1, exon = 2)
#> DNAStringSet object of length 10:
#>      width seq                                              names               
#>  [1]   246 CGGACCATGTGTCAACTTATGCC...AACCACACTCAGGCCACCAACG HLA-DPA1*01:03:01...
#>  [2]   246 CGGACCATGTGTCAACTTATGCC...AACCACACTCAGGCCACCAACG HLA-DPA1*01:03:01...
#>  [3]   246 CGGACCATGTGTCAACTTATGCC...AACCACACTCAGGCCACCAACG HLA-DPA1*01:03:01...
#>  [4]   246 CGGACCATGTGTCAACTTATGCC...AACCACACTCAGGCCACCAACG HLA-DPA1*01:03:01...
#>  [5]   246 CGGACCATGTGTCAACTTATGCC...AACCACACTCAGGCCACCAACG HLA-DPA1*01:03:01...
#>  [6]   246 CGGACCATGTGTCAACTTATGCC...AACCACACTCAGGCCACCAACG HLA-DPA1*01:03:01...
#>  [7]   246 CGGACCATGTGTCAACTTATGCC...AACCACACTCAGGCCACCAACG HLA-DPA1*01:03:01...
#>  [8]   246 CGGACCATGTGTCAACTTATGCC...AACCACACTCAGGCCACCAACG HLA-DPA1*01:03:01...
#>  [9]   246 CGGACCATGTGTCAACTTATGCC...AACCACACTCAGGCCACCAACG HLA-DPA1*01:03:01...
#> [10]   246 CGGACCATGTGTCAACTTATGCC...AACCACACTCAGGCCACCAACG HLA-DPA1*01:03:01...
```

We can access previous releases of the IPD-IMGT/HLA database:

``` r
x <- hlatools::HLAGene("DPB1", db_version = "3.25.0")
x
#> IPD-IMGT/HLA database <3.25.0>; Locus <HLA-DPB1>
#> An object of class 'HLAAllele' for 'HLA-DPB1'
#> DataFrame with 671 rows and 6 columns
#>           allele_name     g_group     p_group      cwd_status    SeqLen
#>           <character> <character> <character>     <character> <integer>
#> 1   HLA-DPB1*01:01:01                                  Common       777
#> 2   HLA-DPB1*01:01:02                                  Common       777
#> 3   HLA-DPB1*01:01:03                         Not CWD defined       264
#> 4   HLA-DPB1*01:01:04                         Not CWD defined       264
#> 5   HLA-DPB1*01:01:05                         Not CWD defined       264
#> ...               ...         ...         ...             ...       ...
#> 667   HLA-DPB1*568:01                         Not CWD defined       264
#> 668   HLA-DPB1*569:01                         Not CWD defined       264
#> 669  HLA-DPB1*570:01N                         Not CWD defined       264
#> 670   HLA-DPB1*571:01                         Not CWD defined       657
#> 671   HLA-DPB1*572:01                         Not CWD defined       546
#>                  FeatureType
#>                  <character>
#> 1   Exon:Exon:Exon:Exon:Exon
#> 2   Exon:Exon:Exon:Exon:Exon
#> 3                       Exon
#> 4                       Exon
#> 5                       Exon
#> ...                      ...
#> 667                     Exon
#> 668                     Exon
#> 669                     Exon
#> 670           Exon:Exon:Exon
#> 671                Exon:Exon
hlatools::db_version(x)
#> [1] '3.25.0'
hlatools::locusname(x)
#> [1] "HLA-DPB1"
x[hlatools::cwd_status(x) == "Common"][1:10]
#> An object of class 'HLAAllele' for 'HLA-DPB1'
#> DataFrame with 10 rows and 6 columns
#>             allele_name     g_group     p_group  cwd_status    SeqLen
#>             <character> <character> <character> <character> <integer>
#> 1     HLA-DPB1*01:01:01                              Common       777
#> 2     HLA-DPB1*01:01:02                              Common       777
#> 3     HLA-DPB1*02:01:02                              Common     11532
#> 4        HLA-DPB1*02:02                              Common     11528
#> 5     HLA-DPB1*03:01:01                              Common     11461
#> 6  HLA-DPB1*04:01:01:01                              Common     11526
#> 7  HLA-DPB1*04:02:01:01                              Common     11516
#> 8     HLA-DPB1*05:01:01                              Common       777
#> 9     HLA-DPB1*06:01:01                              Common       777
#> 10    HLA-DPB1*09:01:01                              Common       777
#>                                                     FeatureType
#>                                                     <character>
#> 1                                      Exon:Exon:Exon:Exon:Exon
#> 2                                      Exon:Exon:Exon:Exon:Exon
#> 3  UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 4  UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 5  UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 6  UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 7  UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 8                                      Exon:Exon:Exon:Exon:Exon
#> 9                                      Exon:Exon:Exon:Exon:Exon
#> 10                                     Exon:Exon:Exon:Exon:Exon
```

These objects come with an additional set of `R6`-methods that implement
an API used mostly in conjunction with the `DR2S` package:

  - `x$get_closest_complete_neighbor(allele, partially = TRUE)`: Find
    the closest full-length allele based on the genetic distance between
    exon sequences.

  - `x$get_reference_sequence(allele)`: Return a `BStringSet` object of
    a full-length reference sequence for `allele`. If `allele` is not
    fully known the missing stretches are taken from the allele returned
    by a call to `get_closest_complete_neighbor()`.

### Session Info

    #> ─ Session info ───────────────────────────────────────────────────────────────
    #>  setting  value                       
    #>  version  R version 4.0.2 (2020-06-22)
    #>  os       Ubuntu 20.04.1 LTS          
    #>  system   x86_64, linux-gnu           
    #>  ui       X11                         
    #>  language en_GB:en                    
    #>  collate  en_GB.UTF-8                 
    #>  ctype    en_GB.UTF-8                 
    #>  tz       Europe/Berlin               
    #>  date     2020-07-31                  
    #> 
    #> ─ Packages ───────────────────────────────────────────────────────────────────
    #>  package              * version  date       lib source        
    #>  assertive.base         0.0-7    2016-12-30 [1] CRAN (R 4.0.0)
    #>  assertive.properties   0.0-4    2016-12-30 [1] CRAN (R 4.0.0)
    #>  assertthat             0.2.1    2019-03-21 [1] CRAN (R 4.0.0)
    #>  backports              1.1.8    2020-06-17 [1] CRAN (R 4.0.1)
    #>  BiocGenerics           0.34.0   2020-04-27 [1] Bioconductor  
    #>  Biostrings             2.56.0   2020-04-27 [1] Bioconductor  
    #>  bitops                 1.0-6    2013-08-17 [1] CRAN (R 4.0.0)
    #>  callr                  3.4.3    2020-03-28 [1] CRAN (R 4.0.0)
    #>  cli                    2.0.2    2020-02-28 [1] CRAN (R 4.0.0)
    #>  codetools              0.2-16   2018-12-24 [1] CRAN (R 4.0.0)
    #>  crayon                 1.3.4    2017-09-16 [1] CRAN (R 4.0.0)
    #>  data.table             1.13.0   2020-07-24 [1] CRAN (R 4.0.2)
    #>  desc                   1.2.0    2018-05-01 [1] CRAN (R 4.0.0)
    #>  devtools               2.3.1    2020-07-21 [1] CRAN (R 4.0.2)
    #>  digest                 0.6.25   2020-02-23 [1] CRAN (R 4.0.0)
    #>  dplyr                  1.0.1    2020-07-31 [1] CRAN (R 4.0.2)
    #>  ellipsis               0.3.1    2020-05-15 [1] CRAN (R 4.0.0)
    #>  evaluate               0.14     2019-05-28 [1] CRAN (R 4.0.0)
    #>  fansi                  0.4.1    2020-01-08 [1] CRAN (R 4.0.0)
    #>  foreach                1.5.0    2020-03-30 [1] CRAN (R 4.0.0)
    #>  fs                     1.4.2    2020-06-30 [1] CRAN (R 4.0.2)
    #>  generics               0.0.2    2018-11-29 [1] CRAN (R 4.0.0)
    #>  GenomeInfoDb           1.24.2   2020-06-15 [1] Bioconductor  
    #>  GenomeInfoDbData       1.2.3    2020-05-05 [1] Bioconductor  
    #>  GenomicRanges          1.40.0   2020-04-27 [1] Bioconductor  
    #>  git2r                  0.27.1   2020-05-03 [1] CRAN (R 4.0.0)
    #>  glue                   1.4.1    2020-05-13 [1] CRAN (R 4.0.0)
    #>  hlatools             * 0.1.5    2020-07-31 [1] local         
    #>  htmltools              0.5.0    2020-06-16 [1] CRAN (R 4.0.1)
    #>  IRanges                2.22.2   2020-05-21 [1] Bioconductor  
    #>  iterators              1.0.12   2019-07-26 [1] CRAN (R 4.0.0)
    #>  knitr                  1.29     2020-06-23 [1] CRAN (R 4.0.2)
    #>  lifecycle              0.2.0    2020-03-06 [1] CRAN (R 4.0.0)
    #>  magrittr               1.5      2014-11-22 [1] CRAN (R 4.0.0)
    #>  memoise                1.1.0    2017-04-21 [1] CRAN (R 4.0.0)
    #>  pillar                 1.4.6    2020-07-10 [1] CRAN (R 4.0.2)
    #>  pkgbuild               1.1.0    2020-07-13 [1] CRAN (R 4.0.2)
    #>  pkgconfig              2.0.3    2019-09-22 [1] CRAN (R 4.0.0)
    #>  pkgload                1.1.0    2020-05-29 [1] CRAN (R 4.0.1)
    #>  prettyunits            1.1.1    2020-01-24 [1] CRAN (R 4.0.0)
    #>  processx               3.4.3    2020-07-05 [1] CRAN (R 4.0.2)
    #>  ps                     1.3.3    2020-05-08 [1] CRAN (R 4.0.0)
    #>  purrr                  0.3.4    2020-04-17 [1] CRAN (R 4.0.0)
    #>  R6                     2.4.1    2019-11-12 [1] CRAN (R 4.0.0)
    #>  Rcpp                   1.0.5    2020-07-06 [1] CRAN (R 4.0.2)
    #>  RCurl                  1.98-1.2 2020-04-18 [1] CRAN (R 4.0.0)
    #>  remotes                2.2.0    2020-07-21 [1] CRAN (R 4.0.2)
    #>  rlang                  0.4.7    2020-07-09 [1] CRAN (R 4.0.2)
    #>  rmarkdown              2.3      2020-06-18 [1] CRAN (R 4.0.2)
    #>  rprojroot              1.3-2    2018-01-03 [1] CRAN (R 4.0.0)
    #>  S4Vectors              0.26.1   2020-05-16 [1] Bioconductor  
    #>  sessioninfo            1.1.1    2018-11-05 [1] CRAN (R 4.0.0)
    #>  stringi                1.4.6    2020-02-17 [1] CRAN (R 4.0.0)
    #>  stringr                1.4.0    2019-02-10 [1] CRAN (R 4.0.0)
    #>  testthat               2.3.2    2020-03-02 [1] CRAN (R 4.0.0)
    #>  tibble                 3.0.3    2020-07-10 [1] CRAN (R 4.0.2)
    #>  tidyselect             1.1.0    2020-05-11 [1] CRAN (R 4.0.0)
    #>  usethis                1.6.1    2020-04-29 [1] CRAN (R 4.0.0)
    #>  vctrs                  0.3.2    2020-07-15 [1] CRAN (R 4.0.2)
    #>  withr                  2.2.0    2020-04-20 [1] CRAN (R 4.0.0)
    #>  xfun                   0.16     2020-07-24 [1] CRAN (R 4.0.2)
    #>  xml2                   1.3.2    2020-04-23 [1] CRAN (R 4.0.0)
    #>  XVector                0.28.0   2020-04-27 [1] Bioconductor  
    #>  yaml                   2.2.1    2020-02-01 [1] CRAN (R 4.0.0)
    #>  zlibbioc               1.34.0   2020-04-27 [1] Bioconductor  
    #> 
    #> [1] /home/gerhard/local/R/Library
    #> [2] /usr/local/lib/R/site-library
    #> [3] /usr/lib/R/site-library
    #> [4] /usr/lib/R/library
