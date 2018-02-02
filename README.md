
<!-- README.md is generated from README.Rmd. Please edit that file -->
hlatools
========

Overview
--------

`hlatools` provides a collection of tools to work with [IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/) data

Installation
------------

Install from GitHub using `devtools`:

``` r
devtools::install_github("gscholf/hlatools")
```

Usage
-----

### Get HLA data

-   `fetch_IMGTHLA()`: Clone the *ANHIG/IMGTHLA* repo on [GitHub](https://github.com/ANHIG/IMGTHLA) to a location specified by the option `hlatools.local_repos`. Defaults to *~/local/db/*.

-   `update_IMGTHLA()`: Send a pull request to *ANHIG/IMGTHLA*.

### Read HLA data

A HLA gene can be read into a `HLAGene` object:

``` r
x <- hlatools::HLAGene("DPB1")
x
#> IPD-IMGT/HLA database <3.31.0>; Locus <HLA-DPB1>
#> An object of class 'HLAAllele' for 'HLA-DPB1'
#> DataFrame with 962 rows and 6 columns
#>              allele_name        g_group     p_group      cwd_status
#>              <character>    <character> <character>     <character>
#> 1   HLA-DPB1*01:01:01:01 DPB1*01:01:01G DPB1*01:01P          Common
#> 2   HLA-DPB1*01:01:01:02 DPB1*01:01:01G DPB1*01:01P Not CWD defined
#> 3   HLA-DPB1*01:01:01:03 DPB1*01:01:01G DPB1*01:01P Not CWD defined
#> 4   HLA-DPB1*01:01:01:04 DPB1*01:01:01G DPB1*01:01P Not CWD defined
#> 5   HLA-DPB1*01:01:02:01 DPB1*01:01:02G DPB1*01:01P          Common
#> ...                  ...            ...         ...             ...
#> 958      HLA-DPB1*689:01           None        None Not CWD defined
#> 959      HLA-DPB1*690:01           None        None Not CWD defined
#> 960     HLA-DPB1*691:01N           None        None Not CWD defined
#> 961      HLA-DPB1*692:01           None        None Not CWD defined
#> 962     HLA-DPB1*693:01N           None        None Not CWD defined
#>        SeqLen                                                  FeatureType
#>     <integer>                                                  <character>
#> 1       11468 UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 2       11468 UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 3       11468 UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 4       11468 UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 5       11469 UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> ...       ...                                                          ...
#> 958       264                                                         Exon
#> 959       264                                                         Exon
#> 960       264                                                         Exon
#> 961       264                                                         Exon
#> 962       263                                                         Exon
```

A number of accessor functions can be used to work with these objects:

``` r
x <- hlatools::HLAGene("DPA1")
hlatools::db_version(x)
#> [1] '3.31.0'
hlatools::locusname(x)
#> [1] "HLA-DPA1"
x1 <- x[hlatools::is_complete(x)][1:10]
hlatools::allele_name(x1)
#>  [1] "HLA-DPA1*01:03:01:01" "HLA-DPA1*01:03:01:02" "HLA-DPA1*01:03:01:03"
#>  [4] "HLA-DPA1*01:03:01:04" "HLA-DPA1*01:03:01:05" "HLA-DPA1*01:03:01:06"
#>  [7] "HLA-DPA1*01:03:01:07" "HLA-DPA1*01:03:01:08" "HLA-DPA1*01:03:01:09"
#> [10] "HLA-DPA1*01:03:01:10"
hlatools::cwd_status(x1)
#>  [1] "Common"          "Not CWD defined" "Not CWD defined"
#>  [4] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [7] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [10] "Not CWD defined"
hlatools::ethnicity(x1)
#>  [1] "Caucasoid:Oriental" "Caucasoid"          "Unknown"           
#>  [4] "Unknown"            "Caucasoid"          "Oriental"          
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
#>   A DNAStringSet instance of length 10
#>      width seq                                         names               
#>  [1]  9775 GGTGGACCTGAAAGAAAGAT...AAAAGGAATTGTTTAAAGTA HLA-DPA1*01:03:01:01
#>  [2]  9775 GGTGGACCTGAAAGAAAGAT...AAAAGGAATTGTTTAAAGTA HLA-DPA1*01:03:01:02
#>  [3]  9775 GGTGGACCTGAAAGAAAGAT...AAAAGGAATTGTTTAAAGTA HLA-DPA1*01:03:01:03
#>  [4]  9775 GGTGGACCTGAAAGAAAGAT...AAAAGGAATTGTTTAAAGTA HLA-DPA1*01:03:01:04
#>  [5]  9757 GGTGGACCTGAAAGAAAGAT...AAAAGGAATTGTTTAAAGTA HLA-DPA1*01:03:01:05
#>  [6]  9479 AAATTCTCCCATCTCTTCCC...AAAAGGAATTGTTTAAAGTA HLA-DPA1*01:03:01:06
#>  [7]  9477 AAATTCTCCCATCTCTTCCC...AAAAGGAATTGTTTAAAGTA HLA-DPA1*01:03:01:07
#>  [8]  9478 AAATTCTCCCATCTCTTCCC...AAAAGGAATTGTTTAAAGTA HLA-DPA1*01:03:01:08
#>  [9]  5265 TTACCCAGCAACAGAGAATG...GATAACTAACTGAGTAGTTA HLA-DPA1*01:03:01:09
#> [10]  5266 TTACCCAGCAACAGAGAATG...GATAACTAACTGAGTAGTTA HLA-DPA1*01:03:01:10
hlatools::exon(x1, exon = 2)
#>   A DNAStringSet instance of length 10
#>      width seq                                         names               
#>  [1]   246 CGGACCATGTGTCAACTTAT...CCACACTCAGGCCACCAACG HLA-DPA1*01:03:01...
#>  [2]   246 CGGACCATGTGTCAACTTAT...CCACACTCAGGCCACCAACG HLA-DPA1*01:03:01...
#>  [3]   246 CGGACCATGTGTCAACTTAT...CCACACTCAGGCCACCAACG HLA-DPA1*01:03:01...
#>  [4]   246 CGGACCATGTGTCAACTTAT...CCACACTCAGGCCACCAACG HLA-DPA1*01:03:01...
#>  [5]   246 CGGACCATGTGTCAACTTAT...CCACACTCAGGCCACCAACG HLA-DPA1*01:03:01...
#>  [6]   246 CGGACCATGTGTCAACTTAT...CCACACTCAGGCCACCAACG HLA-DPA1*01:03:01...
#>  [7]   246 CGGACCATGTGTCAACTTAT...CCACACTCAGGCCACCAACG HLA-DPA1*01:03:01...
#>  [8]   246 CGGACCATGTGTCAACTTAT...CCACACTCAGGCCACCAACG HLA-DPA1*01:03:01...
#>  [9]   246 CGGACCATGTGTCAACTTAT...CCACACTCAGGCCACCAACG HLA-DPA1*01:03:01...
#> [10]   246 CGGACCATGTGTCAACTTAT...CCACACTCAGGCCACCAACG HLA-DPA1*01:03:01...
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

These obects come with an additional set of `R6`-methods that implement an API used mostly in conjunction with the `DR2S` package:

-   `x$get_closest_complete_neighbor(allele)`: Find the closest full-length allele based on the genetic distance between exon 2 sequences.

-   `x$get_reference_sequence(allele)`: Return a `BStringSet` object of a full-length reference sequence for `allele`. If `allele` is not fully known the missing stretches are taken from the allele returned by a call to `get_closest_complete_neighbor()`.

### Session Info

    #> Session info -------------------------------------------------------------
    #>  setting  value                       
    #>  version  R version 3.4.3 (2017-11-30)
    #>  system   x86_64, linux-gnu           
    #>  ui       X11                         
    #>  language en_GB:en                    
    #>  collate  en_GB.UTF-8                 
    #>  tz       Europe/Berlin               
    #>  date     2018-02-02
    #> Packages -----------------------------------------------------------------
    #>  package          * version   date       source                     
    #>  assertive.base     0.0-7     2016-12-30 CRAN (R 3.4.2)             
    #>  backports          1.1.2     2017-12-13 cran (@1.1.2)              
    #>  base             * 3.4.3     2017-12-01 local                      
    #>  BiocGenerics       0.24.0    2017-11-17 Bioconductor               
    #>  Biostrings         2.46.0    2017-11-17 Bioconductor               
    #>  bitops             1.0-6     2013-08-17 CRAN (R 3.4.2)             
    #>  codetools          0.2-15    2016-10-05 CRAN (R 3.3.1)             
    #>  compiler           3.4.3     2017-12-01 local                      
    #>  data.table         1.10.4-3  2017-10-27 cran (@1.10.4-)            
    #>  datasets         * 3.4.3     2017-12-01 local                      
    #>  devtools           1.13.4    2017-11-09 cran (@1.13.4)             
    #>  digest             0.6.14    2018-01-14 cran (@0.6.14)             
    #>  evaluate           0.10.1    2017-06-24 cran (@0.10.1)             
    #>  foreach            1.4.4     2017-12-12 cran (@1.4.4)              
    #>  GenomeInfoDb       1.14.0    2017-11-17 Bioconductor               
    #>  GenomeInfoDbData   0.99.1    2017-11-17 Bioconductor               
    #>  GenomicRanges      1.30.0    2017-11-17 Bioconductor               
    #>  git2r              0.21.0    2018-01-04 cran (@0.21.0)             
    #>  graphics         * 3.4.3     2017-12-01 local                      
    #>  grDevices        * 3.4.3     2017-12-01 local                      
    #>  hlatools         * 0.0.7     2018-01-25 local (gschofl/hlatools@NA)
    #>  htmltools          0.3.6     2017-04-28 CRAN (R 3.4.2)             
    #>  IRanges            2.12.0    2017-11-17 Bioconductor               
    #>  iterators          1.0.9     2017-12-12 cran (@1.0.9)              
    #>  knitr              1.18      2017-12-27 cran (@1.18)               
    #>  magrittr           1.5       2014-11-22 CRAN (R 3.4.2)             
    #>  memoise            1.1.0     2017-04-21 CRAN (R 3.4.2)             
    #>  methods          * 3.4.3     2017-12-01 local                      
    #>  parallel           3.4.3     2017-12-01 local                      
    #>  R6                 2.2.2     2017-06-17 cran (@2.2.2)              
    #>  Rcpp               0.12.15   2018-01-20 cran (@0.12.15)            
    #>  RCurl              1.95-4.10 2018-01-04 cran (@1.95-4.)            
    #>  rmarkdown          1.8       2017-11-17 cran (@1.8)                
    #>  rprojroot          1.3-2     2018-01-03 cran (@1.3-2)              
    #>  S4Vectors          0.16.0    2017-11-17 Bioconductor               
    #>  stats            * 3.4.3     2017-12-01 local                      
    #>  stats4             3.4.3     2017-12-01 local                      
    #>  stringi            1.1.6     2017-11-17 cran (@1.1.6)              
    #>  stringr            1.2.0     2017-02-18 CRAN (R 3.4.2)             
    #>  tools              3.4.3     2017-12-01 local                      
    #>  utils            * 3.4.3     2017-12-01 local                      
    #>  withr              2.1.1     2017-12-19 cran (@2.1.1)              
    #>  xml2               1.2.0     2018-01-24 cran (@1.2.0)              
    #>  XVector            0.18.0    2017-11-17 Bioconductor               
    #>  yaml               2.1.16    2017-12-12 cran (@2.1.16)             
    #>  zlibbioc           1.24.0    2017-11-17 Bioconductor
