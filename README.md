
<!-- README.md is generated from README.Rmd. Please edit that file -->
hlatools
========

A collection of tools to work with [IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/) data

Main functionality
------------------

### Fetch HLA data

A number of functions exist that grab HLA-related data from the internet

-   `fetch_IMGTHLA`: Clone the *ANHIG/IMGTHLA* repo on [GitHub](https://github.com/ANHIG/IMGTHLA) to a location specified by the option `hlatools.local_repos`. Defaults to *~/local/db/*.

-   `update_IMGTHLA`: Send a pull request to *ANHIG/IMGTHLA*.

-   `read_hla_xml(remote = FALSE)`: Read the *hla.xml* file either fetching it from the IPD-IMGT/HLA ftp server or the local *ANHIG/IMGTHLA* clone.

### Parse HLA data

A HLA locus can be read into a `HLAAllele` object:

``` r
doc  <- hlatools::read_hla_xml()
dpb1 <- hlatools::parse_hla_alleles(doc, "DPB1")
dpb1
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
dpb1_cmpl <- dpb1[hlatools::is_complete(dpb1)][1:10]
hlatools::allele_name(dpb1_cmpl)
#>  [1] "HLA-DPB1*01:01:01:01" "HLA-DPB1*01:01:01:02" "HLA-DPB1*01:01:01:03"
#>  [4] "HLA-DPB1*01:01:01:04" "HLA-DPB1*01:01:02:01" "HLA-DPB1*01:01:02:02"
#>  [7] "HLA-DPB1*02:01:02:01" "HLA-DPB1*02:01:02:02" "HLA-DPB1*02:01:02:03"
#> [10] "HLA-DPB1*02:01:02:04"
hlatools::cwd_status(dpb1_cmpl)
#>  [1] "Common"          "Not CWD defined" "Not CWD defined"
#>  [4] "Not CWD defined" "Common"          "Not CWD defined"
#>  [7] "Common"          "Not CWD defined" "Not CWD defined"
#> [10] "Not CWD defined"
hlatools::sequences(dpb1_cmpl)
#>   A DNAStringSet instance of length 10
#>      width seq                                         names               
#>  [1] 11468 TAATCCCTGTAGATGGGCCA...TATAATCTAATACACTTTAA HLA-DPB1*01:01:01:01
#>  [2] 11468 TAATCCCTGTAGATGGGCCA...TATAATCTAATACACTTTAA HLA-DPB1*01:01:01:02
#>  [3] 11468 TAATCCCTGTAGATGGGCCA...TATAATCTAATACACTTTAA HLA-DPB1*01:01:01:03
#>  [4] 11468 TAATCCCTGTAGATGGGCCA...TATAATCTAATACACTTTAA HLA-DPB1*01:01:01:04
#>  [5] 11469 TAATCCCTGTAGATGGGCCA...TATAATCTAATACACTTTAA HLA-DPB1*01:01:02:01
#>  [6] 11469 TAATCCCTGTAGATGGGCCA...TATAATCTAATACACTTTAA HLA-DPB1*01:01:02:02
#>  [7] 11532 TAATCCCTGTAGATGGGCCA...TATAATCTAATACACTTTAA HLA-DPB1*02:01:02:01
#>  [8] 11525 TAATCCCTGTAGATGGGCCA...TATAATCTAATACACTTTAA HLA-DPB1*02:01:02:02
#>  [9] 11516 TAATCCCTGTAGATGGGCCA...TATAATCTAATACACTTTAA HLA-DPB1*02:01:02:03
#> [10] 11528 TAATCCCTGTAGATGGGCCA...TATAATCTAATACACTTTAA HLA-DPB1*02:01:02:04
```

At a slightly higher level a `HLAAllele` can be encapsulated in an `R6`-based `HLAGene` object:

``` r
x <- hlatools::HLAGene("DPB1", db_version = "3.30.0")
x
#> IPD-IMGT/HLA database <3.30.0>; Locus <HLA-DPB1>
#> An object of class 'HLAAllele' for 'HLA-DPB1'
#> DataFrame with 942 rows and 6 columns
#>              allele_name        g_group     p_group      cwd_status
#>              <character>    <character> <character>     <character>
#> 1   HLA-DPB1*01:01:01:01 DPB1*01:01:01G DPB1*01:01P          Common
#> 2   HLA-DPB1*01:01:01:02 DPB1*01:01:01G DPB1*01:01P Not CWD defined
#> 3   HLA-DPB1*01:01:01:03 DPB1*01:01:01G DPB1*01:01P Not CWD defined
#> 4   HLA-DPB1*01:01:01:04 DPB1*01:01:01G DPB1*01:01P Not CWD defined
#> 5   HLA-DPB1*01:01:02:01 DPB1*01:01:02G DPB1*01:01P          Common
#> ...                  ...            ...         ...             ...
#> 938      HLA-DPB1*674:01           None DPB1*04:02P Not CWD defined
#> 939      HLA-DPB1*675:01 DPB1*03:01:01G DPB1*03:01P Not CWD defined
#> 940      HLA-DPB1*676:01 DPB1*03:01:01G DPB1*03:01P Not CWD defined
#> 941      HLA-DPB1*677:01           None DPB1*04:01P Not CWD defined
#> 942      HLA-DPB1*678:01 DPB1*02:01:02G DPB1*02:01P Not CWD defined
#>        SeqLen                                                  FeatureType
#>     <integer>                                                  <character>
#> 1       11468 UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 2       11468 UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 3       11468 UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 4       11468 UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 5       11469 UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> ...       ...                                                          ...
#> 938     11516 UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 939     11461 UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 940     11461 UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 941     11526 UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 942     11532 UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR

hlatools::db_version(x)
#> [1] '3.30.0'
hlatools::locusname(x)
#> [1] "HLA-DPB1"
x[hlatools::cwd_status(x) == "Common"][1:10]
#> An object of class 'HLAAllele' for 'HLA-DPB1'
#> DataFrame with 10 rows and 6 columns
#>             allele_name        g_group     p_group  cwd_status    SeqLen
#>             <character>    <character> <character> <character> <integer>
#> 1  HLA-DPB1*01:01:01:01 DPB1*01:01:01G DPB1*01:01P      Common     11468
#> 2  HLA-DPB1*01:01:02:01 DPB1*01:01:02G DPB1*01:01P      Common     11469
#> 3  HLA-DPB1*02:01:02:01 DPB1*02:01:02G DPB1*02:01P      Common     11532
#> 4  HLA-DPB1*02:02:01:01 DPB1*02:02:01G DPB1*02:02P      Common     11528
#> 5  HLA-DPB1*03:01:01:01 DPB1*03:01:01G DPB1*03:01P      Common     11461
#> 6  HLA-DPB1*04:01:01:01 DPB1*04:01:01G DPB1*04:01P      Common     11526
#> 7  HLA-DPB1*04:02:01:01 DPB1*04:02:01G DPB1*04:02P      Common     11516
#> 8  HLA-DPB1*05:01:01:01 DPB1*05:01:01G DPB1*05:01P      Common     11466
#> 9  HLA-DPB1*06:01:01:01 DPB1*06:01:01G DPB1*06:01P      Common     11461
#> 10    HLA-DPB1*09:01:01           None DPB1*09:01P      Common     11466
#>                                                     FeatureType
#>                                                     <character>
#> 1  UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 2  UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 3  UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 4  UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 5  UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 6  UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 7  UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 8  UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 9  UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 10 UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
```

These obects come with the same set of getter functions as `HLAAllele`s plus an additional set of `R6`-methods that implement an API used mostly in conjunction with the `DR2S` package:

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
    #>  date     2018-01-25
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
