
<!-- README.md is generated from README.Rmd. Please edit that file -->
hlatools
========

A collection of tools to work with IPD-IMGT/HLA data

Main functionality
------------------

### Fetch HLA data

A number of functions exist that grab HLA-related data from the internet

-   `fetch_IMGTHLA`: Clone the *ANHIG/IMGTHLA* repo on GitHub to a location specified by the option `hlatools.local_repos`. Defaults to *~/local/db/*.

-   `update_IMGTHLS`: Send a pull request to *ANHIG/IMGTHLA*.

-   `read_hla_xml(remote = FALSE)`: Read the *hla.xml* file either fetching it from the IMGT/HLA ftp server or the local *ANHIG/IMGTHLA* clone.

### Parse HLA data

A HLA locus can be read into a `HLAAllele` object:

``` r
doc  <- hlatools::read_hla_xml()
dpb1 <- hlatools::parse_hla_alleles(doc, "DPB1")
dpb1
#> An object of class 'HLAAllele'
#> DataFrame with 472 rows and 6 columns
#>           allele_name     g_group     p_group  cwd_status    SeqLen
#>           <character> <character> <character> <character> <integer>
#> 1   HLA-DPB1*01:01:01                                           777
#> 2   HLA-DPB1*01:01:02                                           777
#> 3   HLA-DPB1*01:01:03                                           267
#> 4   HLA-DPB1*01:01:04                                           264
#> 5   HLA-DPB1*01:01:05                                           264
#> ...               ...         ...         ...         ...       ...
#> 468  HLA-DPB1*401:01N                                           264
#> 469   HLA-DPB1*402:01                                           264
#> 470  HLA-DPB1*403:01N                                           264
#> 471   HLA-DPB1*404:01                                           264
#> 472   HLA-DPB1*405:01                                           264
#>                  FeatureType
#>                  <character>
#> 1   Exon:Exon:Exon:Exon:Exon
#> 2   Exon:Exon:Exon:Exon:Exon
#> 3                  Exon:Exon
#> 4                       Exon
#> 5                       Exon
#> ...                      ...
#> 468                     Exon
#> 469                     Exon
#> 470                     Exon
#> 471                     Exon
#> 472                     Exon
```

A number of accessor functions can be used to work with these objects:

``` r
dpb1_cmpl <- dpb1[hlatools::is_complete(dpb1)]
hlatools::allele_name(dpb1_cmpl)
#> [1] "HLA-DPB1*02:01:02"    "HLA-DPB1*02:02"       "HLA-DPB1*03:01:01"   
#> [4] "HLA-DPB1*04:01:01:01" "HLA-DPB1*04:01:01:02" "HLA-DPB1*04:02:01:01"
#> [7] "HLA-DPB1*04:02:01:02"
hlatools::cwd_status(dpb1_cmpl)
#> [1] "" "" "" "" "" "" ""
hlatools::sequences(dpb1_cmpl)
#>   A DNAStringSet instance of length 7
#>     width seq                                          names               
#> [1] 11532 TAATCCCTGTAGATGGGCCAG...TATAATCTAATACACTTTAA HLA-DPB1*02:01:02
#> [2] 11528 TAATCCCTGTAGATGGGCCAG...TATAATCTAATACACTTTAA HLA-DPB1*02:02
#> [3] 11461 TAATCCCTGTAGATGGGCCAG...TATAATCTAATACACTTTAA HLA-DPB1*03:01:01
#> [4] 11526 TAATCCCTGTAGATGGGCCAG...TATAATCTAATACACTTTAA HLA-DPB1*04:01:01:01
#> [5] 11518 TAATCCCTGTAGATGGGCCAG...TATAATCTAATACACTTTAA HLA-DPB1*04:01:01:02
#> [6] 11516 TAATCCCTGTAGATGGGCCAG...TATAATCTAATACACTTTAA HLA-DPB1*04:02:01:01
#> [7] 11516 TAATCCCTGTAGATGGGCCAG...TATAATCTAATACACTTTAA HLA-DPB1*04:02:01:02
```

At a slightly higher level a `HLAAllele` can be encapsulated in an `R6`-based `HLAGene` object:

``` r
x <- hlatools::HLAGene("DPB1", db_version = "3.28.0")
x
#> IMGT/HLA database <3.28.0>; Locus <HLA-DPB1>
#> An object of class 'HLAAllele'
#> DataFrame with 828 rows and 6 columns
#>               allele_name        g_group     p_group      cwd_status
#>               <character>    <character> <character>     <character>
#> 1    HLA-DPB1*01:01:01:01 DPB1*01:01:01G DPB1*01:01P          Common
#> 2    HLA-DPB1*01:01:01:02 DPB1*01:01:01G DPB1*01:01P Not CWD defined
#> 3    HLA-DPB1*01:01:01:03 DPB1*01:01:01G DPB1*01:01P Not CWD defined
#> 4    HLA-DPB1*01:01:01:04 DPB1*01:01:01G DPB1*01:01P Not CWD defined
#> 5    HLA-DPB1*01:01:02:01 DPB1*01:01:02G DPB1*01:01P          Common
#> ...                   ...            ...         ...             ...
#> 824       HLA-DPB1*645:01           None        None Not CWD defined
#> 825       HLA-DPB1*646:01           None        None Not CWD defined
#> 826       HLA-DPB1*647:01 DPB1*04:02:01G DPB1*04:02P Not CWD defined
#> 827 HLA-DPB1*648:01:01:01 DPB1*57:01:01G DPB1*57:01P Not CWD defined
#> 828 HLA-DPB1*648:01:01:02 DPB1*57:01:01G DPB1*57:01P Not CWD defined
#>        SeqLen                                                  FeatureType
#>     <integer>                                                  <character>
#> 1       11468 UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 2       11468 UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 3       11468 UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 4       11468 UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 5       11469 UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> ...       ...                                                          ...
#> 824       264                                                         Exon
#> 825       264                                                         Exon
#> 826       677                                          Exon:Exon:Exon:Exon
#> 827     11466 UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 828     11468 UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR

hlatools::db_version(x)
#> [1] "3.28.0"
hlatools::locusname(x)
#> [1] "HLA-DPB1"
x[hlatools::cwd_status(x) == "Common"]
#> An object of class 'HLAAllele'
#> DataFrame with 40 rows and 6 columns
#>               allele_name        g_group     p_group  cwd_status    SeqLen
#>               <character>    <character> <character> <character> <integer>
#> 1    HLA-DPB1*01:01:01:01 DPB1*01:01:01G DPB1*01:01P      Common     11468
#> 2    HLA-DPB1*01:01:02:01 DPB1*01:01:02G DPB1*01:01P      Common     11469
#> 3    HLA-DPB1*02:01:02:01 DPB1*02:01:02G DPB1*02:01P      Common     11532
#> 4    HLA-DPB1*02:02:01:01 DPB1*02:02:01G DPB1*02:02P      Common     11528
#> 5    HLA-DPB1*03:01:01:01 DPB1*03:01:01G DPB1*03:01P      Common     11461
#> ...                   ...            ...         ...         ...       ...
#> 36         HLA-DPB1*59:01           None        None      Common     11516
#> 37         HLA-DPB1*63:01           None        None      Common     11466
#> 38   HLA-DPB1*81:01:01:01 DPB1*81:01:01G DPB1*81:01P      Common     11528
#> 39   HLA-DPB1*85:01:01:01 DPB1*85:01:01G DPB1*85:01P      Common     11461
#> 40  HLA-DPB1*105:01:01:01 DPB1*04:02:01G DPB1*04:02P      Common     11540
#>                                                      FeatureType
#>                                                      <character>
#> 1   UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 2   UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 3   UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 4   UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 5   UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> ...                                                          ...
#> 36  UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 37  UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 38  UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 39  UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 40  UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
```

These obects come with the same set of getter functions as `HLAAllele`s plus an additional set of `R6`-methods that implement an API used mostly in conjunction with the `DR2S` package:

-   `x$get_closest_complete_neighbor(allele)`: Find the closest full-length allele based on the genetic distance between exon 2 sequences.

-   `x$get_reference_sequence(allele)`: Return a `BStringSet` object of a full-length reference sequence for `allele`. If `allele` is not fully known the missing stretches are taken from the allele returned by a call to `get_closest_complete_neighbor()`.
