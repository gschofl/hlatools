
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
#> DataFrame with 716 rows and 6 columns
#>           allele_name        g_group     p_group      cwd_status    SeqLen
#>           <character>    <character> <character>     <character> <integer>
#> 1   HLA-DPB1*01:01:01 DPB1*01:01:01G DPB1*01:01P          Common       777
#> 2   HLA-DPB1*01:01:02 DPB1*01:01:02G DPB1*01:01P          Common       777
#> 3   HLA-DPB1*01:01:03           None DPB1*01:01P Not CWD defined       264
#> 4   HLA-DPB1*01:01:04           None DPB1*01:01P Not CWD defined       264
#> 5   HLA-DPB1*01:01:05           None DPB1*01:01P Not CWD defined       264
#> ...               ...            ...         ...             ...       ...
#> 712   HLA-DPB1*608:01           None        None Not CWD defined       264
#> 713   HLA-DPB1*609:01           None        None Not CWD defined       264
#> 714   HLA-DPB1*610:01           None        None Not CWD defined       264
#> 715   HLA-DPB1*611:01           None        None Not CWD defined       264
#> 716   HLA-DPB1*612:01           None        None Not CWD defined       264
#>                  FeatureType
#>                  <character>
#> 1   Exon:Exon:Exon:Exon:Exon
#> 2   Exon:Exon:Exon:Exon:Exon
#> 3                       Exon
#> 4                       Exon
#> 5                       Exon
#> ...                      ...
#> 712                     Exon
#> 713                     Exon
#> 714                     Exon
#> 715                     Exon
#> 716                     Exon
```

A number of accessor functions can be used to work with these objects:

``` r
dpb1_cmpl <- dpb1[hlatools::is_complete(dpb1)]
hlatools::allele_name(dpb1_cmpl)
#>  [1] "HLA-DPB1*02:01:02"    "HLA-DPB1*02:02"       "HLA-DPB1*03:01:01"   
#>  [4] "HLA-DPB1*04:01:01:01" "HLA-DPB1*04:01:01:02" "HLA-DPB1*04:01:31"   
#>  [7] "HLA-DPB1*04:02:01:01" "HLA-DPB1*04:02:01:02" "HLA-DPB1*16:01:01"   
#> [10] "HLA-DPB1*40:01"       "HLA-DPB1*45:01"       "HLA-DPB1*46:01:01"   
#> [13] "HLA-DPB1*59:01"       "HLA-DPB1*104:01"      "HLA-DPB1*398:01"     
#> [16] "HLA-DPB1*463:01"      "HLA-DPB1*464:01"      "HLA-DPB1*584:01"
hlatools::cwd_status(dpb1_cmpl)
#>  [1] "Common"          "Common"          "Common"         
#>  [4] "Common"          "Not CWD defined" "Not CWD defined"
#>  [7] "Common"          "Not CWD defined" "Common"         
#> [10] "Well-Documented" "Common"          "Common"         
#> [13] "Common"          "Well-Documented" "Not CWD defined"
#> [16] "Not CWD defined" "Not CWD defined" "Not CWD defined"
hlatools::sequences(dpb1_cmpl)
#>   A DNAStringSet instance of length 18
#>      width seq                                         names               
#>  [1] 11532 TAATCCCTGTAGATGGGCCA...TATAATCTAATACACTTTAA HLA-DPB1*02:01:02
#>  [2] 11528 TAATCCCTGTAGATGGGCCA...TATAATCTAATACACTTTAA HLA-DPB1*02:02
#>  [3] 11461 TAATCCCTGTAGATGGGCCA...TATAATCTAATACACTTTAA HLA-DPB1*03:01:01
#>  [4] 11526 TAATCCCTGTAGATGGGCCA...TATAATCTAATACACTTTAA HLA-DPB1*04:01:01:01
#>  [5] 11518 TAATCCCTGTAGATGGGCCA...TATAATCTAATACACTTTAA HLA-DPB1*04:01:01:02
#>  ...   ... ...
#> [14] 11467 TAATCCCTGTAGATGGGCCA...TATAATCTAATACACTTTAA HLA-DPB1*104:01
#> [15] 11526 TAATCCCTGTAGATGGGCCA...TATAATCTAATACACTTTAA HLA-DPB1*398:01
#> [16] 11468 TAATCCCTGTAGATGGGCCA...TATAATCTAATACACTTTAA HLA-DPB1*463:01
#> [17] 11526 TAATCCCTGTAGATGGGCCA...TATAATCTAATACACTTTAA HLA-DPB1*464:01
#> [18] 11457 TAATCCCTGTAGATGGGCCA...TATAATCTAATACACTTTAA HLA-DPB1*584:01
```

At a slightly higher level a `HLAAllele` can be encapsulated in an `R6`-based `HLAGene` object:

``` r
x <- hlatools::HLAGene("DPB1")
x
#> IMGT/HLA database <3.26.0>; Locus <HLA-DPB1>
#> An object of class 'HLAAllele'
#> DataFrame with 716 rows and 6 columns
#>           allele_name        g_group     p_group      cwd_status    SeqLen
#>           <character>    <character> <character>     <character> <integer>
#> 1   HLA-DPB1*01:01:01 DPB1*01:01:01G DPB1*01:01P          Common       777
#> 2   HLA-DPB1*01:01:02 DPB1*01:01:02G DPB1*01:01P          Common       777
#> 3   HLA-DPB1*01:01:03           None DPB1*01:01P Not CWD defined       264
#> 4   HLA-DPB1*01:01:04           None DPB1*01:01P Not CWD defined       264
#> 5   HLA-DPB1*01:01:05           None DPB1*01:01P Not CWD defined       264
#> ...               ...            ...         ...             ...       ...
#> 712   HLA-DPB1*608:01           None        None Not CWD defined       264
#> 713   HLA-DPB1*609:01           None        None Not CWD defined       264
#> 714   HLA-DPB1*610:01           None        None Not CWD defined       264
#> 715   HLA-DPB1*611:01           None        None Not CWD defined       264
#> 716   HLA-DPB1*612:01           None        None Not CWD defined       264
#>                  FeatureType
#>                  <character>
#> 1   Exon:Exon:Exon:Exon:Exon
#> 2   Exon:Exon:Exon:Exon:Exon
#> 3                       Exon
#> 4                       Exon
#> 5                       Exon
#> ...                      ...
#> 712                     Exon
#> 713                     Exon
#> 714                     Exon
#> 715                     Exon
#> 716                     Exon

hlatools::db_version(x)
#> [1] "3.26.0"
hlatools::locusname(x)
#> [1] "HLA-DPB1"
x[hlatools::cwd_status(x) == "Common"]
#> An object of class 'HLAAllele'
#> DataFrame with 40 rows and 6 columns
#>           allele_name        g_group     p_group  cwd_status    SeqLen
#>           <character>    <character> <character> <character> <integer>
#> 1   HLA-DPB1*01:01:01 DPB1*01:01:01G DPB1*01:01P      Common       777
#> 2   HLA-DPB1*01:01:02 DPB1*01:01:02G DPB1*01:01P      Common       777
#> 3   HLA-DPB1*02:01:02 DPB1*02:01:02G DPB1*02:01P      Common     11532
#> 4      HLA-DPB1*02:02 DPB1*02:02:01G DPB1*02:02P      Common     11528
#> 5   HLA-DPB1*03:01:01 DPB1*03:01:01G DPB1*03:01P      Common     11461
#> ...               ...            ...         ...         ...       ...
#> 36     HLA-DPB1*59:01           None        None      Common     11516
#> 37     HLA-DPB1*63:01           None        None      Common       255
#> 38     HLA-DPB1*81:01           None        None      Common       264
#> 39     HLA-DPB1*85:01           None        None      Common       546
#> 40    HLA-DPB1*105:01 DPB1*04:02:01G DPB1*04:02P      Common       777
#>                                                      FeatureType
#>                                                      <character>
#> 1                                       Exon:Exon:Exon:Exon:Exon
#> 2                                       Exon:Exon:Exon:Exon:Exon
#> 3   UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 4   UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 5   UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> ...                                                          ...
#> 36  UTR:Exon:Intron:Exon:Intron:Exon:Intron:Exon:Intron:Exon:UTR
#> 37                                                          Exon
#> 38                                                          Exon
#> 39                                                     Exon:Exon
#> 40                                      Exon:Exon:Exon:Exon:Exon
```

These obects come with the same set of getter functions as `HLAAllele`s plus an additional set of `R6`-methods that implement an API used mostly in conjunction with the `DR2S` package:

-   `x$get_closest_complete_neighbor(allele)`: Find the closest full-length allele based on the genetic distance between exon 2 sequences.

-   `x$get_reference_sequence(allele)`: Return a `BStringSet` object of a full-length reference sequence for `allele`. If `allele` is not fully known the missing stretches are taken from the allele returned by a call to `get_closest_complete_neighbor()`.
