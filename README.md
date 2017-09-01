
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
```

A number of accessor functions can be used to work with these objects:

``` r
dpb1_cmpl <- dpb1[hlatools::is_complete(dpb1)]
hlatools::allele_name(dpb1_cmpl)
#>   [1] "HLA-DPB1*01:01:01:01"  "HLA-DPB1*01:01:01:02" 
#>   [3] "HLA-DPB1*01:01:01:03"  "HLA-DPB1*01:01:01:04" 
#>   [5] "HLA-DPB1*01:01:02:01"  "HLA-DPB1*01:01:02:02" 
#>   [7] "HLA-DPB1*02:01:02:01"  "HLA-DPB1*02:01:02:02" 
#>   [9] "HLA-DPB1*02:01:02:03"  "HLA-DPB1*02:01:02:04" 
#>  [11] "HLA-DPB1*02:01:02:05"  "HLA-DPB1*02:01:02:06" 
#>  [13] "HLA-DPB1*02:01:02:07"  "HLA-DPB1*02:01:02:08" 
#>  [15] "HLA-DPB1*02:01:02:09"  "HLA-DPB1*02:01:02:10" 
#>  [17] "HLA-DPB1*02:01:02:11"  "HLA-DPB1*02:01:02:12" 
#>  [19] "HLA-DPB1*02:01:02:13"  "HLA-DPB1*02:01:02:14" 
#>  [21] "HLA-DPB1*02:01:02:15"  "HLA-DPB1*02:01:02:16" 
#>  [23] "HLA-DPB1*02:01:02:17"  "HLA-DPB1*02:01:02:18" 
#>  [25] "HLA-DPB1*02:01:02:19"  "HLA-DPB1*02:01:02:20" 
#>  [27] "HLA-DPB1*02:01:02:21"  "HLA-DPB1*02:01:02:22" 
#>  [29] "HLA-DPB1*02:01:03"     "HLA-DPB1*02:01:04"    
#>  [31] "HLA-DPB1*02:01:08"     "HLA-DPB1*02:02:01:01" 
#>  [33] "HLA-DPB1*02:02:01:02"  "HLA-DPB1*02:02:01:03" 
#>  [35] "HLA-DPB1*02:02:01:04"  "HLA-DPB1*03:01:01:01" 
#>  [37] "HLA-DPB1*03:01:01:02"  "HLA-DPB1*03:01:01:03" 
#>  [39] "HLA-DPB1*03:01:01:04"  "HLA-DPB1*03:01:01:05" 
#>  [41] "HLA-DPB1*03:01:01:06"  "HLA-DPB1*03:01:01:07" 
#>  [43] "HLA-DPB1*04:01:01:01"  "HLA-DPB1*04:01:01:02" 
#>  [45] "HLA-DPB1*04:01:01:03"  "HLA-DPB1*04:01:01:04" 
#>  [47] "HLA-DPB1*04:01:01:05"  "HLA-DPB1*04:01:01:06" 
#>  [49] "HLA-DPB1*04:01:01:07"  "HLA-DPB1*04:01:01:08" 
#>  [51] "HLA-DPB1*04:01:01:09"  "HLA-DPB1*04:01:01:10" 
#>  [53] "HLA-DPB1*04:01:01:11"  "HLA-DPB1*04:01:03"    
#>  [55] "HLA-DPB1*04:01:04"     "HLA-DPB1*04:01:31"    
#>  [57] "HLA-DPB1*04:02:01:01"  "HLA-DPB1*04:02:01:02" 
#>  [59] "HLA-DPB1*04:02:01:03"  "HLA-DPB1*04:02:01:04" 
#>  [61] "HLA-DPB1*04:02:01:05"  "HLA-DPB1*04:02:01:06" 
#>  [63] "HLA-DPB1*04:02:01:07"  "HLA-DPB1*04:02:01:08" 
#>  [65] "HLA-DPB1*05:01:01:01"  "HLA-DPB1*05:01:01:02" 
#>  [67] "HLA-DPB1*05:01:01:03"  "HLA-DPB1*05:01:01:04" 
#>  [69] "HLA-DPB1*06:01:01:01"  "HLA-DPB1*06:01:01:02" 
#>  [71] "HLA-DPB1*09:01:01"     "HLA-DPB1*09:01:02"    
#>  [73] "HLA-DPB1*10:01:01:01"  "HLA-DPB1*10:01:01:02" 
#>  [75] "HLA-DPB1*11:01:01"     "HLA-DPB1*13:01:01"    
#>  [77] "HLA-DPB1*13:01:02"     "HLA-DPB1*14:01:01:01" 
#>  [79] "HLA-DPB1*14:01:01:02"  "HLA-DPB1*15:01:01"    
#>  [81] "HLA-DPB1*16:01:01"     "HLA-DPB1*17:01:01:01" 
#>  [83] "HLA-DPB1*17:01:01:02"  "HLA-DPB1*18:01"       
#>  [85] "HLA-DPB1*19:01"        "HLA-DPB1*20:01:01"    
#>  [87] "HLA-DPB1*21:01"        "HLA-DPB1*22:01"       
#>  [89] "HLA-DPB1*23:01:01"     "HLA-DPB1*24:01"       
#>  [91] "HLA-DPB1*25:01"        "HLA-DPB1*26:01:02"    
#>  [93] "HLA-DPB1*27:01"        "HLA-DPB1*28:01"       
#>  [95] "HLA-DPB1*29:01"        "HLA-DPB1*30:01"       
#>  [97] "HLA-DPB1*31:01"        "HLA-DPB1*33:01:01:01" 
#>  [99] "HLA-DPB1*33:01:01:02"  "HLA-DPB1*34:01"       
#> [101] "HLA-DPB1*35:01:01"     "HLA-DPB1*36:01"       
#> [103] "HLA-DPB1*38:01"        "HLA-DPB1*39:01"       
#> [105] "HLA-DPB1*40:01"        "HLA-DPB1*45:01"       
#> [107] "HLA-DPB1*46:01:01"     "HLA-DPB1*47:01:01:01" 
#> [109] "HLA-DPB1*47:01:01:02"  "HLA-DPB1*47:01:01:03" 
#> [111] "HLA-DPB1*48:01"        "HLA-DPB1*49:01:01:01" 
#> [113] "HLA-DPB1*49:01:01:02"  "HLA-DPB1*50:01"       
#> [115] "HLA-DPB1*51:01"        "HLA-DPB1*55:01:01:01" 
#> [117] "HLA-DPB1*55:01:01:02"  "HLA-DPB1*57:01"       
#> [119] "HLA-DPB1*59:01"        "HLA-DPB1*63:01"       
#> [121] "HLA-DPB1*69:01"        "HLA-DPB1*71:01:01"    
#> [123] "HLA-DPB1*72:01:01:01"  "HLA-DPB1*72:01:01:02" 
#> [125] "HLA-DPB1*72:01:01:03"  "HLA-DPB1*76:01"       
#> [127] "HLA-DPB1*77:01"        "HLA-DPB1*78:01"       
#> [129] "HLA-DPB1*80:01"        "HLA-DPB1*81:01:01:01" 
#> [131] "HLA-DPB1*81:01:01:02"  "HLA-DPB1*85:01:01:01" 
#> [133] "HLA-DPB1*85:01:01:02"  "HLA-DPB1*88:01"       
#> [135] "HLA-DPB1*90:01:01"     "HLA-DPB1*91:01:01:01" 
#> [137] "HLA-DPB1*91:01:01:02"  "HLA-DPB1*93:01"       
#> [139] "HLA-DPB1*104:01:01:01" "HLA-DPB1*104:01:01:02"
#> [141] "HLA-DPB1*104:01:01:03" "HLA-DPB1*105:01:01:01"
#> [143] "HLA-DPB1*105:01:01:02" "HLA-DPB1*105:01:01:03"
#> [145] "HLA-DPB1*106:01"       "HLA-DPB1*109:01"      
#> [147] "HLA-DPB1*115:01"       "HLA-DPB1*124:01"      
#> [149] "HLA-DPB1*126:01"       "HLA-DPB1*127:01"      
#> [151] "HLA-DPB1*128:01"       "HLA-DPB1*130:01"      
#> [153] "HLA-DPB1*190:01"       "HLA-DPB1*260:01"      
#> [155] "HLA-DPB1*398:01"       "HLA-DPB1*414:01"      
#> [157] "HLA-DPB1*416:01:01:01" "HLA-DPB1*416:01:01:02"
#> [159] "HLA-DPB1*463:01:01:01" "HLA-DPB1*463:01:01:02"
#> [161] "HLA-DPB1*463:01:01:03" "HLA-DPB1*464:01"      
#> [163] "HLA-DPB1*535:01"       "HLA-DPB1*584:01"      
#> [165] "HLA-DPB1*648:01:01:01" "HLA-DPB1*648:01:01:02"
hlatools::cwd_status(dpb1_cmpl)
#>   [1] "Common"          "Not CWD defined" "Not CWD defined"
#>   [4] "Not CWD defined" "Common"          "Not CWD defined"
#>   [7] "Common"          "Not CWD defined" "Not CWD defined"
#>  [10] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [13] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [16] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [19] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [22] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [25] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [28] "Not CWD defined" "Not CWD defined" "Well-Documented"
#>  [31] "Not CWD defined" "Common"          "Not CWD defined"
#>  [34] "Not CWD defined" "Not CWD defined" "Common"         
#>  [37] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [40] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [43] "Common"          "Not CWD defined" "Not CWD defined"
#>  [46] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [49] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [52] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [55] "Not CWD defined" "Not CWD defined" "Common"         
#>  [58] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [61] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [64] "Not CWD defined" "Common"          "Not CWD defined"
#>  [67] "Not CWD defined" "Not CWD defined" "Common"         
#>  [70] "Not CWD defined" "Common"          "Not CWD defined"
#>  [73] "Common"          "Not CWD defined" "Common"         
#>  [76] "Common"          "Not CWD defined" "Common"         
#>  [79] "Not CWD defined" "Common"          "Common"         
#>  [82] "Common"          "Not CWD defined" "Common"         
#>  [85] "Common"          "Common"          "Common"         
#>  [88] "Well-Documented" "Common"          "Not CWD defined"
#>  [91] "Not CWD defined" "Common"          "Common"         
#>  [94] "Common"          "Common"          "Common"         
#>  [97] "Common"          "Well-Documented" "Not CWD defined"
#> [100] "Common"          "Common"          "Well-Documented"
#> [103] "Well-Documented" "Common"          "Well-Documented"
#> [106] "Common"          "Common"          "Well-Documented"
#> [109] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [112] "Well-Documented" "Not CWD defined" "Not CWD defined"
#> [115] "Common"          "Common"          "Not CWD defined"
#> [118] "Well-Documented" "Common"          "Common"         
#> [121] "Well-Documented" "Not CWD defined" "Well-Documented"
#> [124] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [127] "Not CWD defined" "Well-Documented" "Not CWD defined"
#> [130] "Common"          "Not CWD defined" "Common"         
#> [133] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [136] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [139] "Well-Documented" "Not CWD defined" "Not CWD defined"
#> [142] "Common"          "Not CWD defined" "Not CWD defined"
#> [145] "Well-Documented" "Not CWD defined" "Not CWD defined"
#> [148] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [151] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [154] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [157] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [160] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [163] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [166] "Not CWD defined"
hlatools::sequences(dpb1_cmpl)
#>   A DNAStringSet instance of length 166
#>       width seq                                        names               
#>   [1] 11468 TAATCCCTGTAGATGGGCCA...ATAATCTAATACACTTTAA HLA-DPB1*01:01:01:01
#>   [2] 11468 TAATCCCTGTAGATGGGCCA...ATAATCTAATACACTTTAA HLA-DPB1*01:01:01:02
#>   [3] 11468 TAATCCCTGTAGATGGGCCA...ATAATCTAATACACTTTAA HLA-DPB1*01:01:01:03
#>   [4] 11468 TAATCCCTGTAGATGGGCCA...ATAATCTAATACACTTTAA HLA-DPB1*01:01:01:04
#>   [5] 11469 TAATCCCTGTAGATGGGCCA...ATAATCTAATACACTTTAA HLA-DPB1*01:01:02:01
#>   ...   ... ...
#> [162] 11526 TAATCCCTGTAGATGGGCCA...ATAATCTAATACACTTTAA HLA-DPB1*464:01
#> [163] 11512 TAATCCCTGTAGATGGGCCA...ATAATCTAATACACTTTAA HLA-DPB1*535:01
#> [164] 11457 TAATCCCTGTAGATGGGCCA...ATAATCTAATACACTTTAA HLA-DPB1*584:01
#> [165] 11466 TAATCCCTGTAGATGGGCCA...ATAATCTAATACACTTTAA HLA-DPB1*648:01:0...
#> [166] 11468 TAATCCCTGTAGATGGGCCA...ATAATCTAATACACTTTAA HLA-DPB1*648:01:0...
```

At a slightly higher level a `HLAAllele` can be encapsulated in an `R6`-based `HLAGene` object:

``` r
x <- hlatools::HLAGene("DPB1")
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
