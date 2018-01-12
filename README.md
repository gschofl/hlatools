
<!-- README.md is generated from README.Rmd. Please edit that file -->
hlatools
========

A collection of tools to work with IPD-IMGT/HLA data

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
#> An object of class 'HLAAllele'
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
#>  [29] "HLA-DPB1*02:01:02:23"  "HLA-DPB1*02:01:02:24" 
#>  [31] "HLA-DPB1*02:01:02:25"  "HLA-DPB1*02:01:02:26" 
#>  [33] "HLA-DPB1*02:01:02:27"  "HLA-DPB1*02:01:02:28" 
#>  [35] "HLA-DPB1*02:01:02:29"  "HLA-DPB1*02:01:02:30" 
#>  [37] "HLA-DPB1*02:01:02:31"  "HLA-DPB1*02:01:02:32" 
#>  [39] "HLA-DPB1*02:01:02:33"  "HLA-DPB1*02:01:02:34" 
#>  [41] "HLA-DPB1*02:01:02:35"  "HLA-DPB1*02:01:02:36" 
#>  [43] "HLA-DPB1*02:01:02:37"  "HLA-DPB1*02:01:03"    
#>  [45] "HLA-DPB1*02:01:04"     "HLA-DPB1*02:01:08"    
#>  [47] "HLA-DPB1*02:01:15"     "HLA-DPB1*02:01:20"    
#>  [49] "HLA-DPB1*02:02:01:01"  "HLA-DPB1*02:02:01:02" 
#>  [51] "HLA-DPB1*02:02:01:03"  "HLA-DPB1*02:02:01:04" 
#>  [53] "HLA-DPB1*02:02:01:05"  "HLA-DPB1*02:02:01:06" 
#>  [55] "HLA-DPB1*02:02:01:07"  "HLA-DPB1*03:01:01:01" 
#>  [57] "HLA-DPB1*03:01:01:02"  "HLA-DPB1*03:01:01:03" 
#>  [59] "HLA-DPB1*03:01:01:04"  "HLA-DPB1*03:01:01:05" 
#>  [61] "HLA-DPB1*03:01:01:06"  "HLA-DPB1*03:01:01:07" 
#>  [63] "HLA-DPB1*03:01:01:08"  "HLA-DPB1*03:01:01:09" 
#>  [65] "HLA-DPB1*03:01:01:10"  "HLA-DPB1*04:01:01:01" 
#>  [67] "HLA-DPB1*04:01:01:02"  "HLA-DPB1*04:01:01:03" 
#>  [69] "HLA-DPB1*04:01:01:04"  "HLA-DPB1*04:01:01:05" 
#>  [71] "HLA-DPB1*04:01:01:06"  "HLA-DPB1*04:01:01:07" 
#>  [73] "HLA-DPB1*04:01:01:08"  "HLA-DPB1*04:01:01:09" 
#>  [75] "HLA-DPB1*04:01:01:10"  "HLA-DPB1*04:01:01:11" 
#>  [77] "HLA-DPB1*04:01:01:12"  "HLA-DPB1*04:01:01:13" 
#>  [79] "HLA-DPB1*04:01:01:14"  "HLA-DPB1*04:01:01:15" 
#>  [81] "HLA-DPB1*04:01:01:16"  "HLA-DPB1*04:01:01:17" 
#>  [83] "HLA-DPB1*04:01:01:18"  "HLA-DPB1*04:01:01:19" 
#>  [85] "HLA-DPB1*04:01:01:20"  "HLA-DPB1*04:01:01:21" 
#>  [87] "HLA-DPB1*04:01:01:22"  "HLA-DPB1*04:01:01:23" 
#>  [89] "HLA-DPB1*04:01:01:24N" "HLA-DPB1*04:01:01:25" 
#>  [91] "HLA-DPB1*04:01:01:26"  "HLA-DPB1*04:01:03"    
#>  [93] "HLA-DPB1*04:01:04:01"  "HLA-DPB1*04:01:04:02" 
#>  [95] "HLA-DPB1*04:01:31"     "HLA-DPB1*04:01:33"    
#>  [97] "HLA-DPB1*04:01:34"     "HLA-DPB1*04:01:35"    
#>  [99] "HLA-DPB1*04:01:36"     "HLA-DPB1*04:02:01:01" 
#> [101] "HLA-DPB1*04:02:01:02"  "HLA-DPB1*04:02:01:03" 
#> [103] "HLA-DPB1*04:02:01:04"  "HLA-DPB1*04:02:01:05" 
#> [105] "HLA-DPB1*04:02:01:06"  "HLA-DPB1*04:02:01:07" 
#> [107] "HLA-DPB1*04:02:01:08"  "HLA-DPB1*04:02:01:09" 
#> [109] "HLA-DPB1*04:02:10"     "HLA-DPB1*05:01:01:01" 
#> [111] "HLA-DPB1*05:01:01:02"  "HLA-DPB1*05:01:01:03" 
#> [113] "HLA-DPB1*05:01:01:04"  "HLA-DPB1*05:01:01:05" 
#> [115] "HLA-DPB1*05:01:01:06"  "HLA-DPB1*05:01:01:07" 
#> [117] "HLA-DPB1*05:01:01:08"  "HLA-DPB1*05:01:01:09" 
#> [119] "HLA-DPB1*06:01:01:01"  "HLA-DPB1*06:01:01:02" 
#> [121] "HLA-DPB1*06:01:01:03"  "HLA-DPB1*06:01:04"    
#> [123] "HLA-DPB1*09:01:01"     "HLA-DPB1*09:01:02"    
#> [125] "HLA-DPB1*10:01:01:01"  "HLA-DPB1*10:01:01:02" 
#> [127] "HLA-DPB1*11:01:01"     "HLA-DPB1*13:01:01:01" 
#> [129] "HLA-DPB1*13:01:01:02"  "HLA-DPB1*13:01:01:03" 
#> [131] "HLA-DPB1*13:01:01:04"  "HLA-DPB1*13:01:02"    
#> [133] "HLA-DPB1*14:01:01:01"  "HLA-DPB1*14:01:01:02" 
#> [135] "HLA-DPB1*14:01:01:03"  "HLA-DPB1*15:01:01"    
#> [137] "HLA-DPB1*16:01:01"     "HLA-DPB1*17:01:01:01" 
#> [139] "HLA-DPB1*17:01:01:02"  "HLA-DPB1*18:01"       
#> [141] "HLA-DPB1*19:01:01:01"  "HLA-DPB1*19:01:01:02" 
#> [143] "HLA-DPB1*20:01:01"     "HLA-DPB1*20:01:04"    
#> [145] "HLA-DPB1*21:01"        "HLA-DPB1*22:01"       
#> [147] "HLA-DPB1*23:01:01"     "HLA-DPB1*24:01"       
#> [149] "HLA-DPB1*25:01"        "HLA-DPB1*26:01:02"    
#> [151] "HLA-DPB1*27:01"        "HLA-DPB1*28:01"       
#> [153] "HLA-DPB1*29:01"        "HLA-DPB1*30:01:01:01" 
#> [155] "HLA-DPB1*30:01:01:02"  "HLA-DPB1*31:01"       
#> [157] "HLA-DPB1*33:01:01:01"  "HLA-DPB1*33:01:01:02" 
#> [159] "HLA-DPB1*33:01:01:03"  "HLA-DPB1*33:01:01:04" 
#> [161] "HLA-DPB1*33:01:01:05"  "HLA-DPB1*34:01:01:01" 
#> [163] "HLA-DPB1*34:01:01:02"  "HLA-DPB1*35:01:01"    
#> [165] "HLA-DPB1*36:01"        "HLA-DPB1*38:01"       
#> [167] "HLA-DPB1*39:01:01:01"  "HLA-DPB1*39:01:01:02" 
#> [169] "HLA-DPB1*39:01:01:03"  "HLA-DPB1*39:01:01:04" 
#> [171] "HLA-DPB1*39:01:02"     "HLA-DPB1*40:01:01:01" 
#> [173] "HLA-DPB1*40:01:01:02"  "HLA-DPB1*41:01:01:01" 
#> [175] "HLA-DPB1*41:01:01:02"  "HLA-DPB1*45:01"       
#> [177] "HLA-DPB1*46:01:01"     "HLA-DPB1*47:01:01:01" 
#> [179] "HLA-DPB1*47:01:01:02"  "HLA-DPB1*47:01:01:03" 
#> [181] "HLA-DPB1*48:01"        "HLA-DPB1*49:01:01:01" 
#> [183] "HLA-DPB1*49:01:01:02"  "HLA-DPB1*49:01:01:03" 
#> [185] "HLA-DPB1*50:01"        "HLA-DPB1*51:01"       
#> [187] "HLA-DPB1*55:01:01:01"  "HLA-DPB1*55:01:01:02" 
#> [189] "HLA-DPB1*57:01"        "HLA-DPB1*59:01"       
#> [191] "HLA-DPB1*63:01"        "HLA-DPB1*69:01:01:01" 
#> [193] "HLA-DPB1*69:01:01:02"  "HLA-DPB1*71:01:01"    
#> [195] "HLA-DPB1*72:01:01:01"  "HLA-DPB1*72:01:01:02" 
#> [197] "HLA-DPB1*72:01:01:03"  "HLA-DPB1*76:01"       
#> [199] "HLA-DPB1*77:01"        "HLA-DPB1*78:01"       
#> [201] "HLA-DPB1*80:01"        "HLA-DPB1*81:01:01:01" 
#> [203] "HLA-DPB1*81:01:01:02"  "HLA-DPB1*85:01:01:01" 
#> [205] "HLA-DPB1*85:01:01:02"  "HLA-DPB1*88:01"       
#> [207] "HLA-DPB1*90:01:01"     "HLA-DPB1*91:01:01:01" 
#> [209] "HLA-DPB1*91:01:01:02"  "HLA-DPB1*93:01"       
#> [211] "HLA-DPB1*104:01:01:01" "HLA-DPB1*104:01:01:02"
#> [213] "HLA-DPB1*104:01:01:03" "HLA-DPB1*104:01:01:04"
#> [215] "HLA-DPB1*104:01:01:05" "HLA-DPB1*105:01:01:01"
#> [217] "HLA-DPB1*105:01:01:02" "HLA-DPB1*105:01:01:03"
#> [219] "HLA-DPB1*105:01:01:04" "HLA-DPB1*105:01:01:05"
#> [221] "HLA-DPB1*105:01:01:06" "HLA-DPB1*105:01:01:07"
#> [223] "HLA-DPB1*105:01:01:08" "HLA-DPB1*106:01"      
#> [225] "HLA-DPB1*109:01"       "HLA-DPB1*115:01"      
#> [227] "HLA-DPB1*124:01:01:01" "HLA-DPB1*124:01:01:02"
#> [229] "HLA-DPB1*126:01:01:01" "HLA-DPB1*126:01:01:02"
#> [231] "HLA-DPB1*127:01"       "HLA-DPB1*128:01"      
#> [233] "HLA-DPB1*130:01"       "HLA-DPB1*131:01"      
#> [235] "HLA-DPB1*138:01"       "HLA-DPB1*141:01"      
#> [237] "HLA-DPB1*162:01:02"    "HLA-DPB1*178:01"      
#> [239] "HLA-DPB1*184:01"       "HLA-DPB1*187:01"      
#> [241] "HLA-DPB1*190:01"       "HLA-DPB1*191:01"      
#> [243] "HLA-DPB1*208:01"       "HLA-DPB1*236:01:01"   
#> [245] "HLA-DPB1*244:01"       "HLA-DPB1*260:01"      
#> [247] "HLA-DPB1*328:01N"      "HLA-DPB1*398:01"      
#> [249] "HLA-DPB1*414:01:01:01" "HLA-DPB1*414:01:01:02"
#> [251] "HLA-DPB1*415:01"       "HLA-DPB1*416:01:01:01"
#> [253] "HLA-DPB1*416:01:01:02" "HLA-DPB1*463:01:01:01"
#> [255] "HLA-DPB1*463:01:01:02" "HLA-DPB1*463:01:01:03"
#> [257] "HLA-DPB1*464:01"       "HLA-DPB1*498:01"      
#> [259] "HLA-DPB1*535:01"       "HLA-DPB1*572:01"      
#> [261] "HLA-DPB1*584:01"       "HLA-DPB1*648:01:01:01"
#> [263] "HLA-DPB1*648:01:01:02" "HLA-DPB1*649:01"      
#> [265] "HLA-DPB1*650:01"       "HLA-DPB1*651:01"      
#> [267] "HLA-DPB1*652:01"       "HLA-DPB1*653:01"      
#> [269] "HLA-DPB1*666:01"       "HLA-DPB1*667:01"      
#> [271] "HLA-DPB1*668:01:01:01" "HLA-DPB1*668:01:01:02"
#> [273] "HLA-DPB1*669:01"       "HLA-DPB1*670:01"      
#> [275] "HLA-DPB1*671:01"       "HLA-DPB1*672:01"      
#> [277] "HLA-DPB1*673:01"       "HLA-DPB1*674:01"      
#> [279] "HLA-DPB1*675:01"       "HLA-DPB1*676:01"      
#> [281] "HLA-DPB1*677:01"       "HLA-DPB1*678:01"
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
#>  [28] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [31] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [34] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [37] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [40] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [43] "Not CWD defined" "Not CWD defined" "Well-Documented"
#>  [46] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [49] "Common"          "Not CWD defined" "Not CWD defined"
#>  [52] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [55] "Not CWD defined" "Common"          "Not CWD defined"
#>  [58] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [61] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [64] "Not CWD defined" "Not CWD defined" "Common"         
#>  [67] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [70] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [73] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [76] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [79] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [82] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [85] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [88] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [91] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [94] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#>  [97] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [100] "Common"          "Not CWD defined" "Not CWD defined"
#> [103] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [106] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [109] "Not CWD defined" "Common"          "Not CWD defined"
#> [112] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [115] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [118] "Not CWD defined" "Common"          "Not CWD defined"
#> [121] "Not CWD defined" "Not CWD defined" "Common"         
#> [124] "Not CWD defined" "Common"          "Not CWD defined"
#> [127] "Common"          "Common"          "Not CWD defined"
#> [130] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [133] "Common"          "Not CWD defined" "Not CWD defined"
#> [136] "Common"          "Common"          "Common"         
#> [139] "Not CWD defined" "Common"          "Common"         
#> [142] "Not CWD defined" "Common"          "Not CWD defined"
#> [145] "Common"          "Well-Documented" "Common"         
#> [148] "Not CWD defined" "Not CWD defined" "Common"         
#> [151] "Common"          "Common"          "Common"         
#> [154] "Common"          "Not CWD defined" "Common"         
#> [157] "Well-Documented" "Not CWD defined" "Not CWD defined"
#> [160] "Not CWD defined" "Not CWD defined" "Common"         
#> [163] "Not CWD defined" "Common"          "Well-Documented"
#> [166] "Well-Documented" "Common"          "Not CWD defined"
#> [169] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [172] "Well-Documented" "Not CWD defined" "Not CWD defined"
#> [175] "Not CWD defined" "Common"          "Common"         
#> [178] "Well-Documented" "Not CWD defined" "Not CWD defined"
#> [181] "Not CWD defined" "Well-Documented" "Not CWD defined"
#> [184] "Not CWD defined" "Not CWD defined" "Common"         
#> [187] "Common"          "Not CWD defined" "Well-Documented"
#> [190] "Common"          "Common"          "Well-Documented"
#> [193] "Not CWD defined" "Not CWD defined" "Well-Documented"
#> [196] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [199] "Not CWD defined" "Well-Documented" "Not CWD defined"
#> [202] "Common"          "Not CWD defined" "Common"         
#> [205] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [208] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [211] "Well-Documented" "Not CWD defined" "Not CWD defined"
#> [214] "Not CWD defined" "Not CWD defined" "Common"         
#> [217] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [220] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [223] "Not CWD defined" "Well-Documented" "Not CWD defined"
#> [226] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [229] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [232] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [235] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [238] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [241] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [244] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [247] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [250] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [253] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [256] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [259] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [262] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [265] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [268] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [271] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [274] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [277] "Not CWD defined" "Not CWD defined" "Not CWD defined"
#> [280] "Not CWD defined" "Not CWD defined" "Not CWD defined"
hlatools::sequences(dpb1_cmpl)
#>   A DNAStringSet instance of length 282
#>       width seq                                        names               
#>   [1] 11468 TAATCCCTGTAGATGGGCCA...ATAATCTAATACACTTTAA HLA-DPB1*01:01:01:01
#>   [2] 11468 TAATCCCTGTAGATGGGCCA...ATAATCTAATACACTTTAA HLA-DPB1*01:01:01:02
#>   [3] 11468 TAATCCCTGTAGATGGGCCA...ATAATCTAATACACTTTAA HLA-DPB1*01:01:01:03
#>   [4] 11468 TAATCCCTGTAGATGGGCCA...ATAATCTAATACACTTTAA HLA-DPB1*01:01:01:04
#>   [5] 11469 TAATCCCTGTAGATGGGCCA...ATAATCTAATACACTTTAA HLA-DPB1*01:01:02:01
#>   ...   ... ...
#> [278] 11516 TAATCCCTGTAGATGGGCCA...ATAATCTAATACACTTTAA HLA-DPB1*674:01
#> [279] 11461 TAATCCCTGTAGATGGGCCA...ATAATCTAATACACTTTAA HLA-DPB1*675:01
#> [280] 11461 TAATCCCTGTAGATGGGCCA...ATAATCTAATACACTTTAA HLA-DPB1*676:01
#> [281] 11526 TAATCCCTGTAGATGGGCCA...ATAATCTAATACACTTTAA HLA-DPB1*677:01
#> [282] 11532 TAATCCCTGTAGATGGGCCA...ATAATCTAATACACTTTAA HLA-DPB1*678:01
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

### Session Info

    #> Session info -------------------------------------------------------------
    #>  setting  value                       
    #>  version  R version 3.4.3 (2017-11-30)
    #>  system   x86_64, linux-gnu           
    #>  ui       X11                         
    #>  language en_GB:en                    
    #>  collate  en_GB.UTF-8                 
    #>  tz       Europe/Berlin               
    #>  date     2018-01-12
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
    #>  digest             0.6.13    2017-12-14 cran (@0.6.13)             
    #>  evaluate           0.10.1    2017-06-24 cran (@0.10.1)             
    #>  foreach            1.4.4     2017-12-12 cran (@1.4.4)              
    #>  GenomeInfoDb       1.14.0    2017-11-17 Bioconductor               
    #>  GenomeInfoDbData   0.99.1    2017-11-17 Bioconductor               
    #>  GenomicRanges      1.30.0    2017-11-17 Bioconductor               
    #>  git2r              0.21.0    2018-01-04 cran (@0.21.0)             
    #>  graphics         * 3.4.3     2017-12-01 local                      
    #>  grDevices        * 3.4.3     2017-12-01 local                      
    #>  hlatools         * 0.0.6     2018-01-12 local (gschofl/hlatools@NA)
    #>  htmltools          0.3.6     2017-04-28 CRAN (R 3.4.2)             
    #>  IRanges            2.12.0    2017-11-17 Bioconductor               
    #>  iterators          1.0.9     2017-12-12 cran (@1.0.9)              
    #>  knitr              1.18      2017-12-27 cran (@1.18)               
    #>  magrittr           1.5       2014-11-22 CRAN (R 3.4.2)             
    #>  memoise            1.1.0     2017-04-21 CRAN (R 3.4.2)             
    #>  methods          * 3.4.3     2017-12-01 local                      
    #>  parallel           3.4.3     2017-12-01 local                      
    #>  R6                 2.2.2     2017-06-17 cran (@2.2.2)              
    #>  Rcpp               0.12.14   2017-11-23 cran (@0.12.14)            
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
    #>  xml2               1.1.1     2017-01-24 CRAN (R 3.4.2)             
    #>  XVector            0.18.0    2017-11-17 Bioconductor               
    #>  yaml               2.1.16    2017-12-12 cran (@2.1.16)             
    #>  zlibbioc           1.24.0    2017-11-17 Bioconductor
