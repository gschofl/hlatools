# hlatools 0.1.6

* Update for IPD-IMGT/HLA 3.40.2, 3.42.0 and 3.43.0
* Fix parsing MICA and MIC, which was broken since 3.40.0

# hlatools 0.1.5

* Update for IPD-IMGT/HLA 3.40.0 and 3.41.0

# hlatools 0.1.3

* Complete alleles are now defined as featuring the expected set of exons
  and introns instead of the 2 UTR sequences.
* Update for IPD-IMGT/HLA 3.38.0 and 3.39.0
* As of IPD-IMGT/HLA 3.39.0, the locusname for `MICA` and `MICB` has become `MICA` and `MICB` instead of `HLA-MICA` and `HLA-MICB`.

# hlatools 0.1.3

* Include parsing `MICA` and `MICB`
* Update for IPD-IMGT/HLA 3.37.0

# hlatools 0.1.2

* Update for IPD-IMGT/HLA 3.35.0
* Added `HLAGene` method `x$get_extended_reference_set()`.
* Handle HLA-DRB before IPD-IMGT/HLA release 3.26.0.

# hlatools 0.1.1

* Update for IPD-IMGT/HLA 3.33.0
* Update for git2r 0.23.0

# hlatools 0.1.0

## API changes
* The code table functions `nmdp_table()`, `g_table`, and `generate_nmdp_lookup()`
  return `tibble`s instead of `data.table`s.
* Renamed `fetch_IMGTHLA()` to `clone_IMGTHLA()` and `update_IMGHTHLA()` to
  `pull_IMGTHLA()` (the old functions are aliased to the new).

## New or modified functions
* Added `allele_table()` to parse the *hla_nom.txt* file.
* Added `check_IMGTHLA()` to check if the local repo is up-to-date.
* Added `noutr()` to extract sequences except for UTRs.
  
## Minor improvements and bug fixes
* Added a `db_path` argument to `HLAGene()` and `read_hla_xml()` to allow for
  alternative repository locations.
* Remove package `curl` as a dependency.

## Internal changes
* Distance calculations between alleles are now based on all available exons
  instead of only exon 2.
* Added `HLAGene` method `x$calculate_exon_distance_matrix()`.
* The `HLAGene` method `x$get_closest_complete_neighbor()` gains an argument
  `partially = TRUE` that allows for partial allele matching.

# hlatools 0.0.7

* [feature] `db_version()` returns a `numeric_version` object instead of a `character` string (#3).
* [feature] new `hlatools_version()` method for HLAGene objects (#2).
* [bug] fix `locusname()` method for HLAGene objects (#1).
* [feature] add `exon()`, `intron()`, and `utr()` to extract feature sequences.
* [feature] Partial matching of allele names when subsetting `HLAAllele` objects.
* [feature] Added a `NEWS.md` file to track changes to the package.



