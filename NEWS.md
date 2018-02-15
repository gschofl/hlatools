# hlatools 0.0.7.9000

## API changes
* The code table functions `nmdp_table()`, `g_table`, and `generate_nmdp_lookup()`
  return `tibble`s instead of `data.table`s.

## New or modified functions
* Added `allele_table()` to parse the *hla_nom.txt* file.
* Added `check_IMGTHLA()` to check if the local repo is uo-to-date.
* Renamed `fetch_IMGTHLA()` to `clone_IMGTHLA()` and `update_IMGHTHLA()` to
  `pull_IMGTHLA()` (the old functions are aliased to the new).
* Added `noutr()` to extract sequences except for UTRs.
  
## Minor improvements and bug fixes
* Added `db_path` argument to `HLAGene()` and `read_hla_xml()` to allow for
  alternative repositry locations.
* Remove package `curl` as a dependency.

# hlatools 0.0.7

* [feature] `db_version()` returns a `numeric_version` object instead of a `character` string (#3).
* [feature] new `hlatools_version()` method for HLAGene objects (#2).
* [bug] fix `locusname()` method for HLAGene objects (#1).
* [feature] add `exon()`, `intron()`, and `utr()` to extract feature sequences.
* [feature] Partial matching of allele names when subsetting `HLAAllele` objects.
* [feature] Added a `NEWS.md` file to track changes to the package.



