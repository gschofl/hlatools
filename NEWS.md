# hlatools 0.0.7.9000

* [feature] Added `check_IMGTHLA()` to check if the local repo is uo-to-date.
* Renamed `fetch_IMGTHLA()` to `clone_IMGTHLA()` and `update_IMGHTHLA()` to
  `pull_IMGTHLA()` (the old functions are aliased to the new)

# hlatools 0.0.7

* [feature] `db_version()` returns a `numeric_version` object instead of a `character` string (#3).
* [feature] new `hlatools_version()` method for HLAGene objects (#2).
* [bug] fix `locusname()` method for HLAGene objects (#1).
* [feature] add `exon()`, `intron()`, and `utr()` to extract feature sequences
* [feature] Partial matching of allele names when subsetting `HLAAllele` objects.
* [feature] Added a `NEWS.md` file to track changes to the package.



