MaTra
=====

MaTra is prototype tool written in D to extract MassTrace from a [mzDB](https://github.com/mzdb/pwiz-mzdb) file.
Thanks to mzDB and D language, extraction is much faster (20 fold) than the same algorithm implemented in C++ in OpenMS working on mzML (the PSI standard). 

### Compilation
Compilation with dmd

  `dmd -O -inline -release -noboundscheck extractor.d model.d`
