[![CRAN](http://www.r-pkg.org/badges/version/tcR?style=flat-square)](https://cran.r-project.org/package=tcR)
[![Downloads_all](http://cranlogs.r-pkg.org/badges/grand-total/tcR)](http://www.r-pkg.org/pkg/tcR)
[![Downloads_week](http://cranlogs.r-pkg.org/badges/last-week/tcR)](http://www.r-pkg.org/pkg/tcR)
[![Licence](https://img.shields.io/hexpm/l/plug.svg?style=flat-square)](http://www.apache.org/licenses/LICENSE-2.0)

tcR
===

*The tcR package is no longer supported and current issues will not be fixed. A new package is available that is designed to replace tcR called immunarch.*
*We have solved most of the problems tcR package had and improved the overall pipeline, providing functions for painless repertoire file parsing and publication-ready plot making.*

*The mission of immunarch is to make immune repertoire data analysis as easy and possible - even with R.*
*Please feel free to check it here: https://immunarch.com/*

*We will be happy to help you to integrate the new package into your pipelines. Please do not hesitate to contact us via emails on https://immunarch.com/ or via issues on https://github.com/immunomind/immunarch, should any question arise.*

*Sincerely, immunarch dev team and Vadim I. Nazarov, lead developer*

tcR is a platform designed for TCR and Ig repertoire data analysis in R after preprocessing data with software tools for CDR3 extraction and gene segments aligning (MiTCR, MiXCR, MiGEC, ImmunoSEQ, IMSEQ, etc.). With the power and flexibility of R language and procedures supported by tcR users can perform advanced statistical analysis of TCR and Ig repertoires. The package was published in BMC Bioinformatics, please cite if you use it:

[Nazarov et al., tcR: an R package for T cell receptor repertoire advanced data analysis](http://www.biomedcentral.com/1471-2105/16/175)

See tcR website for more information, manual and examples: [http://imminfo.github.io/tcr/](http://imminfo.github.io/tcr/)

If you have any questions, suggestions or bug reports, feel free to raise an issue here: [https://github.com/imminfo/tcr/issues](https://github.com/imminfo/tcr/issues)

The project was developed mainly in the [Laboratory of Comparative and Functional Genomics](http://labcfg.ibch.ru/lcfg.html).

*Warning!*
tcR internally expects columns with nucleotide and amino acid CDR3 sequences and columns with gene segments to have character class, not factor class. Use `stringsAsFactors=FALSE` parameter if you use R functions for parsing files with tables (.csv, .xls and others).

*Note for installation on Macs with OSX Yosemite (and potentially other versions):  if you receive a compilation error, modify tcR/src/Makvars to:*

```
CXX=clang++
```
