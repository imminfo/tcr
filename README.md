[![Licence](https://img.shields.io/hexpm/l/plug.svg)](http://www.apache.org/licenses/LICENSE-2.0)
[![RepSeq](http://statsarray.com/wp-content/uploads/2014/03/omictools-logo.png)](http://omictools.com/tcr-s9043.html)


tcR
===

tcR is a platform designed for TCR repertoire data analysis in R after preprocessing data with tools for CDR3 extraction and gene segments aligning (MiTCR, MiGEC, IgBLAST, Decombinator, etc.). With the power and flexibility of R language and procedures supported by tcR users can perform advanced statistical analysis of TCR repertoires.

See tcR website for more information, manual and examples: [http://imminfo.github.io/tcr/](http://imminfo.github.io/tcr/)

If you have any questions, suggestions or bug reports, feel free to raise an issue here: [https://github.com/imminfo/tcr/issues](https://github.com/imminfo/tcr/issues)

*Warning!*
tcR internally expects columns with nucleotide and amino acid CDR3 sequences and columns with gene segments to have character class, not factor class. Use `stringsAsFactors=FALSE` parameter if you use R functions for parsing files with tables (.csv, .xls and others).