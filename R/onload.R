.onAttach <- function(libname = find.package("tcR"), pkgname = "tcR"){
  packageStartupMessage("
==================================
==================================
The tcR package is no longer supported and current issues will not be fixed. A new package is available that is designed to replace tcR called immunarch. 
We have solved most of the problems tcR package had and improved the overall pipeline, providing functions for painless repertoire file parsing and publication-ready plot making. 

The mission of immunarch is to make immune repertoire data analysis as easy and possible - even with R. 
Please feel free to check it here: 
https://immunarch.com/

We will be happy to help you to integrate the new package into your pipelines. Please do not hesitate to contact us, should any question arise:
Email: vdm.nazarov at gmail.com
LinkedIn: https://linkedin.com/in/vdnaz

Sincerely, 
  immunarch dev team and Vadim I. Nazarov, lead developer
        
P.S. To suppress this message, just wrap any call to tcR as follows:
suppressPackageStartupMessages(library(tcR))
==================================
==================================")
}