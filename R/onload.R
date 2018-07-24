.onAttach <- function(libname = find.package("tcR"), pkgname = "tcR"){
  packageStartupMessage("
==================================
==================================
We want to help you improve the quality of the research and make a new package with more features added that suit your needs.
Please, complete a survey on immune repertoire analysis and help us understand what features should we definitely implement.
Link to the survey:
https://goo.gl/forms/RcLqJwkUPNfKfaOw2

Thank you so much for your time! Please do not hesitate to contact us, should any question arise:
Email: vdm.nazarov@gmail.com
LinkedIn: https://linkedin.com/in/vdnaz
  
Sincerely, 
  tcR 3.0 Dev Team
        
P.S. To suppress this message, just wrap any call to tcR as follows:
suppressPackageStartupMessages(library(tcR))
==================================
==================================")
}