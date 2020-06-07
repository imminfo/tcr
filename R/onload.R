.onAttach <- function(libname = find.package("tcR"), pkgname = "tcR"){
  packageStartupMessage("
========================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The tcR package WILL SOON BE ORPHANED 
!! AND REMOVED FROM CRAN.
!!
!! A new package is available that is 
!! designed to replace tcR: 
!! immunarch  --  https://immunarch.com/
!!
!! We will be happy to help you to move
!! to the new package. Feel free to contact us:
!! http://github.com/immunomind/immunarch
!!
!! Sincerely, 
!!  immunarch dev team and 
!!  Vadim I. Nazarov, lead developer of tcR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
=======================================")
}