## This code is an update of 6S

Users are assumed to be familiar with 6S radiative transfer model. http://6s.ltdri.org/

The 6S code is rewritten in free form style. I try to keep original source code as much as possible but add necessary declaration of variables and modify all float (real) type variables to double precision type. Three surface BRDF models are added, which are PROSAIL (version 5B), ACRM(version 12.2012), and a snow reflectance model. I modify these code for convenience to be coupled in 6S.

Currently, the latest version of 6S is 6SV2.1.

Several bugs are fixed,

1. subroutine _akalbe_. In _AKTOOL.f90_ line 1656 (_AKTOOL.f_ line 1545), the dimensions of variables _uu_ and _aa_ are changed to 48 to keep consistent with the parameters in subroutine _dakg_.

2. _mainlutaero.f90_ line 2466 (_mainlutaero.f_ line 2508), _pveg = ul_ -> _pveg = uli_

3. subroutine _iso_ In _ISO.f90_ line 38, add the initialization of _xmus_, _xmus = rm(0)_

4. subroutine _disom_. In _DISCOM.f90_ line 60 (_DISCOM.f_ line 64), modify two if statements to avoid memory problems when compiling by Intel Visual Fortran on a Windows platform.

### Attention
There are two main functions, main.f90 and mainlutaero.f90. The Makefile in sixs90 folder is recommend, or just compile them separately in IDE.

Although the code is in free form style, many statements are written using old syntax. You can help improve and extend this code.
	
Aug 30, 2016.