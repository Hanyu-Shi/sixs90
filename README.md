## This code is an update of 6S

* Users are assumed to be familiar with the 6S radiative transfer model. https://salsa.umd.edu/6spage.html

* The 6S code is rewritten in the free-form style. I try to keep the original source code as much as possible but add necessary declarations of variables and modify all float (real) type variables to double precision type. Three surface BRDF models are added: PROSAIL-D, ACRM (version 12.2012), and a snow reflectance model. I modified their original code for convenience to be coupled in 6S.

* Currently, the latest version of 6S is 6SV2.1(2014)

* You may be interested in exploring [sixs90v15](sixs90v15/) if encountering problems


**Usage**

* There are two main functions, main.f90 and mainlutaero.f90. The Makefile in sixs90 folder is recommended, or just compile them separately in IDE.


**Attention**

* Although the code is in the free-form style, many statements are written using old syntax. You can help improve and extend this code.


**Update**

* Roll back to 6SV2.1(2014)

  I noticed that the official 6S code rolled back to the previous version, and I followed their footsteps here. The update on Feb. 11, 2020 is invalid now.
  
  Feb. 26, 2022


* Update to 6S-V2.1 (November, 2018)
   
   * Solar irradiance spectrum is updated
   * A "truncation method" is adopted to avoid negative discriminant. Compare main.f90 and main_v0.f90 for details

  Feb. 11, 2020
  

* This work was present at the 2019 IGARSS conference. 

  [Updates of the 6S radiative transfer model: a case study of 6S+PROSAIL](https://ieeexplore.ieee.org/document/8899146/).

  Nov. 22, 2019
 

* See sixs90v15/ for some updates

  May 24, 2017


* Several bugs are fixed,

  1. subroutine _akalbe_. In _AKTOOL.f90_ line 1656 (_AKTOOL.f_ line 1545), the dimensions of variables _uu_ and _aa_ are changed to 48 to keep consistent with the parameters in subroutine _dakg_.

  2. _mainlutaero.f90_ line 2466 (_mainlutaero.f_ line 2508), _pveg = ul_ -> _pveg = uli_

  3. subroutine _iso_ In _ISO.f90_ line 38, add the initialization of _xmus_, _xmus = rm(0)_

  4. subroutine _disom_. In _DISCOM.f90_ line 60 (_DISCOM.f_ line 64), modify two if statements to avoid memory problems when compiling by Intel Visual Fortran on a Windows platform.

  Aug 30, 2016.


**Reference**

[1] H. Shi and Z. Xiao, “Updates of the 6S radiative transfer model: a case study of 6S+PROSAIL,” in 2019 IEEE International Geoscience and Remote Sensing Symposium (IGARSS), Yokohama, Japan, 2019, pp. 2879–2882. [doi:10.1109/IGARSS.2019.8899146](https://ieeexplore.ieee.org/document/8899146/)

[2] E. F. Vermote, D. Tanre, J. L. Deuze, M. Herman, and J.-J. Morcette, “Second Simulation of the Satellite Signal in the Solar Spectrum, 6S: an overview,” IEEE Trans. Geosci. Remote Sens., vol. 35, no. 3, pp. 675–686, May 1997. [doi:10.1109/36.581987](https://doi.org/10.1109/36.581987)

[3] S. Y. Kotchenova, E. F. Vermote, R. Matarrese, and F. J. Klemm, Jr., “Validation of a vector version of the 6S radiative transfer code for atmospheric correction of satellite data. Part I: Path radiance,” Appl. Opt., vol. 45, no. 26, pp. 6762–6774, 2006. [doi:10.1364/AO.45.006762](https://doi.org/10.1364/AO.45.006762)

[4] S. Y. Kotchenova and E. F. Vermote, “Validation of a vector version of the 6S radiative transfer code for atmospheric correction of satellite data. Part II. Homogeneous Lambertian and anisotropic surfaces,” Appl. Opt., vol. 46, no. 20, pp. 4455–4464, 2007. [doi:10.1364/AO.46.004455](https://doi.org/10.1364/AO.46.004455)

[5] S. Y. Kotchenova, E. F. Vermote, R. Levy, and A. Lyapustin, “Radiative transfer codes for atmospheric correction and aerosol retrieval: intercomparison study,” Appl. Opt., vol. 47, no. 13, pp. 2215–2226, 2008. [doi:10.1364/AO.47.002215](https://doi.org/10.1364/AO.47.002215)

