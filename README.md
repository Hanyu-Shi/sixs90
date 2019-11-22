## This code is an update of 6S

* Users are assumed to be familiar with 6S radiative transfer model. http://6s.ltdri.org/

* The 6S code is rewritten in free form style. I try to keep original source code as much as possible but add necessary declaration of variables and modify all float (real) type variables to double precision type. Three surface BRDF models are added, which are PROSAIL (version 5B), ACRM(version 12.2012), and a snow reflectance model. I modify these code for convenience to be coupled in 6S.

* Currently, the latest version of 6S is 6SV2.1.

Several bugs are fixed,

1. subroutine _akalbe_. In _AKTOOL.f90_ line 1656 (_AKTOOL.f_ line 1545), the dimensions of variables _uu_ and _aa_ are changed to 48 to keep consistent with the parameters in subroutine _dakg_.

2. _mainlutaero.f90_ line 2466 (_mainlutaero.f_ line 2508), _pveg = ul_ -> _pveg = uli_

3. subroutine _iso_ In _ISO.f90_ line 38, add the initialization of _xmus_, _xmus = rm(0)_

4. subroutine _disom_. In _DISCOM.f90_ line 60 (_DISCOM.f_ line 64), modify two if statements to avoid memory problems when compiling by Intel Visual Fortran on a Windows platform.

**Usage**

* There are two main functions, main.f90 and mainlutaero.f90. The Makefile in sixs90 folder is recommend, or just compile them separately in IDE.

**Attention**
* Although the code is in free form style, many statements are written using old syntax. You can help improve and extend this code.
	
Aug 30, 2016.

---

**Update**

* See sixs90v15/ for some updates

  May 24, 2017


* This work was present at the 2019 IGARSS conference. 

  [Updates of the 6S radiative transfer model: a case study of 6S+PROSAIL](https://ieeexplore.ieee.org/document/8899146/).

  Nov. 22, 2019

**Reference**

[1] H. Shi and Z. Xiao, “Updates of the 6S radiative transfer model: a case study of 6S+PROSAIL,” in 2019 IEEE International Geoscience and Remote Sensing Symposium (IGARSS), Yokohama, Japan, 2019, pp. 2879–2882. [doi:10.1109/IGARSS.2019.8899146](https://ieeexplore.ieee.org/document/8899146/)

[2] E. F. Vermote, D. Tanre, J. L. Deuze, M. Herman, and J.-J. Morcette, “Second Simulation of the Satellite Signal in the Solar Spectrum, 6S: an overview,” IEEE Trans. Geosci. Remote Sens., vol. 35, no. 3, pp. 675–686, May 1997. [doi:10.1109/36.581987](https://doi.org/10.1109/36.581987)

[3] S. Y. Kotchenova, E. F. Vermote, R. Matarrese, and F. J. Klemm, Jr., “Validation of a vector version of the 6S radiative transfer code for atmospheric correction of satellite data. Part I: Path radiance,” Appl. Opt., vol. 45, no. 26, pp. 6762–6774, 2006. [doi:10.1364/AO.45.006762](https://doi.org/10.1364/AO.45.006762)

[4] S. Y. Kotchenova and E. F. Vermote, “Validation of a vector version of the 6S radiative transfer code for atmospheric correction of satellite data. Part II. Homogeneous Lambertian and anisotropic surfaces,” Appl. Opt., vol. 46, no. 20, pp. 4455–4464, 2007. [doi:10.1364/AO.46.004455](https://doi.org/10.1364/AO.46.004455)

[5] S. Y. Kotchenova, E. F. Vermote, R. Levy, and A. Lyapustin, “Radiative transfer codes for atmospheric correction and aerosol retrieval: intercomparison study,” Appl. Opt., vol. 47, no. 13, pp. 2215–2226, 2008. [doi:10.1364/AO.47.002215](https://doi.org/10.1364/AO.47.002215)

