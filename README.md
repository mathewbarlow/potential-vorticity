# potential-vorticity

Python code (<b>gfs_pv_1.2.py</b>) for dynamic tropopause (DT) calculations: DT pressure, DT potential temperature (theta), PV on 
the 330K isentropic surface, and a PV and theta cross-section at the latitude where the tropopause is lowest in the 
domain. The date and time need to be set in the beginning of the code; the domain can be changed there as well.

<b>gfs_pv_1.2_3D.py</b> is the same code but with additional (poor) 3D plots

This code has a DOI and is citable:  <a href="https://zenodo.org/badge/latestdoi/110735652"><img src="https://zenodo.org/badge/110735652.svg" alt="DOI"></a>

The data source is the online GFS analysis, and the date and time need to be set within the period of available data.

The program can take a few minutes to run because it is accessing data over the internet.

Samples of the four figures that are generated are included as PNG files:

<img align="left" width="300" height="300" src="image2_DT_theta.png">
<img align="left" width="300" height="300" src="image1_DT_p.png">
<img align="left" width="300" height="300" src="image3_pv330.png">
<img align="left" width="300" height="300" src="image4_pv_xsec.png">
