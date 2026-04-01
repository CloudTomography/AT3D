## OPAC Aerosols

I have imported the [OPAC](1) aerosols, both optical properties and profiles from libRadtran to here. I opted to convert the optical properties rather than reprocess them. The conversion was done using the script `convert_libradtran_phase_to_at3d.py`.

The conversion is necessary because libRadtran's phase functions are stored in legendre expansions for each polarization component rather than in the Wigner function expansions used by SHDOM.

The tables available here only include the spherical aerosols so the mineral species shouldn't really be trusted. Converting those tables would require the math for the conversion of the extra phase function components to Wigner coefficients.

The `aerosol/size_distr.cfg` file contains the distribution parameters of the aerosol size distributions. The original libRadtran data obtained [here](2) also contains the index of refraction data.

The standard OPAC aerosol mixture profiles used by libRadtran have fixed mass density ratios. This is different to the original OPAC profiles as it does not maintain aerosol composition as humidity varies for those aerosol species that humidify. As such, in AT3D, the aerosol module rescales these profiles so that they have fixed particle number concentration mixing ratios, similar to the original OPAC aerosols.

If you have any questions, please contact me.

 -Jesse Loveridge (jesserl2@illinois.edu)

 [1]: https://doi.org/10.1175/1520-0477(1998)079<0831:OPOAAC>2.0.CO;2
 [2]: http://www.libradtran.org/doku.php?id=download
