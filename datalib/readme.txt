files for IERS convetions

files seem to be up to date as of apr 2021  dav


https://iers-conventions.obspm.fr/chapter5.php


IAU 2006/2000 CIO:
Table 5.2a - IAU 2006/2000A expression of the X coordinate of the CIP in the GCRS. [4]
             iau06xtab5.2.a.dat
             coefficients a0xs, a0xsi

Table 5.2b - IAU 2006/2000A expression of the Y coordinate of the CIP in the GCRS. [4]
             iau06ytab5.2.b.dat
             coefficients a0ys, a0ysi

Table 5.2d - IAU 2006/2000A expression of the quantity s(t) + XY/2.[4]
             iau06stab5.2.d.dat
             coefficients a0ss, a0ssi

Table 5.2e - IAU 2006/2000A expression for Greenwich Sidereal Time GST. [4]
             iau06gsttab5.2.e.dat
             coefficients a0gst, a0gsti

Difference
Table 5.2f - Expression for the difference between the IAU 2006 (or IAU 2000A_R06) and IAU 2000 
             expression for various quantities (X, Y, Dpsi, Deps, s+XY/2). [4]
             differences of IAU 2006/2000 and 2000a ? 

IAU 2006/2000 Classical approach
Table 5.3a - Nutation in longitude (IAU 2000_R06 expression) derived from the IAU 2000A lunisolar 
             and planetary components with slight IAU 2006 adjustments. [4]
             iau06nlontab5.3.a.dat
             coefficients a0nl, a0nli

Table 5.3b - Nutation in obliquity (IAU 2000_R06 expression) derived from the IAU 2000A lunisolar 
             and planetary components with slight IAU 2006 adjustments. [4]
             iau06nobltab5.3.b.dat
             coefficients a0npl, a0npli

These aren't used anymore?                          
             iau00an.dat   2000a version of nutation 678
             iau00apl.dat  2000a version of planetary long 687
             
Notes:
The 2000a nutation and planetary nutation only seem to be online embedded in the SOFA source code, and 
they're different from those originaly specified in IERS TN-32. these new coefficients are in iau00annew.dat 
and iau00aplnew.dat, but not used yet.
http://www.iausofa.org/2021_0125_C/sofa/nut00a.c                 
             
It "looks" like the 678 687 terms are still used by nut00a? and that they use the 2006/2000
values? Much better option to simply use the series aproach with interpolation!



                          