star-formation-histories
========================

Leitner2012 star formation history code.

This repository exists to support use of the star formation histories compiled in 
http://adsabs.harvard.edu/abs/2012ApJ...745..149L 

A somewhat cleaned up and simplified version of the code I used is in simplemsi.tar
Please let me know if you have any issues at all. The readme file explains how to compile that code, what the required options are, and what the output file contains.

Star formation histories are also directly available in sfrs_snl.tar and ms_snl.tar (see the enclosed readme). 

NOTE: Table 1 headings are not quite right. Specifically, A11 should have been 10^11Msun times the values listed in the table (psi is a SFR not an SSFR). 

The following (supermongo) may help clarify as well: sfr_anfitR below is the SFR for a galaxy given a constant recycling rate and the Karim data (in supermongo syntax). Adjust Mgal0 for different z=0 masses; msi_z is the redshift array you choose to use; A0 (in Msun/yr) is A11 in yr^-1 in the paper times 10^11Msun.

<pre><code>
set Mgal0=10**10
set Rfact=0.45
set beta=-0.35
set alpha=3.45
set M0=10**11.0
set dadt=0.064/1e9
set A0= 3.24
set r11=(Mgal0/M0)**(1+beta)
set m_anfit = M0*( ((Mgal0/M0)**-beta +beta/dadt*(1-Rfact)/(alpha-1)*A0/M0*((1+msi_z)**(alpha-1)-1)))**(-1/beta)
set sfr_anfitR = A0*( m_anfit/M0 )**(1+beta)*(1+msi_z)**alpha
</code></pre>

