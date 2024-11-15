# SNAPI - Spectropolarimetric Non-LTE Analytically Powered Inversion code

SNAPI is a code intended for modeling and interpreting spectropolarimetric observations in solar spectral lines. It is written in C++ and aimed at massively parallel calculations.

## Installation 

The code requires working c++ compiler (gcc is used for the development), autotools (autoconf) and installation or at least download of fft3w and cfitsio. The following assumes installation using debian like linux, and a gcc. The code works fine on WSL2 in Windows, and it should work at MacOS, too. We did not test it yet on M chips, so that might require some meddling. Do not hesitate to reach out. 

### Installing cfitsio 

Download cfitsio from https://heasarc.gsfc.nasa.gov/fitsio/

The installation instructions are given in the README file of cfitsio, but for most Unix-like systems unpacking, and running the following from the unpacked directory should work: 

./configure

make 

make install 

### Installing fft3w 

Similarly to the previous one, you will want to download it from: https://www.fftw.org/download.html

And do:

./configure

make

make install

### Installing snapi 

Clone snapi from this get repo: 

https://github.com/ivanzmilic/snapi

Then: 

cd snapi 

Try: 

autoconf 

Follow that up by: 

./configure 

Typical reason for this to fail is that the script cannot find cfitsio. As the code only really needs headers, try with: 

./configure --with-cfitsio=/path/where/you/unpacked/cfitsio

If that goes through fine, the only thing that remains to be done is: 

make clean

make (optionally with, say, -j8 to compile in parallel)

After that, you are ready to use the code. 

### Testing the installation 

Open three instances of terminal, or use screen command in terminal to keep things neat. Try the following:

1) Enter snapi directory in the first terminal. Run master thread: 

./master/imaster -v 

2) Enter snapi/cfg directory in the second terminal. Run

../jsub/jsuv -c -cfg synth.cfg 

If that goes through fine, do: 

tail -f invlog.00001 

3) Enter snapi/slave (apologies for terminology, will be changed) and run

./islave 

4) if you go back to the second terminal, the job should be finished within a fraction of a second

5) the spectrum is put into the falc_6300.dat ascii file in the /snapi directory. Inspect the structure of the file and plot the values of the second column (intensity), against the values of the first column (wavelength, in angstroms) to check if everything is fine. 

### Congrats, now you can start doing more fancy stuff!

### Modeling 

The code reas in a 4D cube in the pyana (.f0) format. The cube has the dimensions 12 x NX x NY x NZ, where 12 is the number of physical parameters that serve as input while NX x NY x NY are physical dimensions. 

The 12 parameters, by index are: 
0 - log optical depth at 500 nm
1 - height [cm]
2 - Temperature [K]
3 - Total gas pressure [erg / cm^3]
4 - Total electron pressure [erg / cm ^3] - not needed, will actually be re-calculated by the code
5 - Mass density - - not needed, will actually be re-calculated by the code
6 - Magnitude of the magnetic field 
7 - Microturbulent velocity 
8 - LOS velocity 
9 - Continuum opacity - not needed, will actually be re-calculated by the code
10 - Inclination of the magnetic field 
11 - Azimuth of the magnetic field 

A "mini" version of such an atmosphere will be implemented, where only one vertical coordinate, gas pressure, temperature and velocities are provided - for spectroscopic calculations

It can handle non-LTE spectral lines using formalism of Rybicki and Hummer (ALI, 1991).
