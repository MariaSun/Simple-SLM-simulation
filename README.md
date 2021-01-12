# Simple-SLM-simulation

This code provides simplistic representation of the Fraunhofer diffraction on a DXD screen.

Function propFF(u1,L1,lambda_0,z) does far-field diffraction propagation on a square screen u1 of dimentions L1=D illuminated by light of the wavelength lambda_0 (credit to Jonathan George for converting this function from Matlab to python, see Computational Fourier optics by D. Voelz).

Function screen(amplitude, phase) returns a complex-valued DXD modulation matrix that corresponds to supplied amplitude and phase DXD arrays.

Current version performs a simulation for the diffraction propagation that generates an OAM mode with vorticity 10.
