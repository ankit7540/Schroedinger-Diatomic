# Usage

Load all the procedures in IgorPro first, then press compile. After this run `ini_hm(); ini_alpha()` in the command browser. This will open a total of three panels showing various buttons for defining the input and running the calculation.

![panels](./../img/panels.png)

 - **Panel (a)**: prepare H section is for H2 and isotopologues.
The other sections (normalization, compute expectation values and performing these in batches) are general in nature and can be used for other systems as well.

 - **Panel (b)**: Prepare_H_matrix is for general diatomic system, for a given potential and inter-nuclear distance available as 1D waves.
The other sections are for viewing the results after computation, i.e. print out the energy levels for various transitions, etc..

 - **Panel ( c)**: For processing polarizability data available as 2D over wavelength and inter-nuclear distance. Functions for the interpolation of polarizability over specific wavelengths are included (based on cublic-spline interpolation). Functions for the computation of expectation values are included, to perform the computation via numerical integration. For the evaluaton of expectation values of matrix elements, i.e. computation of ro-vibrational matrix elements of polarizability invariants for specific wavelength, refer to faster implementations in python/FORTRAN used in the repository. (https://github.com/ankit7540/H2-PolarizabilityMatrixElements)

This documentation is incomplete. Use Issues to to raise questions  and/or to request further documentation.
