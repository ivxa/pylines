# pylines

This code computes the **streamlines** of **2D Relativistic Hydrodynamical simulations** (RHD). The expected input file is a Fortran 90 binary file containing the physical variables from the hydrodynamical simulation (e.g., the output from Ratpenat hydrodynamical code). We also provide the **branch `pylines-txt`** which accepts an arbitrary text file as an input file. The output files are the plots and the text files, one file per current line containing the physical variables along the line.

Usage:

1. Edit `read_fortran.py` (or `read_dat.py` for `pylines-txt`) according to the output and set-up of your hydrodynamical simulation (e.g., injector's location).
2. Edit `setup.py`:
    * `input_file`: Fortran binary file (or text file for `pylines-txt`) containing the physical output of the hydrodynamical simulation
    * `nick_name`: Nick name of the simulation
    * `output_file`: File name of the output file
    * `test_rhd`: Enable/Disable some testing (requieres your own test definition)
    * `tstep`: Time step (1 for maximum precision, >1 values for faster computations)
    * `itemax`: Maximum number of time steps for each current line calculation
    * `resamp`: Resample all the current lines to a fixed number of cells for each current line
    * `nlines`: Maximum number of current lines (the transversal surface is conserved)
    * `CGS_units`: Enable/Disable CGS unit conversion (only available in the master branch)
    * `fb0`: Scaling factor related to the magnetic field estimation
    * `tr0`: Above this tracer value the line is discarded (discard lines with high values of mixing)
    * `int_method`: Interpolation method (0: no interpolation, 1: bilinear, 2: one dimensional interpolation)
    * `int_test`: Enable/Disable interpolation testing
    * `make_plots`: Enable/Disable all the plots
    * `plots_path`: Plots path
    * `plot_maps`: Enable/Disable maps plots
    * `plot_lines`: Enable/Disable lines plots
    * `plot_profiles`: Enable/Disable profile plots
3. Execute: `python pylines.py`

This software has been developed during my PhD thesis (University of Barcelona) to
analyze RHD simulations of some amazing astrophysical objects, such as gamma-ray binaries, AGNs, and microqusars. The most part of the code has been written during my research stay in ISAS/JAXA (Tokyo) with the help and supervision of Dmitry Khangulyan. The hydrodynamical simulations used as an input for this code have been carried on in collaboration with Manel Perucho. All of this work have been done in collaboration with my supervisor Valentí Bosch-Ramon. Marc Ribó, my second supervisor, has also collaborated providing the observational connection.


## Author
Xavier Paredes-Fortuny (xparedesfortuny@gmail.com)


## License
See the [LICENSE](LICENSE.md) file for license rights and limitations (MIT).
