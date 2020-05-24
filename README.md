# ronchigram-matlab
By Noah Schnitzer, Suk Hyun Sung @ [Hovden Lab](http://hovdenlab.com/)

A MATLAB library for quick simulation of the STEM probe and electron Ronchigram, and training of a CNN to assess probe quality from the Ronchigram.

See also: http://ronchigram.com/ ([source](https://github.com/sukhsung/ronchigram/)) for a similar JS/Wasm project.

MATLAB scripts and functions are organized into folders but are mutually dependent -- all subfolders must be added to the path to work properly. `assessment` and `simulation` functions are well tested for reasonable inputs, `dataset_building` scripts are not and should be used cautiously/as examples.

To get started, check out `misc/example.m`


## Description of scripts and functions:

- `assessment`: functions to calculate heuristics for aberration functions
	- `indiv_p4_calculator.m`: Calculate the individual aberration phase shift convergence angle for a single aberration
	- `par_strehl_calculator.m`: Calculate the Strehl ratio converence angle (parallel) for aberrations
	- `pi4_calculator.m`: Calculate the total aberration phase shift convergence angle for a single aberration
	- `probe_sizer.m`: Wraps `resolution_test.m` to calculate 50% probe current diameter for given aberration functions and convergence angles
	- `resolution_test.m`: Poorly named probe size assessment, should be used with `effprobe` option
	- `strehl_calculator.m`: Calculate the Strehl ratio converence angle for aberrations
- `CNN`: scripts to train and test CNN on simulated Ronchigrams
	- `lim_net.m`: Train CNN used in paper
	- `transfer_learning.m`: Transfer learn on Alexnet
- `dataset_building`: scripts and functions to build data sets. Note many parameters are embedded in scripts, have only been tested for limitied inputs.
	- `aberration_series.m`: Calculates heuristics as specific aberrations are varied. Used to e.g. find defocus to balance other aberrations.
	- `dataset_generator.m`: Generates a dataset with a ~ uniform distibution of convergence angles. 
	- `defocus_distribution.m`: More slowly generates a dataset with ~ uniform distribution of convergence angles and defocus set to compensate other aberrations.
	- `distribution_generator.m`: Subroutine to generate aberrations and dynamically scale to try to get uniform CA distribtion.
- `misc`: miscellaneous scripts and functions, e.g. utilities and examples
	- `colordef.m`: Colors from paper
	- `example.m`: Brief walkthrough for Ronchigram simulation and heuristic calculation
	- `get_aberration.m`: Calculates the phase shift for a specific aberration
	- `normalize_data.m`: Normalizes data to min 0 max 1
	- `px_to_ang.m `: Calculates the scale factor Å/px for a given accelerating voltage and simulation dimension
	- `radial_average.m`: Calculates a radial average
- `simulation`: functions for STEM probe and Ronchigram simulation
	- `aberration_generator.m`: Generates aberrations with random magnitude and angle out to 5th order. Relative scale of magnitudes is based off observations from AC-STEM, Kirkland 2011 Ultramicrosc.
	- `aperture_mask.m`: Generates a binary disc of given radius
	- `calculate_aberration_function.m`: Calculates phase shift for an aberration function
	- `calculate_probe.m`: Calculates STEM probe given ab phase shift, convergence angle
	- `shifted_ronchigram.m`: Calculates a Ronchigram for an aberration function, shift in Ronch center relative to obj aperture, and convergence angle
	- `shifted_ronchigram_o.m`: Legacy shifted_ronchigram to match parameters (e.g. grating resize factor) of trained networks.