#################################################################
# TestSuite for Drops						#
#################################################################
Timo Henrich , SC RWTH Aachen
Version 0.1
12.09.2007

This application is used to run and test application ,written for Drops, under
varying configuration an number of used processors.
application that are meant to be used with TestSuite should take their configuration 
from a seperate config-file and write their results (also error-messages) straightly to screen.

## Installation ##
TestSuite is installed using a makefile. A simple 
	make all 
in the applications directory should be enough.

Note: If you don't use 'mpirun' to run your MPI-applications you have to change the "mprunner"-attribute in your configuration-file.

## Configuration ##
All tests that should be performed are defined in a special formed xml-configuration
file.
Look config.example for more information.

## Usage ##
After you have written a configuration-file you can start using TestSuite.
Simply type:
	TestSuite [config-file] ([output-file])
Where config-file is the mentioned xml-file and output-file an optional parameter.
If output-file is set, all results a written in HTML-format to the specified location.

