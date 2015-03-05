Fast and Efficient Black Box Optimization using the Parameter-less Population Pyramid
==

To compile you will need C++11.  We use gcc version 4.8.3 for our complilation.

Our build system uses Makefiles to build.  You can compile the release version
by changing directory to Release and calling "make P3".

All of the source code is available in the 'src' directory.

To run an experiment, call the executable with command line arguments for configuration.
This will run the default test configuration:

Release/P3 config/default.cfg

The command line accepts any number of configuration files and configuration flags given
in the form "-key value" where key is the name of the option to change and value is the
value to change it to.  Arguments override in last given order.  For example:

Release/P3 config/default.cfg config/p3.cfg -problem NearestNeighborNK -length 50

This line will use the default configuration, replacing configuration values with those
found in p3.cfg, setting the problem to Nearest Neighbor NK with a genome size of 50.
All experiments combined "default.cfg" "tune.cfg" and the solver specific configuration file.

Also, see main.cpp for implementation related details.

All of the processed output files are contained in data/csv_files.tar.gz.
All data analysis was performed using these files and the scripts found
in the data folder.
