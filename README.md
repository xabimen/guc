---
title: Get Unit Cell
description: Program that finds smaller unit cells.
author: Xabier M. Aretxabaleta
created:  2017 Oct 19
Last update: 2017 Oct 19

---


# Get Unit Cell

This program finds a smaller unit cell from a structure in a supercell. It is very usefull to use it with structures predicted with evolutionary algorithms.

Download
--------

To download create a directory and use the following command in terminal:

	$ git clone https://github.com/xabimen/guc

Compiling
---------

To compile the program there is a script called 'compile_guc.sh'. First, give the script executing permissions with the following command.

	$ chmod -x compile_guc.sh

Finally, execute the script by

	$ ./compile_guc.sh

The executable file 'guc' will be created

Executing
---------

To execute the program use the following command:

	$ ./guc [input file name] [number of different types of atoms] [tolerance]
	
#### *Report Bugs to xabierm.mendez@ehu.eus*
