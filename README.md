# convection-in-fortran

This numerical modeling project simulates low Prandt number convections. It is coded in modern FORTRAN. 

## Makefile

The Makefile allows running the simulation and to launch the visualization program coded in Python.
- `make exec` : compiles and runs the executable file.
- `make plot` : compiles and runs the executable file. In addition, it runs the `visualization.py` program to generate a movie of the convection system.
- `make vis`  : only runs `visualization.py`. 
- `make clean` : cleans the folder from unecessary files

## Input

The input is a `.txt` files that contain all the necessary constants. Change `Pr = n` from 0.01 to 10.0 to observe different behaviors. 

## Output

The output of the FORTRAN program is a `.txt` files that compiles every selected time steps. In addition, a `.csv` file is generated and contains some information about the program.

## Visualization

The visualization is done using a Python program. From the text file, the program extracts every time steps into a matrix and plotted into a matplotlib contourf plot. Every plot is saved. Then the program compiles all the plots into a movie in format `.avi` (can be changed).
