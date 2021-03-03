# GPROPT
Generalized PROgram for OPTimization

GPROPT is a program coded in Fortran intended to provide a generic tool to minimize a user-defined objective function.
Several classical methods are implemented along with a genetic algorithm.

This program has been coded initially in 2015, in the framework of an assignment for a doctoral program.
I have decided now to make it available under the MIT license.
My expectation is that it becomes useful to anyone else.
Contributions to improve it are also welcomed.
A pull-request system is required to merge for the master branch.

## Documentation
The report of the assignment has been included in the repository: `doc/Igor_OPT.pdf`.
There you can find some theoretical background of the methods implemented, an overview of the program, and application examples.

## Compilation
A `Makefile` has been included to compile the code with Intel Fortran compiler `ifort`.
It has been tested in `Ubuntu OS`.
It is ready to compile with other compilers and in other OS, but some adjustments will be needed.
I hope to fix these issues in following versions.

You just need to open a terminal in the directory `src/` and type `make`.
To compile it in release mode (optimized) just type `make optm`.
Objects are stored in a directory `obj/` and the executables will appear in `bin/`.

Alternatively you can use any IDE suitable for Fortran (e.g. create a Fortran project in Visual Studio).

## Usage
When running the executable `GPROPT(_debug)`, the number of design variables and the optimization method are prompted.
For uni-dimensional problems the following options are available:
```
Uni-dimensional problem
 Select Optimization algorithm:
 1 - Polynomial approximation
 2 - Golden Section Method
 3 - Genetic Algorithm
 0 - QUIT

```
For a general multi-dimensional problem, the following prompt is displayed:
```
Select Optimization algorithm:
 1 - Powell Method
 2 - Steepest Descent
 3 - Fletcher-Reeves Conjugate Direction
 4 - Polak-Ribiere Conjugate Direction
 5 - Quasi-Newton Methods (Variable Metrics)
 6 - Newton Method
 7 - Genetic Algorithm
 0 - QUIT

```
Additional options will be displayed according to the method selected.

The output is written for the terminal a an results file `<method>.res`.

### User-defined problem
The user must define an appropriate interface to evaluate the objective function value in `src/evalfunc.f90` 
(or `src/evalfunc1.f90` for the golden section method).
Some examples of objective functions are available in this routine.

#### Examples
The results of the examples addressed in the assignment report are available in `examples/`.

## Contacts

Send me an e-mail to ilopes@fe.up.pt if you have any question or comment.
