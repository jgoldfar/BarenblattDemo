# README #

[![Build Status](https://travis-ci.org/jgoldfar/BarenblattDemo.svg?branch=master)](https://travis-ci.org/jgoldfar/BarenblattDemo)

This repository contains a very basic implementation of a finite difference scheme for nonlinear diffusion-type equations based on methods given in A. A. Samarskii's "The Theory of Difference Schemes".
Results from this scheme match the qualitative picture (modulo noise due to approximation error) which we are deriving through theoretical work, and have been presented as poster presentations:

* [U. G. Abdulla and J. Goldfarb (presenter), *Numerical Analysis of Interface Evolution for the Nonlinear Degenerate Diffusion-Convection Equation*, MAA MathFest (2011)](http://www.maa.org/mathfest-2011-program)

* [U. G. Abdulla and J. Goldfarb (presenter) and N. Mertins, *Analysis of Interfaces for Nonlinear Diffusion-Convection Equations*, SIAM AN12 (2012)](http://meetings.siam.org/sess/dsp_programsess.cfm?SESSIONCODE=15153)

Work on the NDCE is ongoing, but this particular implementation has been superseded by other codes that I hope to make available some other time.
In particular, though, this repository exists to make the results of those posters reproducible, and to explore options for modifying and visualizing the output of these codes.

These codes are considered to be version 0.1.
See below for the set-up instructions, such as they are, and the Markdown files under `doc` for more details on the codes themselves.

### How do I get set up? ###

* Cloning this repository and should be enough to bring this package into your Python current scope to begin working on it.

* Dependencies are minimal; the code is tested on Python 2.7, 3.4, and 3.6. PyPy{2,3} should be supported as long as plotting (from Python) is not required.
Running

```shell
pip install -r requirements.txt
```
should install a known-good version of `NumPy`, as well as the plotting library `matplotlib` and the test framework `nose`.
All other necessary implementations are included.

### Contribution guidelines ###

Contribution is encouraged: unit tests (written in nose) are particularly in need.
All code is reviewed by the project maintainer.
Improvements and changes can be discussed through a pull request and/or corresponding issue.

### Roadmap

* v1.0: Release current version of codes with full documentation and "modern", portable CLI utilities.

* v1.1: Modernize/refactor codes into more idiomatic NumPy.

* v1.2: Refactor performance-critical codes into C or Cython.

### Who do I talk to? ###

* Jonathan Goldfarb <jgoldfar@gmail.com>
