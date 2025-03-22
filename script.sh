#!/bin/bash

gfortran constants.f90 setup.f90 isotope_io.f90 integrator_mod.f90 integrator.f90 -o integrator
./integrator

