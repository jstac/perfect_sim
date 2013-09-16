
.. _perfect_sim:

******************************************************************************
Perfect Simulation of Stationary Equilibria
******************************************************************************

This page collects files and computer code for the paper **Perfect Simulation of Stationary Equilibria**
by Kazuo Nishimura and John Stachurski.

Publication Details
-----------------------

| Perfect Simulation of Stationary Equilibria
| Kazuo Nishimura and John Stachurski
| **Journal of Economic Dynamics and Control**, 34, 577--584, 2010


Abstract
----------


Using a variation of the coupling from the past technique, this paper
develops algorithms which generate independent observations from
the stationary distributions of various dynamic economic models. These variates
can be used for calibration, calculation of steady state phenomena, and
simulation-based estimation.  As an application, we demonstrate how to
generate exact samples from the stationary distribution of an incomplete
markets model routinely calibrated by macroeconomists.  Our
implementation generates 100,000 independent draws from the stationary
distribution in less than 3 seconds.



Code
--------

C code for replicating results in the paper can be found in the repository.
C or Fortran are significantly faster than high-level languages with
this algorithm.  However, for illustrative purposes, I've also included a
Python implementation below.  There are two files.  The first computes the
policy function and the second implements perfect sampling.
