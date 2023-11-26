# HEFT Code
This file is intended to provide an overview of how I would use this code to perform a typical HEFT analysis.
I will provide a brief overview of my typical workflow, with more detail provided for each step afterwards.

# Work Flow
1. Run 'newProject.py' to easily set up all the config files, and skip to step 4.
To do the initial setup manually, go to step 2.
2. Make a folder for the new project, and link the makefile from 'heftCode/makefile' in this folder.
3. Enter the details of the system into 'HEFT.config'.
4. Set up an initial guess for the fit parameters in \ttt{allFits.params}.
   If I'm trying to replicate the results of a paper, I will then get the fit parameters from that paper and put them into \ttt{allFits.params}, then skip to step TODO.
5. Get some scattering data to fit to, generally the WI08 solution from [SAID](https://gwdac.phys.gwu.edu/).
   This should be in a file called \ttt{dataInf.in} (see the example in this folder)
6. Setup \ttt{HEFTFitting.config} with the bounds on each category of fit parameter, and which parameters will be active in the fitting process.
7. If I want to run a bounded fit, make and run \ttt{fitBQ.x}.
   This program uses Powell's BOBYQA algorithm to minimise a $\chi^2$ between the HEFT phase shift and inelasticity, and the scattering data.
   If I don't want to use bounds, make and run \ttt{fit.x}, which uses minfun (Powell's NEWUOA algorithm).
8. Make sure \ttt{HEFTInfinite.config} is correct, with the correct energy range.
9. Set \ttt{iParamChoice} in \ttt{allFits.params} to the new parameter set found by the fit procedure, then make and run \ttt{inf.x} to calculate the scattering observables for this parameter set.
   I usually plot these using a file called \ttt{PhasevE.py}, which is copied from another project I've worked on.
10. To search for pole positions, you can make and run \ttt{poles.x}.
	This does a grid search from $ 1.0 - 0.01i $ to around $ 2.0 - 0.2i $ GeV, and should pick up the positions of any poles in the $ T $-matrix.
11. In a finite-volume, we're typically interested in how the eigenvalues and eigenvectors of the finite Hamiltonian vary with $m_{\pi}^{2}$.
	Finite-volume config options such as the $m_{\pi}^2$ range to be calculated are found in \ttt{HEFTFinite.config}.
12. To find the slope of the bare mass(es), I need lattice data to fit to.
	I store this in a file called \ttt{slopeFitMasses.data}.
13. To fit the bare mass slope(s), I make and run \ttt{fitBare.x}.
	The config file \ttt{HEFTFinite.config} contains the guess for the bare slope, and the order of the bare mass extrapolation.
14. To calculate the finite-volume energy spectrum, there are two options.
	The usual case is that the lattice size $L$ is varying with $m_{\pi}^{2}$, so I run \ttt{lqcdFin.x}.
	If I want to keep $L$ constant and just vary $m_{\pi}^{2}$, I run \ttt{mpiFin.x}.
15. To plot the spectrum, I have a plot script called \ttt{EvMpi.py}.
	Typically I'll have a few of these for different sets of lattice data, such as in the odd-parity nucleon case where I have \ttt{EvMpi\_3fm.py}, \ttt{EvMpi\_2fm.py}, and \ttt{EvMpi\_D200.py}.
16. To plot the eigenvectors, I have a script called \ttt{EvecvMpi.py}.
17. **Extra:** There are also a few other things I might be interested in, such as comparing multiple parameter sets, such as in the $\Delta(1232)$ regularisation and multiple bare state papers.
  I can search for poles with a varying set of parameters with \ttt{multiFitPoles.x}.
  I calculate the finite-volume spectrum for a varying set of parameters with \ttt{multiFitFin.x}.

## 1. newProject.py
todo

## 3. HEFT.config
here.
'   1  # NStar_2b3c config
   2
   3  n_channels   3
   4  n_bare       2
   5
   6  # Channel                   partialWave
   7  pion          Nucleon       S
   8  eta           Nucleon       S
   9  Kaon          Lambda        S
  10
  11  OnshellChannel    1
  12  useCustomMasses   F
  13  particles_3fm.in
  14
  15  # Bare state labels
  16  N1
  17  N2
  18
  19  # g (2-1 potential), v (2-2),   u (regulator)
  20  B                    B          A'
end
