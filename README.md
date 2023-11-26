# HEFT Code
This file is intended to provide an overview of how I would use this code to perform a typical HEFT analysis.
I will provide a brief overview of my typical workflow, with more detail provided for each step afterwards.

## Folder Structure
The folder `src` contains all of todo

# Work Flow
1. Run `newProject.py` to easily set up all the config files, and skip to step 4.
To do the initial setup manually, go to step 2.
2. Make a folder for the new project, and link the makefile from `heftCode/makefile` in this folder.
3. Enter the details of the system into `HEFT.config`.
4. Set up an initial guess for the fit parameters in `allFits.params`.
   If I'm trying to replicate the results of a paper, I will then get the fit parameters from that paper and put them into `allFits.params`, then skip to step TODO.
5. Get some scattering data to fit to, generally the WI08 solution from [SAID](https://gwdac.phys.gwu.edu/).
   This should be in a file called `dataInf.in` (see the example in this folder)
6. Setup `HEFTFitting.config` with the bounds on each category of fit parameter, and which parameters will be active in the fitting process.
7. If I want to run a bounded fit, make and run `fitBQ.x`.
   This program uses Powell's BOBYQA algorithm to minimise a $\chi^2$ between the HEFT phase shift and inelasticity, and the scattering data.
   If I don't want to use bounds, make and run `fit.x`, which uses minfun (Powell's NEWUOA algorithm).
8. Make sure `HEFTInfinite.config` is correct, with the correct energy range.
9. Set `iParamChoice` in `allFits.params` to the new parameter set found by the fit procedure, then make and run `inf.x` to calculate the scattering observables for this parameter set.
   I usually plot these using a file called `PhasevE.py`, which is copied from another project I've worked on.
10. To search for pole positions, you can make and run `poles.x`.
	This does a grid search from $ 1.0 - 0.01i $ to around $ 2.0 - 0.2i $ GeV, and should pick up the positions of any poles in the $ T $-matrix.
11. In a finite-volume, we're typically interested in how the eigenvalues and eigenvectors of the finite Hamiltonian vary with $m_{\pi}^{2}$.
	Finite-volume config options such as the $m_{\pi}^2$ range to be calculated are found in `HEFTFinite.config`.
12. To find the slope of the bare mass(es), I need lattice data to fit to.
	I store this in a file called `slopeFitMasses.data`.
13. To fit the bare mass slope(s), I make and run `fitBare.x`.
	The config file `HEFTFinite.config` contains the guess for the bare slope, and the order of the bare mass extrapolation.
14. To calculate the finite-volume energy spectrum, there are two options.
	The usual case is that the lattice size $L$ is varying with $m_{\pi}^{2}$, so I run `lqcdFin.x`.
	If I want to keep $L$ constant and just vary $m_{\pi}^{2}$, I run `mpiFin.x`.
15. To plot the spectrum, I have a plot script called `EvMpi.py`.
	Typically I'll have a few of these for different sets of lattice data, such as in the odd-parity nucleon case where I have `EvMpi_3fm.py`, `EvMpi_2fm.py`, and `EvMpi_D200.py`.
16. To plot the eigenvectors, I have a script called `EvecvMpi.py`.
17. **Extra:** There are also a few other things I might be interested in, such as comparing multiple parameter sets, such as in the $\Delta(1232)$ regularisation and multiple bare state papers.
  I can search for poles with a varying set of parameters with `multiFitPoles.x`.
  I calculate the finite-volume spectrum for a varying set of parameters with `multiFitFin.x`.

## 1. newProject.py
todo

## 3. HEFT.config
I had the general goal over my PhD of being able to do an entire new HEFT study having to change as little as possible in the actual code.
This goal was achieved to some degree, where a majority of the configuration specific to a system is stored in a file called `HEFT.config`.
An example of this file for the $S_{11}$ 2b3c (odd-parity nucleons) system is shown here.

**HEFT.config for $S_{11}$ 2b3c:**
```
 1  # NStar_2b3c config
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
20  B                    B          A
```
1.
    Lines 3 and 4 simply denote how many two-particle channels are present, and how many bare basis states are present in the system.
2.
    Lines 7 through 8 describe each two-particle channel in the system.
    The details of each hadron are stored in `heftCode/particles.in}, such as the mass, and the mass slope (how it varies in $m_{\pi}^{2}$).
    The partial wave label tells the program which partial wave the channel is interacting in, and is used for things such as choosing whether to start at $k=0$ or $k=1$ in a finite-volume, or having factors of $k^{l}$ for angular momentum $l$ in the potentials.
3.
    Line 11 tells the program which of those channel is on-shell in the scattering process (starting from 1).
    Generally this is going to be the $ \pi N $ channel, though in the $\Lambda^{*}(1405)$ case it would be a $\bar{K}N$ channel.
4.
    Lines 12 and 13 tell the program which file to use for the particle masses.
    By default, it uses `heftCode/particles.in`.
    When I'm doing a finite-volume study where I want to have slightly larger masses at the physical point, I set `useCustomMasses` to `T`, and give the file name of the new particles file to use.
5.
    Lines 16 and 17 are just the labels I give to each bare state, which is helpful for printing/writing data to files.
6.
    Line 20 is where I set which potential is used for each interaction.
    There are three customisation options here, which describe the functions used for $G_{\alpha}^{B_{0}}(k)$, $V_{\alpha\beta}(k,k')$, and the regulator $u(k,\Lambda)$.
    The labels `A`, `B`, etc. correspond to functions in `heftCode/heft.f90`, with names such as `u_k_A`, `g_k_B`, and `f_k_A`.
    Note that `f_k` represents the separable part of $V_{\alpha\beta}(k,k') = v_{\alpha\beta}\, f_{\alpha}(k)\, f_{\beta}(k')$.


## 4. allFits.params
With the general properties of the system of interest set out in \ttt{HEFT.config}, we're able to generate the file which will contain all of the parameter sets for the system.
These parameters are the bare mass(es) $m_{B_{0}}^{(0)}$, the couplings between bare basis states and two-particle basis states (2-1 couplings) $g_{\alpha}^{B_{0}}$, the regulator parameters for 2-1 interactions $\Lambda_{\alpha}^{B_{0}}$, the couplings between two-particle basis states (2-2 couplings) $v_{\alpha\beta}$, and the regulator parameters for 2-2 interactions $\Lambda_{v,\alpha}$.

**allFits.params for $P_{33}$ 1b1c:**
```
1      nParam        5
2      iChoice       1
3      paramEnds  1 2 3 4 5
4      n    m_bare   g_Delta   Lambda    v_piN     Lambda_v     chi2    Notes
5      1    1.4000   0.10000   0.8000   -0.10000   0.8000        0.00   Default
6      2    1.3588   0.17615   0.8000   -0.02863   0.8000      236.81   AAA
7      3    1.4073   0.00000   0.8000   -0.00290   8.0000    24373.49   AAA
8      4    1.3850   0.14057   0.8000   -0.03066   0.8000      244.65   AAB
9      0    0.0000   0.00000   0.0000    0.00000   0.0000        0.00
```

-
  `nParam` refers to the total number of parameters.
  For $n_{b}$ bare basis states, and $n_{c}$ scattering channels, considering all five categories of fit parameters gives

   $n_{\text{param}} = n_{b} + n_{b} n_{c} + n_{b} n_{c} + \frac{1}{2} n_{c}(n_{c}+1) + n_{c} = n_{b} + 2 n_{b} n_{c} + \frac{1}{2} n_{c}(n_{c}+1) + n_{c}$

-
  `iChoice` chooses which parameter set to use (given by the `n` column)
-
  `paramEnds` denotes which column each parameter category ends on.
  This is trivial for the 1b1c case, but for the 2b3c case (see `allFits_N2b3c.params`) this would be `2  8  14  20  23`.
  This refers to $m_{B_{0}}^{(0)}$ ending on column 2, $g_{\alpha}^{B_{0}}$ ending on column 8, $\Lambda_{\alpha}^{B_{0}}$ ending on column 14 etc.
-
  For lines 5 through 8, each of the numbered rows are a parameter set, with parameters in the order $m_{B_{0}}^{(0)}$, $g_{\alpha}^{B_{0}}$, $\Lambda_{\alpha}^{B_{0}}$, $v_{\alpha\beta}$, $\Lambda_{v, \alpha}$.
  Following this is the $\chi^{2}$ for the fit, then any notes about the fit I've added.
  Usually this is just which set of potentials and regulators that fit corresponds with.
-
  The last row is always just a row of zeros, which makes it each to detect the end of file when reading/writing.
  As an example the third `A` corresponds with a dipole regulator, while a `B` in that spot would be a Guassian regulator.
  Note that these are just a note for myself, and do not set what the program uses.
  The actual set of potentials/regulators are set in `HEFT.config}.
