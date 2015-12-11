Agile-IDSA
==========

A simple supernova model
------------------------

The Agile-IDSA code provides tools to run a rudimentary and approximate
model of a core-collapse supernova with neutrino transport in
spherical symmetry through the phases of stellar collapse, bounce,
and early postbounce evolution. The code is not fool-proof, requires
undocumented customisation to the local computational infrastructure,
and is meant for people who work, or get started to work, in the
field of supernova theory. Of course we are pleased if Agile-IDSA
turns out to be useful also for investigators of supernova input
physics, or if parts of the code can be extended to broader
astrophysical application.

The separate hydrodynamics part Agile (Adaptive Grid with Implicit
Leap Extrapolation) has been used in general relativistic spherically
symmetric models of stellar core collapse and postbounce evolution
with Boltzmann neutrino transport (e.g. Liebendörfer et al. 2004,
ApJS 150, 263), while the neutrino transport part based on the
isotropic diffusion source approximation (IDSA) has been used in
multidimensional models of stellar core collapse and postbounce
evolution (e.g. Whitehouse & Liebendörfer 2008, PoS(NIC X)243; Suwa
et al. 2011, ApJ 738, 165). Here we put the two pieces together to
a standalone and physically simple spherically symmetric supernova
model and possible toolkit for further projects. The code implements
physics outlined in the following publications:

Hydrodynamics: Liebendörfer, Rosswog & Thielemann 2002, ApJS 141, 229.
Neutrino parameterisation: Liebendörfer 2005, ApJ 633, 1042.
Neutrino transport: Liebendörfer, Whitehouse & Fischer 2009, ApJ 698, 1174.
Weak interactions: Mezzacappa & Bruenn 1993, ApJ 405, 637.
Equation of state: Lattimer & Swesty 1991, Nuc. Phys. A 535, 331. K=220 MeV.

Distribution
------------

Current version: agile-idsa_110930.tar.gz. In our git repository, this
corresponds to commit c10bdcf8ff3404f36be399eab00e2d138c49f7d0.

Agile-IDSA is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Acknowledgement
---------------

This project was mainly funded by the Swiss National Science
Foundation under grant No 20-47252.96, 20-53798.98, PP002-106627/1
and PP00P2_124879/1. We equally acknowledge the stimulating and
supportive research-environment at the Canadian Institute for
Theoretical Astrophysics in Toronto and within the ESF Research
Networking Programme 'CompStar'.
