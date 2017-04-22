                         _____________________

                                COSMO.EL

                          Francesco Montanari
                         _____________________


Table of Contents
_________________

1 Description
2 Download
3 Installation
4 Usage
5 Support
6 Contributing


cosmo.el --- Cosmological Calculator

Copyright (C) 2017 Francesco Montanari

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program.  If not, see [http://www.gnu.org/licenses/].


1 Description
=============

  This package provides a cosmological calculator integrated into the
  Emacs text editor.

  Cosmological quantities are computed for a Lambda-CDM model. Such a
  framework describes a homogeneous and isotropic universe containing a
  cosmological constant (Lambda) and a Cold Dark Matter (CDM) component,
  besides ordinary species. The model is characterized by the following
  parameters:
  H_0: Hubble parameter (expansion rate) today
  Omega_m0: Fractional matter density parameter today.
  Omega_Lambda: Cosmological constant density parameter. At the moment we
                consider only flat Lambda-CDM models, for which
                `Omega_Lambda = 1 - Omega_m0'.

  All cosmological quantities are computed at a given redshift value:
  redshift: Gravitational redshift of photons frequency due to the
            expansion of the Universe.


2 Download
==========

  The source code can be downloaded from the [git repository page].

  Alternatively, the git repository can be cloned:
  ,----
  | git clone git@gitlab.com:montanari/cosmo-el.git
  `----


  [git repository page] https://gitlab.com/montanari/cosmo-el


3 Installation
==============

  Copy the `cosmo-el' file into your Emacs configuration
  directory. E.g., on GNU systems:

  ,----
  | mkdir ~/.emacs.d/cosmo-el
  | cp cosmo.el ~/.emacs.d/cosmo-el
  `----

  Open your Emacs configuration file (e.g., `~\/.emacs' on GNU systems),
  and add the following lines:

  ,----
  | ;; Cosmological calculator
  | (add-to-list 'load-path "~/.emacs.d/cosmo-el/")
  | (require 'cosmo)
  `----


4 Usage
=======

  Open Emacs and type `M-x cosmo-command', where the `command' name
  varies depending on the desired computation. For example:

  cosmo-set-params: Set cosmological parameters.

  cosmo-calculator: Compute cosmology and display summary table in a new
                    buffer. The following illustrates a possible output:
                    ,----
                    | Cosmology calculator.
                    |
                    | Input Parameters
                    | ----------------
                    | - Redshift:                       	0.1
                    | - Hubble constant, now [km/s/Mpc]:	70
                    | - Matter fractional density, now: 	0.3
                    | <...>
                    |
                    | Cosmography at required redshift
                    | -------------------------------
                    | - Hubble parameter [km/s/Mpc]:	73.39325582095401
                    | <...>
                    `----

  A complete list of interactive commands can be obtained by typing `M-x
  cosmo- TAB'. Documentation is available through `C-h f cosmo-command',
  where `command' should be adapted to the particular command.


5 Support
=========

  Bugs and issues are tracked through the [git repository page]. Please
  see [this page] about how to report bugs effectively.


  [git repository page] https://gitlab.com/montanari/cosmo-el

  [this page] http://www.chiark.greenend.org.uk/~sgtatham/bugs.html


6 Contributing
==============

  /The project is still at an early stage/. Recommendations (and
  contributions) aimed to improve the source code are highly
  appreciated. New feature suggestions are also welcome, but at this
  point priority will be given to reach an idiomatic and extensible
  code.

  Contributions can be submitted as patches. See [this page] for an
  example of good patches contributions.

  More substantial contributions should proceed through git
  [Integration-Manager Workflow]. See [this page] for an example of a
  complete working session.


  [this page] http://orgmode.org/worg/org-contribute.html#patches

  [Integration-Manager Workflow]
  https://git-scm.com/book/en/v2/Distributed-Git-Distributed-Workflows

  [this page]
  https://www.gnu.org/software/gnuastro/manual/html_node/Contributing-to-Gnuastro.html
