"""test_cosmo -- Comparison values to test cosmo.el

Copyright (C) 2017 Francesco Montanari

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from astropy.cosmology import LambdaCDM
import numpy as np


def sexp_print(reds, vals, name=' '):
    """Print redshift list and respective values list as sexp."""
    vals_str = '(' + ' '.join(map(str, vals)) + ')'
    print('{}:\n  {}'.format(name, vals_str))

    reds_str = '(' + ' '.join(map(str, reds)) + ')'
    print('redshifts:\n  {}'.format(reds_str))

    print('')


def print_curves(reds, cosmo):
    """Print cosmological functions given redshifts and a cosmology."""
    sexp_print(reds, cosmo.efunc(reds), name='efunc')
    sexp_print(reds, cosmo.inv_efunc(reds), name='inv-efunc')
    sexp_print(reds, cosmo.H(reds).value,
               name='hubble [km/s/Mpc]')
    sexp_print(['-'], [cosmo.hubble_distance.value],
               name='hubble-distance [Mpc]')
    sexp_print(['-'], [cosmo.hubble_time.value],
               name='hubble-time [Gyr]')
    sexp_print(reds, cosmo.comoving_distance(reds).value,
               name='los-comoving-distance [Mpc]')
    sexp_print(reds, cosmo.comoving_transverse_distance(reds).value,
               name='transverse-comoving-distance [Mpc]')
    sexp_print(reds, cosmo.luminosity_distance(reds).value,
               name='luminosity-distance [Mpc]')
    sexp_print(reds, cosmo.angular_diameter_distance(reds).value,
               name='angular-diameter-distance [Mpc]')


def print_open_cosmo():
    """Print cosmological functions for an open cosmology."""
    H0 = 70.
    Om0 = 0.31
    Ode0 = 0.7
    cosmo = LambdaCDM(H0=H0, Om0=Om0, Ode0=Ode0)

    print('* Open cosmology:')
    print('H0={} [km/s/Mpc], Om0={}, Ode0={}, Orel0={}, Ok0={}\n'
          .format(H0, Om0, Ode0, cosmo.Ogamma(0.)+cosmo.Onu(0.),
                  cosmo.Ok(0.)))

    redshifts = (0., 0.1, 10., 1000.)
    print_curves(redshifts, cosmo)


if __name__ == '__main__':
    print_open_cosmo()
