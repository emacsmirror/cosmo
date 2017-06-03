"""TODO: instead of producing Elisp code, create a comparison
table. Ideally, load Elisp code from within this script to create the
respective table, and assert if the two tables are almost equal. At
worst this should be doable by calling a Elisp script from python, as
a shell command.

"""
from astropy.cosmology import LambdaCDM

LISPFMT = """{NAME}
------
(setq redshifts '({REDS}))
(setq {NAME} '({VAL}))

"""

def _list_to_str(list_):
    """Return a string containing the list arguments without enclosing
    parenthesis.
    """
    list_str = str(list_)
    for par in ('(', ')', '[', ']', '{', '}'):
        list_str = list_str.replace(par,'')
    return list_str

def efunc(redshifts, cosmo):
    reds = _list_to_str(redshifts)
    val = _list_to_str(cosmo.efunc(redshifts))
    print(LISPFMT.format(REDS=reds, NAME='efunc',
                         VAL=val))

def inv_efunc(reds, cosmo):
    reds = _list_to_str(redshifts)
    val = _list_to_str(cosmo.inv_efunc(redshifts))
    print(LISPFMT.format(REDS=reds, NAME='inv-efunc',
                         VAL=val))

def comoving_distance(reds, cosmo):
    reds = _list_to_str(redshifts)
    val = _list_to_str(cosmo.comoving_distance(redshifts).value)
    print(LISPFMT.format(REDS=reds, NAME='comoving-los-distance',
                         VAL=val))


if __name__ == '__main__':
    # Remember to setq photons+neutrinos (Omega_rel*h^2 = 4.18e-05) in
    # cosmo.el tests.
    cosmo = LambdaCDM(H0=70., Om0=0.3, Ode0=0.7)
    redshifts = (0., 0.1, 10., 1000.)
    efunc(redshifts, cosmo)
    inv_efunc(redshifts, cosmo)
    comoving_distance(redshifts, cosmo)
