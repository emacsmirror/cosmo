;;; cosmo.el --- Cosmological Calculator    -*- lexical-binding: t; -*-

;; Copyright (C) 2017 Francesco Montanari
;;
;; Author: Francesco Montanari <fmnt@fmnt.info>
;; Maintainer: Francesco Montanari <fmnt@fmnt.info>
;; Created: 22 April 2017
;; Version: 0.1
;; Keywords: tools
;; Homepage: https://gitlab.com/montanari/cosmo-el
;;
;; This file is not part of GNU Emacs.
;;
;; This program is free software: you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or
;; (at your option) any later version.
;;
;; This program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.
;;
;; You should have received a copy of the GNU General Public License
;; along with this program.  If not, see <http://www.gnu.org/licenses/>.

;;; Commentary:

;; This package provides a cosmological calculator for Lambda-CDM
;; models.  Such a framework describes a homogeneous and isotropic
;; universe containing a cosmological constant (Lambda) and a Cold
;; Dark Matter (CDM) component, besides ordinary species.  The
;; model is characterized by the following parameters:
;;
;; - H_0 :: Hubble parameter (expansion rate) today.
;; - Omega_m0 :: Matter density parameter today.
;; - Omega_Lambda :: Cosmological constant density parameter.
;; - Omega_r0 :: Relativistic species (e.g., photons plus
;;               neutrinos) density parameter today.
;; - Omega_k0 :: Curvature density parameter today.  This
;;               parameter is derived from the others above
;;               according to Friedmann's equation
;;               Omega_m0 + Omega_Lambda + Omega_r0 + Omega_k0 = 1.
;;
;; All cosmological quantities are computed at a given redshift
;; value:
;;
;; - redshift :: Gravitational redshift of photons frequency due to the
;;               expansion of the Universe.
;;
;; Definitions follow Hogg (1999)
;; <https://arxiv.org/abs/astro-ph/9905116>.
;;
;; Names with "--" are for functions and variables that are meant to
;; be for internal use only.

;;; Bugs:

;; - None known.

;;; Todo:

;; - Suggest default parameters when reading them with the related
;;   command; set the to default values if none is entered.
;;
;; - (Consider the following only if performance becomes critical.) At
;;   a fixed redshift, only cosmo-get-los-comoving-distance perform
;;   the integral. Other distances are defined in terms of this
;;   one. However, each time that other distances are required this
;;   distance is called once again. This is not necessary, it can be
;;   called just once for a given redshift (or use memoization).

;;; Code:

;;; Global variables.

(defgroup cosmo nil
  "Cosmological calculator."
  :group 'convenience)

(defcustom cosmo-int-prec 1e-3
  "Fractional accuracy for numerical integrals."
  :group 'cosmo)

(defcustom cosmo-int-maxsteps 20
  "Maximum number of steps in the numerical integral algorithm."
  :group 'cosmo)

(defvar cosmo-pedia-text
  ;; Text appearing in the *Cosmopedia* buffer.
  "(`q` to quite)

Units system: hbar = c = k_Boltzmann = 1.

Distances relations
-------------------
- Comoving distance (transverse): D_M
- Angular diameter distance:      D_A = D_M / (1+z)
- Luminosity distance:            D_L = (1+z) D_M = (1+z)^2 D_A

Conversion factors, units
-------------------------
- 1GeV = 1.6022e-3 erg
       = 1.1605e13 K
       = 1.7827e-24 g
       = 5.0684e13 1/cm
       = 1.5192 1/s
- 1 pc = 3.2612 light years
       = 3.0856e18 cm
- 1 Mpc = 1e6pc ~ 3e24 cm ~ 1e14 s
- 1 AU = 1.4960e13 cm
- 1 Jy = 1e-23 erg/cm^2/s/Hz
       = 2.4730e-48 GeV^3
- 1 yr ~ pi * 1e17 s

Important constants
-------------------
- Hubble constant: H0 = 100h km/s/Mpc
                      = 2.1332e-42 h GeV
- Hubble time, distance: 1/H0 = 3.0856e17/h s
                               = 9.7776e9/h yr
                               = 2997.9/h Mpc
                               = 9.2503e27/h cm")

;; Hash table containing all independent cosmological parameters.
(defvar cosmo--params
  (let ((table (make-hash-table :test #'equal)))
    (puthash "H0 [Km/s/Mpc]" 70.0 table) ; Hubble today km/s/Mpc.
    (puthash "omatter" 0.3 table)        ; Matter density today.
    (puthash "olambda" 0.7 table)        ; Curvature density today.
    (puthash "orel" 0.0 table)           ; Relativistic density today.
    table)
  "Table containing Lambda-CDM cosmological parameters.")

;;; Handle input.

(defun cosmo--string-number-p (string)
  "Test whether STRING represents a number."
  (if (string-match "\\`[-+]?[0-9]+\\.?[0-9]*\\'" string)
      t
    nil))

(defun cosmo--read-param (name)
  "Read parameter NAME from minibuffer and convert it to a number."
  (let ((value (read-from-minibuffer (format "Enter %s: " name))))
    (if (cosmo--string-number-p value)
        (string-to-number value)
      (error "Error: parameter must be a number"))))

(defun cosmo--put-param (name)
  "Read parameter NAME from minibuffer and add it to the parameter table."
  (puthash name (cosmo--read-param name) cosmo--params))

(defun cosmo--check-param (name value)
  "Check the validity of NAME (a cosmological parameter) VALUE."
  (cond ((or (string= name "omatter")
             (string= name "olambda")
             (string= name "orel"))
         (unless (>= value 0.0)
           (error "Error: density parameter must be positive")))))

(defun cosmo-set-params ()
  "Change the values of cosmological parameters."
  (interactive)
  (maphash (lambda (key _value)
             (cosmo--put-param key)
             (cosmo--check-param key (gethash key cosmo--params)))
           cosmo--params))

;;; Numerical utilities.

(defun cosmo-sinh (x)
  "Hyperbolic sine of real arguments X."
  (* 0.5 (- (exp x) (exp (- x)))))

(defun cosmo-asinh (x)
  "Inverse hyperbolic sine of real arguments X."
  (log (+ x (sqrt (1+ (* x x))))))

(defun cosmo--trapzd (func lo hi sum niter)
  "Extended trapezoidal rule to integrate FUNC from LO to HI.
Compute the refinement to an old SUM result, previous to the
current NITER stage of refinement.  The function must be called
for NITER = 1, 2, 3, ... in order.  Boundaries LO and HI are
always converted to floats."
  (let ((lo (float lo))                 ; Always assume float boundaries.
        (hi (float hi)))
    ;; NITER=1 corresponds simply to integrate(func(x), x, lo, hi)
    (if (= niter 1)
        (* 0.5 (- hi lo) (+ (funcall func lo) (funcall func hi)))
      ;; Subsequent calls NITER = 2, 3, ... (in order) improve the
      ;; accuracy by adding 2^(NITER-2) interior points.
      (let* ((iter (expt 2 (- niter 2)))
             (step (/ (- hi lo) iter))
             (xval (+ lo (* 0.5 step)))
             (sumf 0.0)
             (j 0))
        (while (< j iter)
          (setq sumf (+ sumf (funcall func xval)))
          (setq xval (+ xval step))
          (setq j (1+ j)))
        (* 0.5 (+ sum (* (- hi lo) (/ sumf iter))))))))

(defun cosmo-qsimp (func lo hi &optional eps jmax)
  "Simpson's integration.
Argument FUNC integrand.
Argument LO lower bound.
Argument HI upper bound.
Optional argument EPS fractional accuracy.
Optional argument JMAX maximum number of steps."
  (let ((eps (or eps 1e-6))             ; Set defaults.
        (jmax (or jmax 20))
        (nsum 0.0)                      ; N-th iteration sum.
        (nsumt 0.0)
        (osum 0.0)                      ; Old sum.
        (osumt 0.0)
        (converged nil)
        (j 1))
    (while (and (<= j jmax) (not converged))
      (setq nsumt (cosmo--trapzd func lo hi nsum j))
      (setq nsum (/ (- (* 4.0 nsumt) osumt) 3.0))
      (if (and (> j 5)                  ; Avoid spurious early convergence.
               (or (< (abs (- nsum osum)) (* eps (abs osum)))
                   (and (= nsum 0.0) (= osum 0.0))))
          (setq converged t))
      (setq osum nsum)
      (setq osumt nsumt)
      (setq j (1+ j)))
    (if converged
        nsum
      (error "Error: the integral did not converge, try to
      increase the number of steps"))))

;;; Compute cosmological functions.

(defun cosmo-get-ocurvature ()
  "Get curvature density parameter today from Friedmann equations."
  (let ((omatter (gethash "omatter" cosmo--params))
        (olambda (gethash "olambda" cosmo--params))
        (orel (gethash "orel" cosmo--params)))
    (- 1.0 omatter olambda orel)))

(defun cosmo-efunc (redshift)
  "E(z) function at a given REDSHIFT."
  (let ((omatter (gethash "omatter" cosmo--params))
        (olambda (gethash "olambda" cosmo--params))
        (orel (gethash "orel" cosmo--params))
        (ocurvature (cosmo-get-ocurvature))
        (zp1 (+ 1 redshift)))
    (sqrt (+ (* orel (expt zp1 4.0))
             (* omatter (expt zp1 3.0))
             (* ocurvature (expt zp1 2.0))
             olambda))))

(defun cosmo-inv-efunc (redshift)
  "Inverse E(z) function at a given REDSHIFT."
  (/ 1.0 (cosmo-efunc redshift)))

(defun cosmo-get-hubble (redshift)
  "Hubble parameter [Km/s/Mpc] for Lambda-CDM at a given REDSHIFT."
  (let ((H0 (gethash "H0 [Km/s/Mpc]" cosmo--params)))
    (* H0 (cosmo-efunc redshift))))

(defun cosmo-hubble ()
  "Display Hubble parameter in mini-buffer."
  (interactive)
  (let ((z (cosmo--read-param "redshift")))
    (message (format "%s km/s/Mpc" (cosmo-get-hubble z)))))

(defun cosmo-get-hubble-distance ()
  "Hubble distance c/H0 [Mpc] for Lambda-CDM."
  (let ((H0 (gethash "H0 [Km/s/Mpc]" cosmo--params)))
    (/ 3.0e5 H0)))

(defun cosmo-hubble-distance ()
  "Display Hubble distance c/H0 [Mpc] in mini-buffer."
  (interactive)
  (message (format "%s Mpc" (cosmo-get-hubble-distance))))

(defun cosmo-get-hubble-time ()
  "Hubble time 1/H0 [Gyr] for Lambda-CDM."
  (let ((H0 (gethash "H0 [Km/s/Mpc]" cosmo--params)))
    (/ 9.78e2 H0)))

(defun cosmo-hubble-time ()
  "Display Hubble distance 1/H0 [yr] in mini-buffer."
  (interactive)
  (message (format "%s Gyr" (cosmo-get-hubble-time))))

(defun cosmo-get-los-comoving-distance (redshift)
  "Line-of-sight comoving distance [Mpc] for Lambda-CDM at a given REDSHIFT."
  (let ((DH (cosmo-get-hubble-distance))
        (int (cosmo-qsimp #'cosmo-inv-efunc 0.0 redshift
                          cosmo-int-prec cosmo-int-maxsteps)))
    (* DH int)))

(defun cosmo-los-comoving-distance ()
  "Display line-of-sight comoving distance in mini-buffer."
  (interactive)
  (let ((z (cosmo--read-param "redshift")))
    (message (format "%s Mpc" (cosmo-get-los-comoving-distance z)))))

(defun cosmo-get-transverse-comoving-distance (redshift)
  "Line-of-sight comoving distance [Mpc] for Lambda-CDM at a given REDSHIFT."
  (let* ((DH (cosmo-get-hubble-distance))
         (DC (cosmo-get-los-comoving-distance redshift))
         (ocurvature (cosmo-get-ocurvature))
         (sqrt-ok (sqrt (abs (cosmo-get-ocurvature))))
         (DH-over-sqrtok (/ DH sqrt-ok)))
    (cond ((> ocurvature 0)
           (* DH-over-sqrtok (cosmo-sinh (/ DC DH-over-sqrtok))))
          ((= ocurvature 0)
           DC)
          ((< ocurvature 0)
           (* DH-over-sqrtok (sin (/ DC DH-over-sqrtok)))))))

(defun cosmo-transverse-comoving-distance ()
  "Display transverse comoving distance in mini-buffer."
  (interactive)
  (let ((z (cosmo--read-param "redshift")))
    (message (format "%s Mpc"
                     (cosmo-get-transverse-comoving-distance z)))))

(defun cosmo-get-angular-diameter-distance (redshift)
  "Angular diameter distance [Mpc] for Lambda-CDM at a given REDSHIFT."
  (let* ((DM (cosmo-get-transverse-comoving-distance redshift)))
    (/ DM (1+ redshift))))

(defun cosmo-angular-diameter-distance ()
  "Display angular diameter distance in mini-buffer."
  (interactive)
  (let ((z (cosmo--read-param "redshift")))
    (message (format "%s Mpc"
                     (cosmo-get-angular-diameter-distance z)))))

(defun cosmo-get-luminosity-distance (redshift)
  "Luminosity distance [Mpc] for Lambda-CDM at a given REDSHIFT."
  (let* ((DM (cosmo-get-transverse-comoving-distance redshift)))
    (* DM (1+ redshift))))

(defun cosmo-luminosity-distance ()
  "Display luminosity distance in mini-buffer."
  (interactive)
  (let ((z (cosmo--read-param "redshift")))
    (message (format "%s Mpc"
                     (cosmo-get-luminosity-distance z)))))

(defun cosmo-get-parallax-distance (redshift)
  "Parallax distance [Mpc] for Lambda-CDM at a given REDSHIFT."
  (let* ((DH (cosmo-get-hubble-distance))
         (DM (cosmo-get-transverse-comoving-distance redshift))
         (dM (/ DM DH))
         (ocurvature (cosmo-get-ocurvature)))
    (/ DM (+ dM (sqrt (1+ (* ocurvature dM dM)))))))

(defun cosmo-parallax-distance ()
  "Display parallax distance in mini-buffer."
  (interactive)
  (let ((z (cosmo--read-param "redshift")))
    (message (format "%s Mpc"
                     (cosmo-get-parallax-distance z)))))

(defun cosmo--get-comoving-volume-nonflat (DM DH ocurvature)
  "Return the comoving volume for non-vanishing curvature.
Argument DM comoving distance (transverse).
Argument DH Hubble distance.
Argument OCURVATURE curvature density parameter."
  (let* ((DM-over-DH (/ DM DH))
         (sqrt-ok (sqrt (abs (cosmo-get-ocurvature))))
         (pref (* 2.0 float-pi (/ (expt DH 3.0) ocurvature)))
         (func (cond ((> ocurvature 0.0)
                      #'cosmo-asinh)
                     ((< ocurvature 0.0)
                      #'asin)
                     (t (error "Error: wrong curvature parameter value")))))
    (setq term1 (* DM-over-DH
                   (sqrt (1+ (* ocurvature (expt DM-over-DH 2.0))))))
    (setq term2 (/ (funcall func (* sqrt-ok DM-over-DH)) sqrt-ok))
    (* pref (- term1 term2))))

(defun cosmo-get-comoving-volume (redshift)
  "Comoving volume [Mpc^3] for Lambda-CDM at a given REDSHIFT."
  (let* ((DH (cosmo-get-hubble-distance))
         (DM (cosmo-get-transverse-comoving-distance redshift))
         (ocurvature (cosmo-get-ocurvature)))
    (if (= ocurvature 0.0)
        (* (/ 4.0 3.0) float-pi (expt DM 3.0))
        (cosmo--get-comoving-volume-nonflat DM DH ocurvature))))

(defun cosmo-comoving-volume ()
  "Display comoving volume in mini-buffer."
  (interactive)
  (let ((z (cosmo--read-param "redshift")))
    (message (format "%s Mpc^3"
                     (cosmo-get-comoving-volume z)))))

;;; Handle output.

(defun cosmo--write-calc-header ()
  "Write header for the cosmological calculator summary buffer."
  (let ((head "Cosmology calculator.\n\n")
        (help "(`q` to quite)\n\n"))
    (insert (propertize help 'font-lock-face 'italic))
    (insert head)))

(defun cosmo--write-calc (redshift H0 omatter olambda orel hubble
                                   los-dist transverse-dist
                                   luminosity-dist angular-dist)
  "Format and insert cosmological table in buffer.
Argument REDSHIFT redshift.
Argument H0 Hubble parameter today.
Argument OMATTER matter density parameter.
Argument OLAMBDA cosmological constant density parameter.
Argument OREL density parameter.
Argument HUBBLE Hubble parameter at given redshift.
Argument LOS-DIST line-of-sight comoving distance at given redshift.
Argument TRANSVERSE-DIST transverse comoving distance at given redshift.
Argument LUMINOSITY-DIST luminosity distance at given redshift.
Argument ANGULAR-DIST angular diameter distance at given redshift."
  ;; Input parameters.
  (cosmo--write-calc-header)
  (insert "Input Parameters\n"
          "----------------\n"
          (format "- Redshift:                                 %s\n"
                  redshift)
          (format "- Hubble constant, now [km/s/Mpc]:          %s\n"
                  H0)
          (format "- Matter fractional density, now:           %s\n"
                  omatter)
          (format "- Cosmological constant fractional density: %s\n"
                  olambda)
          (format "- Relativistic fractional density, now:     %s\n"
                  orel)
          "\n")
  ;; Derived parameters.
  (insert "Derived parameters\n"
          "------------------\n"
          (format "- Curvature fractional density: %s\n"
                  (cosmo-get-ocurvature))
          (format "- Hubble distance [Mpc]:        %s\n"
                  (cosmo-get-hubble-distance))
          (format "- Hubble time [Gyr]:            %s\n"
                  (cosmo-get-hubble-time))
          "\n")
  ;; Cosmological functions.
  (insert "Cosmography at required redshift\n"
          "--------------------------------\n"
          (format "- Hubble parameter [km/s/Mpc]:             %s\n"
                  hubble)
          (format "- Comoving distance (line-of-sight) [Mpc]: %s\n"
                  los-dist)
          (format "- Comoving distance (transverse) [Mpc]:    %s\n"
                  transverse-dist)
          (format "- Angular diameter distance [Mpc]:         %s\n"
                  angular-dist)
          (format "- Luminosity distance [Mpc]:               %s\n"
                  luminosity-dist))
  nil)

(defun cosmo-calculator ()
  "Compute cosmology and display summary table in a new buffer."
  (interactive)
  (let* ((cosmo-buffer "*Cosmo*")
         (redshift (cosmo--read-param "redshift"))
         (omatter (gethash "omatter" cosmo--params))
         (olambda (gethash "olambda" cosmo--params))
         (orel (gethash "orel" cosmo--params))
         (H0 (gethash "H0 [Km/s/Mpc]" cosmo--params))
         (hubble (cosmo-get-hubble redshift))
         (los-dist (cosmo-get-los-comoving-distance redshift))
         (transverse-dist (cosmo-get-transverse-comoving-distance redshift))
         (luminosity-dist (cosmo-get-luminosity-distance redshift))
         (angular-dist (cosmo-get-angular-diameter-distance redshift)))
    (with-output-to-temp-buffer cosmo-buffer
      (pop-to-buffer cosmo-buffer)
      (cosmo--write-calc redshift H0 omatter olambda orel hubble
                         los-dist transverse-dist luminosity-dist
                         angular-dist))))

(defun cosmo-pedia ()
  "Display a reference to basic cosmological definitions."
  (interactive)
  (let* ((cosmo-buffer "*Cosmopedia*"))
    (with-output-to-temp-buffer cosmo-buffer
      (pop-to-buffer cosmo-buffer)
      (insert cosmo-pedia-text))))

(provide 'cosmo)

;;; cosmo.el ends here
