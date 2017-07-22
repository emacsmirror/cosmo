;;; cosmo.el --- Cosmological Calculator    -*- lexical-binding: t; -*-

;; Copyright (C) 2017 Francesco Montanari
;;
;; Author: Francesco Montanari <fmnt@fmnt.info>
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

;; This package provides a cosmological calculator Lambda-CDM
;; models. Such a framework describes a homogeneous and isotropic
;; universe containing a cosmological constant (Lambda) and a Cold
;; Dark Matter (CDM) component, besides ordinary species. The
;; model is characterized by the following parameters:
;;
;; - H_0 :: Hubble parameter (expansion rate) today
;; - Omega_m0 :: Matter density parameter today.
;; - Omega_Lambda :: Cosmological constant density parameter.
;; - Omega_r0 :: Relativistic species (e.g., photons plus
;;               neutrinos) density parameter today.
;; - Omega_k0 :: Curvature density parameter today. This
;;               parameter is derived from the others above
;;               according to Friedmann's equation
;;               =Omega_m0 + Omega_Lambda + Omega_r0 + Omega_k0 = 1=.
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

;; - Add all distances from Hogg 1999.
;;
;; - Simpson's rule performs well for standard cosmologies and for z <
;;   1000. If non-standard cosmologies or very large redshifts are
;;   required, more steps may be required. They should be set
;;   automatically, the user should only pass a max number of
;;   steps.
;;
;; - Consider using Calc as a library for quadrature, special
;;   functions (sinh) and maybe to plot.
;;
;; - Suggest default parameters when reading them with the related
;;   command; set the to default values if none is entered.


;;; Code:

;;; Define cosmological parameters.

;; Hash table containing all independent cosmological parameters.
(defvar cosmo--params
  (let ((table (make-hash-table :test #'equal)))
    (puthash "H0 [Km/s/Mpc]" 70.0 table) ; Hubble today km/s/Mpc.
    (puthash "omatter" 0.3 table)        ; Matter density today.
    (puthash "olambda" 0.7 table)        ; Curvature density today.
    (puthash "orel" 0.0 table)           ; Relativistic density today.
    table)
  "Table containing Lambda-CDM cosmological parameters.")

;; Derived cosmological parameter.

(defun cosmo-get-ocurvature ()
  "Get curvature density parameter today from Friedmann equations."
  (let ((omatter (gethash "omatter" cosmo--params))
        (olambda (gethash "olambda" cosmo--params))
        (orel (gethash "orel" cosmo--params)))
    (- 1.0 omatter olambda orel)))

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

;;; Numerical utilities. (May be replaced by calling calc-eval.)

(defun cosmo-sinh (x)
  "Hyperbolic sine of real arguments X."
  (* 0.5 (- (exp x) (exp (- x)))))

(defun cosmo-simps (func a b &optional n)
  "Simpson rule.

Integrate a function FUNC of one argument from A to B in N steps.
The values A and B will be considered as float.

Example:
\(cosmo-simps #'\(lambda \(x\) x\) 0.0 1.0\)"
  (let ((a (float a))                   ; Extremes must be floats.
        (b (float b))
        (n (or n 50)))
    (loop with h = (/ (- b a) n)
          with sum1 = (funcall func (+ a (/ h 2)))
          with sum2 = 0
          for i from 1 below n
          do (incf sum1 (funcall func (+ a (* h i) (/ h 2))))
          do (incf sum2 (funcall func (+ a (* h i))))
          finally (return (* (/ h 6)
                             (+ (funcall func a)
                                (funcall func b)
                                (* 4 sum1)
                                (* 2 sum2)))))))


;;; Compute cosmological functions.

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

(defun cosmo-get-los-comoving-distance (redshift)
  "Line-of-sight comoving distance [Mpc] for Lambda-CDM at a given REDSHIFT."
  (let ((DH (cosmo-get-hubble-distance))
        (int (cosmo-simps #'cosmo-inv-efunc 0.0 redshift 1000)))
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

(defun cosmo-get-hubble-time ()
  "Hubble time 1/H0 [Gyr] for Lambda-CDM."
  (let ((H0 (gethash "H0 [Km/s/Mpc]" cosmo--params)))
    (/ 9.78e2 H0)))

(defun cosmo-hubble-time ()
  "Display Hubble distance 1/H0 [yr] in mini-buffer."
  (interactive)
  (message (format "%s Gyr" (cosmo-get-hubble-time))))

;;; Handle output.

(defun cosmo--write-calc-header ()
  "Write header for the cosmological calculator summary buffer."
  (let ((head "Cosmology calculator.\n\n")
        (help "(`q` to quite)\n\n"))
    (insert (propertize help 'font-lock-face 'italic))
    (insert head)))

(defun cosmo--write-calc (redshift H0 omatter olambda orel hubble los-dist)
  "Format and insert cosmological table in buffer.
Argument REDSHIFT redshift.
Argument H0 Hubble parameter today.
Argument OMATTER matter density parameter.
Argument OLAMBDA cosmological constant density parameter.
Argument OREL density parameter.
Argument HUBBLE Hubble parameter at given redshift.
Argument LOS-DIST line-of-sight comoving distance at given redshift."
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
          "\n")
  ;; Cosmological functions.
  (insert "Cosmography at required redshift\n"
          "--------------------------------\n"
          (format "- Hubble parameter [km/s/Mpc]:           %s\n"
                  hubble)
          (format "- Line-of-sight comoving distance [Mpc]: %s\n"
                  los-dist))
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
         (los-dist (cosmo-get-los-comoving-distance redshift)))
    (with-output-to-temp-buffer cosmo-buffer
      (pop-to-buffer cosmo-buffer)
      (cosmo--write-calc redshift H0 omatter olambda orel hubble los-dist))))


;;; Cosmopedia

(defun cosmo-pedia ()
  "Display a reference to basic cosmological definitions."
  (interactive)
  (let* ((cosmo-buffer "*Cosmopedia*"))
    (with-output-to-temp-buffer cosmo-buffer
      (pop-to-buffer cosmo-buffer)
      (insert
       "* Distances relations\n\n"
       "  - D_M :: Comoving distance (transverse).\n"
       "  - D_A = D_M / (1+z) :: Angular diameter distance.\n"
       "  - D_L = (1+z) D_M = (1+z)^2 D_A :: Luminosity distance.\n\n"
       "* Conversion factors, units\n\n"
       "  *Units system*: hbar = c = k_Boltzmann = 1.\n\n"
       "  - 1 pc = 3.2612 light years = 3.0856e18 cm\n"
       "  - 1 Mpc = 1e6pc ~ 3e24 cm ~ 1e14 s\n\n"
       "* Important constants\n\n"
       "  - Hubble constant :: H0 = 100h km/s/Mpc\n"
       "                          = 2.1332e-42 h GeV\n"
       "  - Hubble time, distance :: 1/H0 = 3.0856e17/h s\n"
       "                                  = 9.7776e9/h yr\n"
       "                                  = 2997.9/h Mpc\n"
       "                                  = 9.2503e27/h cm\n"))))

;;; Unit test.

(eval-when-compile
  (require 'cl)

  (defun cosmo-and-reduce (seq)
    "Reduce sequence SEQ by applying the `and` special form."
    (reduce #'(lambda (x y) (and x y)) seq))

  (defun cosmo-almost-eq (num1 num2 &optional reltol)
    ;; TODO: Decide how to handle the case in which the arguments are
    ;; zero; for instance, both a relative and absolute tolerance can
    ;; be specified and the largest one is considered.
    "Return t if the two numbers differ by less than ABSTOL."
    (or reltol (setq reltol 1e-8))       ; Default relative tolerance
      (> reltol (abs (1- (/ num1 num2)))))

  (defun cosmo-test-string-number-p ()
    "Test string representing numbers."
    (let ((numbers '("1" "+2" "-30" "1.2" "+30.4" "-5.60")))
      (assert
       (cosmo-and-reduce (mapcar #'cosmo--string-number-p numbers)))))

  (defun cosmo-test-string-notnumber-p ()
    "Test strings not representing numbers."
    (let ((notnumbers '("a" "+abc" "-30..0" "10.0..2")))
      (assert
       (not (cosmo-and-reduce (mapcar #'cosmo--string-number-p notnumbers))))))

  (defun cosmo-test-set-default ()
    "Set default test cosmological parameters."
        (puthash "H0 [Km/s/Mpc]" 70.0 cosmo--params)
        (puthash "omatter" 0.31 cosmo--params)
        (puthash "olambda" 0.7 cosmo--params)
        (puthash "orel" 8.52444340102e-05 cosmo--params)
        nil)

  (defun cosmo-test-efunc ()
    (assert (cosmo-almost-eq (cosmo-efunc 1000.0) 19912.4772226 1e-3)))

  (defun cosmo-test-inv-efunc ()
    (assert (cosmo-almost-eq (cosmo-inv-efunc 1000.0) 5.02197686819e-05 1e-3)))

  (defun cosmo-test-hubble ()
    (assert (cosmo-almost-eq (cosmo-get-hubble 1000.0) 1393873.40558 1e-3)))

  (defun cosmo-test-hubble-distance ()
    (assert (cosmo-almost-eq (cosmo-get-hubble-distance) 4282.7494 1e-3)))

  (defun cosmo-test-hubble-time ()
    (assert (cosmo-almost-eq (cosmo-get-hubble-time) 13.9684603096 1e-3)))

  (defun cosmo-test-los-comoving-distance ()
    (assert (cosmo-almost-eq
             (cosmo-get-los-comoving-distance 1000.0) 13454.7229832 1e-3)))

  (defun cosmo-test-transverse-comoving-distance ()
    (assert (cosmo-almost-eq
             (cosmo-get-transverse-comoving-distance 1000.0)
             13232.6210034 1e-3)))

  (defun cosmo-test-luminosity-distance ()
    (assert (cosmo-almost-eq
             (cosmo-get-luminosity-distance 1000.0) 13245853.6244 1e-3)))

  (defun cosmo-test-angular-diameter-distance ()
    (assert (cosmo-almost-eq
             (cosmo-get-angular-diameter-distance 1000.0) 13.2194016018 1e-3)))

  (cosmo-test-string-number-p)
  (cosmo-test-string-notnumber-p)
  (cosmo-test-set-default)
  (cosmo-test-efunc)
  (cosmo-test-inv-efunc)
  (cosmo-test-hubble)
  (cosmo-test-hubble-distance)
  (cosmo-test-hubble-time)
  (cosmo-test-los-comoving-distance)
  (cosmo-test-transverse-comoving-distance)
  (cosmo-test-luminosity-distance)
  (cosmo-test-angular-diameter-distance))

(provide 'cosmo)

;;; cosmo.el ends here
