;;; cosmo.el --- Cosmological Calculator    -*- lexical-binding: t; -*-
;;
;; Copyright (C) 2017 Francesco Montanari
;;
;;
;; Author: Francesco Montanari <fmnt@fmnt.info>
;; Created: 22 April 2017
(defconst cosmo:version "0.1") ;; Version:
;; Keywords: tools
;; Homepage: http://fmnt.info/software/cosmo_el.html
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
;;
;;
;;; Commentary:
;;
;; This program provides a cosmological calculator for Lambda-CDM
;; cosmological models.  Such a framework describes a homogeneous and
;; isotropic universe containing a cosmological constant (Lambda) and
;; a Cold Dark Matter (CDM) component, besides ordinary species.
;;
;; Several interactive commands are provided to set the cosmological
;; parameters and to compute cosmological functions at given redshift
;; values.  Definitions follow Hoggs (1999)
;; <https://arxiv.org/abs/astro-ph/9905116>.
;;
;;; Code:


;; Hash table containing all independent cosmological parameters
(defvar cosmo--params
  (let ((table (make-hash-table :test 'equal)))
    (puthash "H0 [Km/s/Mpc]" 70.0 table) ; Hubble today km/s/Mpc
    (puthash "omatter" 0.3 table)        ; Matter today
    (puthash "olambda" 0.7 table)        ; Curvature today
    (puthash "oradiation" 0.0 table)     ; Radiation today
    table)
  "Table containing Lambda-CDM cosmological parameters.")


;; Derived cosmological parameter
(defun cosmo--get-ocurvature ()
  "Get curvature density parameter today from Friedmann equations."
  (let ((omatter (gethash "omatter" cosmo--params))
        (olambda (gethash "olambda" cosmo--params))
        (oradiation (gethash "oradiation" cosmo--params)))
    (- 1.0 omatter olambda oradiation)))


;; Numerical utilities
(defun cosmo--sinh (x)
  "Hyperbolic sine of real arguments X."
  (* 0.5 (- (exp x) (exp (- x)))))


(defun cosmo--trapz (func a b &optional nstep)
  "Trapezoidal rule.

Integrate a function FUNC of one argument from A to B in NSTEP
equally spaced steps. The values A and B will be considered as
float.

Example:
\(cosmo--trapz '\(lambda \(x\) x\) 0.0 1.0\)"
  (let* ((a (float a))                  ; Extremes must be floats
         (b (float b))
         (nstep (or nstep 50))
         (step (/ (- b a) nstep))
         (sum 0.0))
    (setq sum (+ sum (* 0.5 (funcall func a))))
    (dotimes (i (- nstep 1))
      (setq sum (+ sum (funcall func (+ a (* step (1+ i)))))))
    (setq sum (+ sum (* 0.5 (funcall func b))))
    (* step sum)))


;; Read parameters
(defun cosmo--read-param (name)
  "Read parameter NAME from minibuffer and convert it to number."
  (string-to-number (read-from-minibuffer (format "Enter %s: " name))))


(defun cosmo--put-param (name)
  "Read parameter NAME from minibuffer and add it to the
parameter table."
  (puthash name (cosmo--read-param name) cosmo--params))


(defun cosmo--check-param (name value)
  "Check the validity of NAME (a cosmological parameter) VALUE."
  (cond ((or (string= name "omatter")
             (string= name "olambda")
             (string= name "oradiation"))
         (unless (>= value 0.0)
           (error "Error: density parameter must be positive")))))


(defun cosmo-set-params ()
  "Change the values of cosmological parameters."
  (interactive)
  (maphash (lambda (key value)
             (cosmo--put-param key)
             (cosmo--check-param key (gethash key cosmo--params)))
           cosmo--params))


;; (defun cosmo-set-default ()
;;   "Set cosmological parameters to the default values."
;;   (interactive)
;;     (clrhash cosmo--params)
;;     (puthash "H0 [Km/s/Mpc]" 70.0 cosmo--params)
;;     (puthash "omatter" 0.3 cosmo--params)
;;   nil)


;; Compute cosmological functions
(defun cosmo--efunc (redshift)
  "E(z) function at a given REDSHIFT."
  (let ((omatter (gethash "omatter" cosmo--params))
        (olambda (gethash "olambda" cosmo--params))
        (oradiation (gethash "oradiation" cosmo--params))
        (ocurvature (cosmo--get-ocurvature))
        (zp1 (+ 1 redshift)))
    (sqrt (+ (* oradiation (expt zp1 4.0))
             (* omatter (expt zp1 3.0))
             (* ocurvature (expt zp1 2.0))
             olambda))))


;; Compute cosmological functions
(defun cosmo--inv-efunc (redshift)
  "Inverse E(z) function at a given REDSHIFT."
  (cosmo--efunc redshift))


(defun cosmo--get-hubble (redshift)
  "Hubble parameter [Km/s/Mpc] for Lambda-CDM at a given
REDSHIFT."
  (let ((H0 (gethash "H0 [Km/s/Mpc]" cosmo--params))
        (zp1 (+ 1 redshift)))
    (* H0 (cosmo--efunc redshift))))


(defun cosmo-hubble ()
  "Display Hubble parameter in mini-buffer."
  (interactive)
  (let ((z (cosmo--read-param "redshift")))
    (message (format "%s km/s/Mpc" (cosmo--get-hubble z)))))


(defun cosmo--get-hubble-distance ()
  "Hubble distance [1/Mpc] for Lambda-CDM."
  (let ((H0 (gethash "H0 [Km/s/Mpc]" cosmo--params)))
    (/ 3.0e5 H0)))


;; (defun cosmo-hubble-distance ()
;;   "Display Hubble parameter in mini-buffer."
;;   (interactive)
;;   (message (format "%s 1/Mpc" (cosmo--get-hubble-distance))))


(defun cosmo--get-los-comoving-distance (redshift)
  "Line-of-sight comoving distance [Mpc] for Lambda-CDM at a
given REDSHIFT."
  (let ((DH (cosmo--get-hubble-distance))
        (int (cosmo--trapz 'cosmo--inv-efunc 0.0 redshift)))
    (* DH int)))


(defun cosmo-los-comoving-distance ()
  "Display line-of-sight comoving distance in mini-buffer."
  (interactive)
  (let ((z (cosmo--read-param "redshift")))
    (message (format "%s Mpc" (cosmo--get-los-comoving-distance z)))))


;; Write output
(defun cosmo--write-calc-header ()
  "Write header for the cosmological calculator summary buffer."
  (let ((head "Cosmology calculator.\n\n")
        (help (concat "(SPACE to scroll-down, BACKSPACE to scroll-up, "
                      "`q` to quit.)\n\n")))
    ;; (insert (propertize help 'font-lock-face 'italic))
    (insert head)))


(defun cosmo--write-calc (redshift H0 omatter olambda oradiation hubble)
  "Format and insert cosmological table in buffer.
Argument REDSHIFT redshift.
Argument H0 Hubble parameter today.
Argument OMATTER matter density parameter.
Argument OLAMBDA cosmological constant density parameter.
Argument ORADIATION density parameter.
Argument HUBBLE Hubble parameter at given redshift."
  ;; Input parameters
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
          (format "- Radiation fractional density, now:        %s\n"
                  oradiation)
          "\n")
  ;; Derived parameters
  (insert "Derived parameters\n"
          "------------------\n"
          (format "- Curvature fractional density: %s\n"
                  (cosmo--get-ocurvature))
          "\n")
  ;; Cosmological functions
  (insert "Cosmography at required redshift\n"
          "--------------------------------\n"
          (format "- Hubble parameter [km/s/Mpc]: %s\n"
                  hubble))
  nil)


(defun cosmo-calculator ()
  "Compute cosmology and display summary table in a new buffer."
  (interactive)
  (let* ((redshift (cosmo--read-param "redshift"))
         (omatter (gethash "omatter" cosmo--params))
         (olambda (gethash "olambda" cosmo--params))
         (oradiation (gethash "oradiation" cosmo--params))
         (H0 (gethash "H0 [Km/s/Mpc]" cosmo--params))
         (hubble (cosmo--get-hubble redshift)))
    (switch-to-buffer-other-window "*cosmo*")
    (erase-buffer)
    (cosmo--write-calc redshift H0 omatter olambda oradiation hubble)
    (beginning-of-buffer)
    ;; (special-mode) ; Problem: the buffer won't update, it must be
    ;;                ; killed.
    (other-window 1)))


(provide 'cosmo)


;;; cosmo.el ends here
