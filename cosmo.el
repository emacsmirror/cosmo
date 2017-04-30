;;; cosmo.el --- Cosmological Calculator
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
;; a Cold Dark Matter (CDM) component, besides ordinary
;; species.
;;
;; Several interactive commands are provided to set the cosmological
;; parameters and to compute cosmological functions at given redshift
;; values.
;;
;;; Code:


(require 'subr-x)


;; Hash table containing all independent cosmological parameters
(defvar cosmo--params
  (let ((table (make-hash-table :test 'equal)))
    (puthash "H0 [Km/s/Mpc]" 70.0 table) ; Hubble today km/s/Mpc
    (puthash "omatter" 0.3 table)        ; Matter today
    table)
  "Table containing LCDM cosmological parameters.")


;; Derived cosmological parameter
(defun cosmo--get-olambda ()
  "Get cosmological constant density parameter according to flat LCDM."
  (- 1. (gethash "omatter" cosmo--params)))


;; (defun cosmo-set-default ()
;;   "Set cosmological parameters to the default values."
;;   (interactive)
;;     (clrhash cosmo--params)
;;     (puthash "H0 [Km/s/Mpc]" 70.0 cosmo--params)
;;     (puthash "omatter" 0.3 cosmo--params)
;;   nil)


(defun cosmo--read-param (name)
  "Read parameter NAME from minibuffer and convert it to number."
  (string-to-number (read-from-minibuffer (format "Enter %s: " name))))


(defun cosmo--check-param (name value)
  "Check the validity of NAME (a cosmological parameter) VALUE."
  (cond ((string= name "omatter")
         (unless (and (> value 0.0) (< value 1.0))
           (error "Error: omatter must be positive and less than 1")))))


(defun cosmo--put-param (name)
  "Read parameter NAME from minibuffer and add it to the parameter table."
  (puthash name (cosmo--read-param name) cosmo--params))


(defun cosmo-set-params ()
  "Change the values of cosmological parameters."
  (interactive)
  (let ((params (hash-table-keys cosmo--params)))
    (dolist (param params)
      (cosmo--put-param param)
      (cosmo--check-param param (gethash param cosmo--params)))))


(defun cosmo--get-hubble (redshift)
  "Compute Hubble parameter for flat ΛCDM at a given REDSHIFT."
  (let ((omatter (gethash "omatter" cosmo--params))
        (olambda (cosmo--get-olambda))
        (H0 (gethash "H0 [Km/s/Mpc]" cosmo--params))
        (zp1 (+ 1 redshift)))
    (* H0 (sqrt (+ (* omatter (expt zp1 3)) olambda)))))


(defun cosmo-hubble ()
  "Display Hubble parameter in mini-buffer."
  (interactive)
  (let ((z (cosmo--read-param "redshift")))
    (message (format "%s km/s/Mpc" (cosmo--get-hubble z)))))


(defun cosmo--write-calc (redshift H0 omatter hubble)
  "Format and insert cosmological table in buffer.
Argument REDSHIFT redshift.
Argument H0 Hubble parameter today.
Argument OMATTER matter density parameter.
Argument HUBBLE Hubble parameter at given redshift."
  ;; Input parameters
  (insert "Cosmology calculator.\n\n"
	  "Input Parameters\n"
	  "----------------\n"
	  (format "- Redshift:                       \t%s\n"
		  redshift)
	  (format "- Hubble constant, now [km/s/Mpc]:\t%s\n"
		  H0)
	  (format "- Matter fractional density, now: \t%s\n"
		  omatter)
	  "\n")
  ;; Derived parameters
  (insert "Derived parameters\n"
	  "------------------\n"
	  (format "- Cosmological constant fractional density: \t%s\n"
		  (cosmo--get-olambda))
	  "\n")
  ;; Cosmological functions
  (insert "Cosmography at required redshift\n"
	  "--------------------------------\n"
	  (format "- Hubble parameter [km/s/Mpc]:\t%s\n"
		  hubble))
  nil)


(defun cosmo-calculator ()
  "Compute cosmology and display summary table in a new buffer."
  (interactive)
  (let* ((redshift (cosmo--read-param "redshift"))
         (omatter (gethash "omatter" cosmo--params))
         (H0 (gethash "H0 [Km/s/Mpc]" cosmo--params))
         (hubble (cosmo--get-hubble redshift)))
    (switch-to-buffer-other-window "*cosmo*")
    (erase-buffer)
    (cosmo--write-calc redshift H0 omatter hubble)
    (beginning-of-buffer)
    (other-window 1)))


(provide 'cosmo)

;;; cosmo.el ends here
