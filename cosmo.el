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


(require 'subr-x)


;; Define in this table all independent cosmological parameters
(defconst cosmo--params (make-hash-table :test 'equal)
  "Table containing LCDM cosmological parameters")

;; For flat LCDM olambda is a derived parameter
(defconst cosmo--olambda 0.7
  "Cosmological constant density parameter.")


(defun cosmo-set-default ()
  (interactive)
  "Set cosmological parameters table to default values"
  (clrhash cosmo--params)
  (puthash "H0 [Km/s/Mpc]" 70.0 cosmo--params) ; Hubble today km/s/Mpc
  (puthash "omatter" 0.3 cosmo--params)        ; Matter today
  )


(defun cosmo--read-param (name)
  "Read parameter from minibuffer and convert it to number."
  (string-to-number (read-from-minibuffer (format "Enter %s: " name))))


(defun cosmo--set-olambda ()
  "Set the cosmological constant density parameter for flat LCDM."
  (setq cosmo--olambda (- 1.0 (gethash "omatter" cosmo--params))))


(defun cosmo--check-param (name value)
  "Check the validity of a cosmological parameter value."
  (cond ((string= name "omatter")
         (unless (and (> value 0.0) (< value 1.0))
           (error "Error: omatter must be positive and less than 1")))))


(defun cosmo--put-param (name)
  "Read parameter from minibuffer and add it to the parameter table."
  (puthash name (cosmo--read-param name) cosmo--params))


(defun cosmo-set-params ()
  "Set cosmological parameters."
  (interactive)
  (let ((params (hash-table-keys cosmo--params)))
    (dolist (param params)
      (cosmo--put-param param)
      (cosmo--check-param param (gethash param cosmo--params)))))


(defun cosmo--get-hubble (redshift)
  "Compute Hubble parameter for flat ΛCDM."
  (let ((omatter (gethash "omatter" cosmo--params))
        (H0 (gethash "H0 [Km/s/Mpc]" cosmo--params))
        (zp1 (+ 1 redshift)))
    (* H0 (sqrt (+ (* omatter (expt zp1 3)) (- 1 omatter))))))


(defun cosmo-hubble ()
  "Display Hubble parameter in mini-buffer."
  (interactive)
  (let ((z (cosmo--read-param "redshift")))
    (message (format "%s km/s/Mpc" (cosmo--get-hubble z)))))


(defun cosmo--write-calc (redshift H0 omatter hubble)
  "Format and insert cosmological table in buffer."
  (insert "Cosmology calculator.\n\n"
          "Input Parameters\n"
          "----------------\n"
          (format "- Redshift:                       \t%s\n"
                  redshift)
          (format "- Hubble constant, now [km/s/Mpc]:\t%s\n"
                  H0)
          (format "- Matter fractional density, now: \t%s\n"
                  omatter)
          "\n"
          "Cosmography at required redshift\n"
          "-------------------------------\n"
          (format "- Hubble parameter [km/s/Mpc]:\t%s\n"
                  hubble)))


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
    (other-window 1)))


;; Set default values when loading the package
(cosmo-set-default)

(provide 'cosmo)
