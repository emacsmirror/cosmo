;;; test_cosmo.el --- Unit test for cosmo.el    -*- lexical-binding: t; -*-

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

;; This package provides an unit test for cosmo.el, a cosmological
;; calculator.

;;; Todo:

;; - Refactor tests, now they are too cumbersome and repetitive.

;;; Code:


(require 'cl)
(require 'cosmo)

;;; Test utilities.

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

(defun cosmo-test-set-default (H0 omatter olambda orel)
  "Set default test cosmological parameters.
Argument H0 Hubble parameter [Km/s/Mpc].
Argument OMATTER matter density parameter today.
Argument OLAMBDA cosmological constant density parameter.
Argument OREL relativistic density parameter today."
  (puthash "H0 [Km/s/Mpc]" H0 cosmo--params)
  (puthash "omatter" omatter cosmo--params)
  (puthash "olambda" olambda cosmo--params)
  (puthash "orel" orel cosmo--params)
  nil)

(defmacro cosmo-measure-time (&rest body)
  "Measure the time it takes to evaluate BODY."
  `(let ((time (current-time)))
     ,@body
     (message "%.06f seconds" (float-time (time-since time)))))

;;; Tests independent of cosmology.

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

(cosmo-test-string-number-p)
(cosmo-test-string-notnumber-p)

;;; Open cosmology.

(cosmo-test-set-default 70.0 0.31 0.7 8.52444340102e-05)

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

(defun cosmo-test-transverse-comoving-distance-open ()
  (assert (cosmo-almost-eq
           (cosmo-get-transverse-comoving-distance 1000.0)
           13232.6210034 1e-3)))

(defun cosmo-test-luminosity-distance ()
  (assert (cosmo-almost-eq
           (cosmo-get-luminosity-distance 1000.0) 13245853.6244 1e-3)))

(defun cosmo-test-angular-diameter-distance ()
  (assert (cosmo-almost-eq
           (cosmo-get-angular-diameter-distance 1000.0) 13.2194016018 1e-3)))

(defun cosmo-test-comoving-volume ()
  (assert (cosmo-almost-eq
           (cosmo-get-comoving-volume 1000.0) 1.00014515316e+13 3e-3)))

(cosmo-test-efunc)
(cosmo-test-inv-efunc)
(cosmo-test-hubble)
(cosmo-test-hubble-distance)
(cosmo-test-hubble-time)
(cosmo-test-los-comoving-distance)
(cosmo-test-transverse-comoving-distance-open)
(cosmo-test-luminosity-distance)
(cosmo-test-angular-diameter-distance)
(cosmo-test-comoving-volume)

;;; Close cosmology.

;; If the transverse comoving distance passes this test, and the other
;; functions already passed the open cosmology test, then they are
;; also consistent with a close cosmology. This is because different
;; curvature were coded in the transverse distance, and other
;; distances are defined as functions of this one.

(cosmo-test-set-default 70.0 0.27 0.7 8.52444340102e-05)

(defun cosmo-test-transverse-comoving-distance-close ()
  (assert (cosmo-almost-eq
           (cosmo-get-transverse-comoving-distance 1000.0)
           14832.933576 1e-3)))

(defun cosmo-test-comoving-volume ()
  (assert (cosmo-almost-eq
           (cosmo-get-comoving-volume 1000.0) 1.24288502093e+13 3e-3)))

(cosmo-test-transverse-comoving-distance-close)
(cosmo-test-comoving-volume)

;;; Flat cosmology.

;; If the transverse comoving distance passes this test, and the other
;; functions already passed the open cosmology test, then they are
;; also consistent with a close cosmology. This is because different
;; curvature were coded in the transverse distance, and other
;; distances are defined as functions of this one.

(cosmo-test-set-default 70.0 0.3 0.7 0.0)

(defun cosmo-test-transverse-comoving-distance-flat ()
  (assert (cosmo-almost-eq
           (cosmo-get-transverse-comoving-distance 1000.0) 13660.5292969 1e-3)))

(defun cosmo-test-comoving-volume ()
  (assert (cosmo-almost-eq
           (cosmo-get-comoving-volume 1000.0) 1.06780313213e+13 3e-3)))

(cosmo-test-transverse-comoving-distance-flat)
(cosmo-test-comoving-volume)


;;; Benchmarks

;; Benchmarks: numerical integral.
;; (cosmo-test-set-default 70.0 0.31 0.7 8.52444340102e-05) ; Open
;; (cosmo-measure-time (cosmo-get-angular-diameter-distance 0.1))
;; (cosmo-measure-time (cosmo-get-angular-diameter-distance 1000.0))

;; (cosmo-test-set-default 70.0 0.27 0.7 8.52444340102e-05) ; Close
;; (cosmo-measure-time (cosmo-get-angular-diameter-distance 0.1))
;; (cosmo-measure-time (cosmo-get-angular-diameter-distance 1000.0))

;; (cosmo-test-set-default 70.0 0.3 0.7 0.0) ; Flat
;; (cosmo-measure-time (cosmo-get-angular-diameter-distance 0.1))
;; (cosmo-measure-time (cosmo-get-angular-diameter-distance 1000.0))

;; Benchmarks: hyperbolic sine.
;; (cosmo-measure-time (cosmo-sinh 0.5))
;; (cosmo-measure-time (calc-eval "sinh(0.5)"))
;; (cosmo-measure-time (calcFunc-sinh '(float 5 -1)))
