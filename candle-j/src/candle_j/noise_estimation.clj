(ns candle-j.noise-estimation
  (:import [bit.math BitMath])
  (:require [tech.v3.tensor :as dtt]
            [tech.v3.datatype :as dtype]
            [tech.v3.datatype.functional :as dfn]
            [tech.v3.datatype.statistics :as stats]
            [tech.v3.datatype.convolve :as dt-conv]))

(def analysis-filter (dtype/->array
                      [0 0 -0.08838834764832 -0.08838834764832 0.69587998903400 -0.69587998903400 0.08838834764832 0.08838834764832 0.01122679215254 -0.01122679215254]))

(defn zero-clamping [image]
  (let [min (stats/min image)]
    (if (< min 0)
      (dfn/- image min)
      image)))

(defn get-padded-dimension [dimension]
  (BitMath/lp2 dimension))

(defn get-padded-dimensions [dimensions]
  (hash-map :padded-height (get-padded-dimension (nth dimensions 2))
            :padded-width (get-padded-dimension (nth dimensions 1))
            :padded-depth (get-padded-dimension (nth dimensions 0))))

(defn zero-padding [image]
  (let [dimensions (dtype/shape image)
        padded-dimensions (-> image
                              (dtype/shape)
                              (get-padded-dimensions))]
    (dtype/set-value!
     (dtt/clone (dtt/const-tensor 0
                                  [(:padded-depth padded-dimensions)
                                   (:padded-width padded-dimensions)
                                   (:padded-height padded-dimensions)])
                {:datatype :float64})
     [(range (dimensions 0)) (range (dimensions 1)) (range (dimensions 2))] image)))

(defn get-permutation-array [dimension]
  (case dimension
    0 [0 1 2]
    1 [2 0 1]
    2 [1 2 0]))

(defn get-inverse-permutation-array [dimension]
  (case dimension
    0 [0 1 2]
    1 [1 2 0]
    2 [2 0 1]))

(defn upfirdn [image]
  image)

(defn down-sample-image [image]
  image)

(defn analysis-filter-along-dimension [image dimension]
  (let [rotated-and-transposed-image (-> dimension
                                         (get-permutation-array)
                                         (->> (dtt/transpose image))
                                         (dtt/rotate [0 0 -5]))]
    rotated-and-transposed-image) 
  )

(defn estimate-noise [image])

(def test-1 (dtt/->tensor (partition 3 (partition 3 (range 1 28)))))


(def test-2 (dtt/->tensor [[[1 2 3 0] [4 5 6 0] [7 8 9 0] [0 0 0 0]] 
                           [[10 11 12 0] [13 14 15 0] [16 17 18 0][0 0 0 0]]
                           [[19 20 21 0] [22 23 24 0] [25 26 27 0][0 0 0 0]]
                           [[0 0 0 0] [0 0 0 0] [0 0 0 0][0 0 0 0]]]))

(dt-conv/convolve1d test-1 analysis-filter)
(dfn/equals (zero-padding test-1) test-2)
(dtt/const-tensor 0 [4 4 4])
