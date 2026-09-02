(ns candle-j.noise-estimation
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

(defn- next-pow-2 ^long [^Integer val]
  (if (zero? val) 0 (loop [retval 1]
                      (if (and (< retval val) (not= val 0))
                        (recur (bit-shift-left retval 1))
                        retval))))

(defn get-padded-dimension [dimension]
  (next-pow-2 dimension))

(defn get-padded-dimensions [dimensions]
  (mapv get-padded-dimension dimensions))

(defn zero-padding [image]
  (let [dimensions (dtype/shape image)
        padded-dimensions (get-padded-dimensions dimensions)]
    (if (= dimensions padded-dimensions)
      image  ; Return original image if no padding needed
      (let [new-tensor (dtt/new-tensor padded-dimensions {:datatype (dtype/elemwise-datatype image)})]
        (dtt/tensor-copy! image (dtt/select new-tensor (range (dimensions 0)) (range (dimensions 1)) (range (dimensions 2))))
        new-tensor))))


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
    rotated-and-transposed-image))


(defn estimate-noise [image])
