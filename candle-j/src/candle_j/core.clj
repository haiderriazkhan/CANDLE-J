(ns candle-j.core
  (:import [ij IJ ImageStack ImagePlus]
           [ij.process StackStatistics ]
           [ij.gui GenericDialog]
           [ij.plugin Filters3D]
           [ij.io FileSaver])
  (:require [candle-j.noise-estimation :refer [zero-padding]]
            [candle-j.data-copying :refer [get-pixels-from-image]]
            [criterium.core :as criterium]
            [tech.v3.tensor :as dtt]
            [tech.v3.datatype :as dtype]
            [tech.v3.datatype.statistics :as dts]
            [tech.v3.datatype.functional :as dtf]
            [clojure.core.matrix :as matrix]))

(def inputImage (IJ/openImage))
(def width (.getWidth inputImage))
(def height (.getHeight inputImage))
(def depth (.getNSlices inputImage))
(def dimensions (.getNDimensions inputImage))


(print width height depth dimensions)
(new StackStatistics inputImage)


(def filterParams (let [gd (new GenericDialog "Options")]
                    (.addMessage gd "Filter Parameters")
                    (.addNumericField gd "Smoothing Parameter" 0.1 2)
                    (.addNumericField gd "Patch radius" 2 1)
                    (.addNumericField gd "Search volume radius" 3 1)
                    (.showDialog gd)
                    (if (.wasCanceled gd) (IJ/log "User canceled dialog!")
                        {:beta (.getNextNumber gd)
                         :patch-radius (.getNextNumber gd)
                         :search-radius (.getNextNumber gd)})))

(print (filterParams :patch-radius))

(def medianFilteredImage (-> inputImage
                              (.getStack)
                              (Filters3D/filter (Filters3D/MEDIAN) 2 2 2)
                              (->> (new ImagePlus "Median Filtered Image"))))

(new StackStatistics medianFilteredImage)

(defn anscombe-transform [image] 
  (IJ/run image "32-bit" "")
  (IJ/run image "Add...", "value=0.375 stack")
  (IJ/run image "Square Root", "stack")
  (IJ/run image "Multiply...", "value=2 stack"))

(anscombe-transform inputImage)
(anscombe-transform medianFilteredImage)


(def inputImageStack (.getStack inputImage))
(def medianImageStack (.getStack medianFilteredImage))
(.getVoxel inputImageStack 270 200 6)

(new StackStatistics inputImage)

;; [z x y]
;; (def inputImageArray (dtt/->tensor (get-pixels-from-image inputImageStack) {:datatype :float64}))
;; (dtt/mget inputImageArray 6 270 200)

;; (criterium/quick-bench (zero-padding (dtt/->tensor (get-pixels-from-image inputImageStack) {:datatype :float64})))

;; (dtf/max image-array)

;; (criterium/quick-bench (dts/min (dtt/transpose (dtt/->tensor (get-pixels-from-image-2 inputImageStack)) [2 0 1])))
