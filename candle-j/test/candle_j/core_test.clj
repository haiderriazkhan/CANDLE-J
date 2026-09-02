(ns candle-j.core-test
  (:require [clojure.test :refer :all]
            [candle-j.noise-estimation :refer :all]
            [tech.v3.tensor :as dtt]
            [tech.v3.datatype.functional :as dfn]))

(deftest lp2-test
  (testing "If x is a power of 2, the function simply returns x. 
            Otherwise, the function returns the next greatest power of 2. Undefined for zero."
    (are [input expected] (= (get-padded-dimension (int input)) expected)
      0 0
      1 1
      2 2
      3 4
      5 8
      6 8
      7 8
      8 8
      511 512
      512 512
      513 1024
      1024 1024
      1025 2048
      2147483647 2147483648)))

(deftest padded-dimensions-test
  (testing "The function returns a map of padded dimensions, 
            given a vector of unpadded dimensions [depth, width, height]."
    (is (= [8 16 32]
           (get-padded-dimensions [5 9 17])))))

(deftest zero-clamping-test
  (testing "The function returns a matrix with all elements >= 0, by shifting the 
            matrix up by the minimum element, if the minimum element is < 0."
    (is (dfn/equals (dtt/->tensor (partition 2 (partition 2 (range 1 9)))) 
                    (zero-clamping (dtt/->tensor (partition 2 (partition 2 (range 1 9)))))))
     (is (dfn/equals (dtt/->tensor (partition 2 (partition 2 (range 0 8))))
                    (zero-clamping (dtt/->tensor (partition 2 (partition 2 (range -1 7)))))))
    (is (dfn/equals (dtt/->tensor [[[0 3] [4 5]] [[6 7] [8 1]]]) 
                    (zero-clamping (dtt/->tensor [[[-1 2] [3 4]] [[5 6] [7 0]]]))))))

(deftest zero-padding-test
  (testing "The function returns a matrix with padded dimensions, 
            given a matrix with unpadded dimensions [depth, width, height]."
    (let [input (dtt/->tensor
                 (partition 2 (partition 2 (range 1 9)))
                 {:datatype :float32})]
      (is (identical? input (zero-padding input))))
    (is (dfn/equals (dtt/->tensor [[[1 2 3 0] [4 5 6 0] [7 8 9 0] [0 0 0 0]] 
                                        [[10 11 12 0] [13 14 15 0] [16 17 18 0] [0 0 0 0]]
                                        [[19 20 21 0] [22 23 24 0] [25 26 27 0] [0 0 0 0]]
                                        [[0 0 0 0] [0 0 0 0] [0 0 0 0] [0 0 0 0]]] {:datatype :float32})
                (zero-padding (dtt/->tensor (partition 3 (partition 3 (range 1 28))) {:datatype :float32}))))))
