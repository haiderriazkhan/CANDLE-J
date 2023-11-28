(ns candle-j.core-test
  (:require [clojure.test :refer :all]
            [candle-j.noise-estimation :refer :all]
            [tech.v3.tensor :as dtt]
            [tech.v3.datatype.functional :as dfn]))

(deftest lp2-test
  (testing "If x is a power of 2, the function simply returns x. 
            Otherwise, the function returns the next greatest power of 2."
    (is (= 0 (get-padded-dimension 0)))
    (is (= 1 (get-padded-dimension 1)))
    (is (= 2 (get-padded-dimension 2)))
    (is (= 4 (get-padded-dimension 3)))
    (is (= 8 (get-padded-dimension 5)))
    (is (= 8 (get-padded-dimension 6)))
    (is (= 8 (get-padded-dimension 7)))
    (is (= 8 (get-padded-dimension 8)))
    (is (= 512 (get-padded-dimension 511)))
    (is (= 512 (get-padded-dimension 512)))
    (is (= 1024 (get-padded-dimension 513)))
    (is (= 1024 (get-padded-dimension 1024)))
    (is (= 2048 (get-padded-dimension 1025)))
    (is (= 0 (get-padded-dimension 2147483647)))))

(deftest padded-dimensions-test
  (testing "The function returns a map of padded dimensions, 
            given a vector of unpadded dimensions [depth, height, width]."
    (is (= {:padded-width 16
            :padded-height 16
            :padded-depth 8}
           (get-padded-dimensions [5 9 16])))))

(deftest zero-clamping-test
  (testing "The function returns a matrix with all elements >= 0, by shifting the 
            matrix up by the minimum element, if the minimum element is < 0."
    (is (dfn/equals (dtt/->tensor [[[1 2] [3 4]] [[5 6] [7 0]]]) 
                (zero-clamping (dtt/->tensor [[[1 2] [3 4]] [[5 6] [7 0]]]))))
    (is (dfn/equals (dtt/->tensor [[[0 3] [4 5]] [[6 7] [8 1]]]) 
                (zero-clamping (dtt/->tensor [[[-1 2] [3 4]] [[5 6] [7 0]]]))))))

(deftest zero-padding-test
  (testing "The function returns a matrix with padded dimensions, 
            given a matrix with unpadded dimensions [depth, height, width]."
    (is (dfn/equals (dtt/->tensor [[[1 2] [3 4]] [[5 6] [7 8]]])
                (zero-padding (dtt/->tensor [[[1 2] [3 4]] [[5 6] [7 8]]]))))
    (is (dfn/equals (dtt/->tensor [[[1 2 3 0] [4 5 6 0] [7 8 9 0] [0 0 0 0]] 
                                        [[10 11 12 0] [13 14 15 0] [16 17 18 0] [0 0 0 0]]
                                        [[19 20 21 0] [22 23 24 0] [25 26 27 0] [0 0 0 0]]
                                        [[0 0 0 0] [0 0 0 0] [0 0 0 0] [0 0 0 0]]])
                (zero-padding (dtt/->tensor [[[1 2 3] [4 5 6] [7 8 9]]
                                             [[10 11 12] [13 14 15] [16 17 18]]
                                             [[19 20 21] [22 23 24] [25 26 27]]]))))))
