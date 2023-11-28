(ns candle-j.data-copying)

(defn get-pixels-from-image [imgStack]
  (loop [i 1 v (transient [])]
    (if (<= i (.size imgStack))
      (recur (inc i) (conj! v (-> imgStack
                                  (.getProcessor i)
                                  (.getFloatArray))))
      (persistent! v))))
