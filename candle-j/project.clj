(defproject candle-j "0.1.0-SNAPSHOT"
  :description "FIXME: write description"
  :url "http://example.com/FIXME"
  :license {:name "EPL-2.0 OR GPL-2.0-or-later WITH Classpath-exception-2.0"
            :url "https://www.eclipse.org/legal/epl-2.0/"}

  :dependencies [[org.clojure/clojure "1.11.1"]
                 [net.imagej/ij "1.52g"]
                 [net.mikera/vectorz-clj "0.48.0"]
                 [criterium "0.4.6"]
                 [cnuernber/dtype-next "10.108"]]
  :java-source-paths ["java_src"]
  :repl-options {:init-ns candle-j.core})
