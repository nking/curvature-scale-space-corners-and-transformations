SUBDIRS = \
	src/main/c++/algorithms/curves \
	src/test/c++/algorithms/curves 

## target for running the test
runtest: src/test/c++/algorithms/curves/GEVSimilarityParameters_test
	./src/test/c++/algorithms/curves/GEVSimilarityParameters_test

valgrind: all 
        valgrind --leak-check=yes --leak-check=full --read-var-info=yes --track-origins=yes --show-reachable=yes -v ./src/test/c++/algorithms/curves/GEVSimilarityParameters_test
