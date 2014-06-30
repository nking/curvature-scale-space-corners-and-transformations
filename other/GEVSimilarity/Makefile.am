SUBDIRS = \
	src/main/c++/algorithms/curves \
	src/test/c++/algorithms/curves 

## target for running the test
runtest: src/test/c++/algorithms/curves/GEVSimilarityParameters_test
	./src/test/c++/algorithms/curves/GEVSimilarityParameters_test

runtest0: src/test/c++/algorithms/curves/ParametersKey_test
	./src/test/c++/algorithms/curves/ParametersKey_test

runtest00: src/test/c++/algorithms/curves/ParametersEncoder_test
	./src/test/c++/algorithms/curves/ParametersEncoder_test

runtest000: src/test/c++/algorithms/curves/Sorter_test
	./src/test/c++/algorithms/curves/Sorter_test

runtest0000: src/test/c++/algorithms/curves/MinUnknownUniverseCover_test
	./src/test/c++/algorithms/curves/MinUnknownUniverseCover_test

runtests:
	$(MAKE) runtest; \
	$(MAKE) runtest0; \
	$(MAKE) runtest00; \
	$(MAKE) runtest000; \
	$(MAKE) runtest0000

valgrind: 
	valgrind --leak-check=yes --leak-check=full --read-var-info=yes --track-origins=yes --show-reachable=yes -v ./src/test/c++/algorithms/curves/GEVSimilarityParameters_test

COVERDIR = docs/coverage

$(COVERDIR):
	if [ ! -d $(COVERDIR) ]; then mkdir -p $(COVERDIR); fi

cov-clean:
	for p in $(SUBDIRS); do \
	   rm -f $$p/*.gcda $$p/*.gcno ; \
	done ; \
	rm -f *.info 

cov-clean-reports:
	rm -rf $(COVERDIR)/* 
    
coverage: $(COVERDIR)
	$(MAKE) clean && \
	$(MAKE) cov-clean && \
	$(MAKE) cov-clean-reports && \
	$(MAKE) CXXFLAGS="-fprofile-arcs -ftest-coverage" check && \
	$(MAKE) runtests && \
	$(MAKE) cov

cov:
	lcov -t 'GEVSimilarityCover' -o gcov0.info -c -d src/main/c++/algorithms/curves; \
	lcov -r gcov0.info /usr/include/\* -o gcov.info; \
	genhtml -o docs/coverage gcov.info
