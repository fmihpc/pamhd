all:
	+@make tests/divergence_executables tests/divergence_tests -C ../..

t: test
test:
	@$(MAKE) tests/divergence_tests -C ../..

c: clean
clean:
	@make tests/divergence_clean -C ../..
