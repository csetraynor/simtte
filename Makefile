.PHONY: test test-ode test-legacy test-slow check document smoke

test:
	Rscript dev/run-tests.R

test-ode:
	Rscript dev/run-tests.R ode

test-legacy:
	Rscript dev/run-tests.R legacy

test-slow:
	Rscript dev/run-tests.R all --slow

check:
	Rscript dev/run-tests.R --check

document:
	Rscript -e 'devtools::document()'

smoke:
	Rscript dev/smoke-ode.R $(ARGS)
