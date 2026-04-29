# maass-form-core Makefile

build:
	sage -python setup.py build_ext --inplace $* 2>&1

sdist:
	sage -python setup.py sdist

install: clean sdist
	sage -pip install --no-build-isolation dist/maass_form_core-*.gz

test:
	sage -t src/*

tox:
	sage -pip install tox
	sage --python -m tox src -c tox.ini

clean:
	rm -rf src/*.c
	rm -rf src/*/*.c
	rm -rf src/*.so
	rm -rf *.egg-info
	rm -rf build
	rm -rf dist

.PHONY: build sdist install test clean tox
