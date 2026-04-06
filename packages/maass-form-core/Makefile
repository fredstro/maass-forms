# maass-form-core Makefile

# Warn if no virtual environment is active
define check_venv
	@if [ -z "$$VIRTUAL_ENV" ]; then \
		echo ""; \
		echo "WARNING: No virtual environment is active."; \
		echo "It is recommended to run this in the 'venv_test' virtual environment."; \
		echo ""; \
		read -p "Continue without a virtual environment? [y/N] " answer; \
		case "$$answer" in \
			[yY]*) ;; \
			*) echo "Aborted."; exit 1 ;; \
		esac; \
	fi
endef

build:
	$(call check_venv)
	sage -python setup.py build_ext --inplace $* 2>&1

sdist:
	$(call check_venv)
	sage -python setup.py sdist

install: clean sdist
	$(call check_venv)
	sage -pip install --no-build-isolation dist/maass_form_core-*.gz

test:
	$(call check_venv)
	sage -t src/*

tox:
	$(call check_venv)
	sage -pip install tox
	sage --python -m tox src -c tox.ini

clean:
	$(call check_venv)
	rm -rf src/*.c
	rm -rf src/*/*.c
	rm -rf src/*.so
	rm -rf *.egg-info
	rm -rf build
	rm -rf dist

.PHONY: build sdist install test clean tox
