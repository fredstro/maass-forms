# Maass Forms Monorepo - Top-level Makefile
# Orchestrates builds and tests across all packages.


# Warn if no virtual environment is active (skip for docker targets)
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

PACKAGES = packages/maass-form-core packages/maass-forms-hilbert packages/maass-forms-klein

# Detect the best available Python/pip.
# Prefer sage -python/sage -pip when sage is available (needed for Cython builds).
# Fall back to plain python/pip for environments where sage is on the Python path
# but not installed as a standalone command (e.g. passagemath via pip).
SAGE := $(shell command -v sage 2>/dev/null)
ifdef SAGE
  PYTHON = sage -python
  PIP = sage -pip
else
  PYTHON = python
  PIP = pip
endif

.PHONY: all install pip-install build test tox clean lint format help

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

all: install  ## Build and install all packages

install:  ## Install all packages via sub-Makefiles (requires sage)
	cd packages/maass-form-core && make install
	cd packages/maass-forms-hilbert && make install
	cd packages/maass-forms-klein && make install

pip-install:  ## Install all packages via pip (auto-detects sage or plain python)
	$(call check_venv)
	@# Core is installed with build isolation so pip creates a temp env with
	@# passagemath (from build-system.requires) for Cython compilation.
	@# This also installs passagemath into the environment as a runtime dep.
	$(PIP) install packages/maass-form-core
	@# Hilbert and klein use --no-build-isolation since maass_form_core (needed
	@# for cimport at build time) is now installed but not available on PyPI.
	$(PIP) install --no-build-isolation packages/maass-forms-hilbert
	$(PIP) install --no-build-isolation packages/maass-forms-klein

build:  ## Build Cython extensions in all packages
	@for pkg in $(PACKAGES); do \
		echo "=== Building $$pkg ===" && \
		cd $$pkg && make build && cd ../..; \
	done

test:  ## Run tests in all packages
	@for pkg in $(PACKAGES); do \
		echo "=== Testing $$pkg ===" && \
		cd $$pkg && make test && cd ../..; \
	done

tox:  ## Run full tox test suite in all packages
	@for pkg in $(PACKAGES); do \
		echo "=== Tox $$pkg ===" && \
		cd $$pkg && make tox && cd ../..; \
	done

lint:  ## Run ruff check across all packages
	ruff check packages/

format:  ## Run ruff format across all packages
	ruff format packages/

clean:  ## Clean build artifacts in all packages
	@for pkg in $(PACKAGES); do \
		echo "=== Cleaning $$pkg ===" && \
		cd $$pkg && make clean && cd ../..; \
	done

docker:  ## Build Docker containers for all packages
	@for pkg in $(PACKAGES); do \
		if [ -f "$$pkg/Dockerfile" ]; then \
			echo "=== Docker build $$pkg ===" && \
			cd $$pkg && make docker && cd ../..; \
		fi \
	done
