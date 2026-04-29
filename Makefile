# Maass Forms Monorepo - Top-level Makefile
# Orchestrates builds and tests across all packages.

PACKAGES = packages/maass-form-core packages/maass-forms-hilbert packages/maass-forms-klein

.PHONY: all install build test tox clean lint format help

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

all: install  ## Build and install all packages

install:  ## Install all packages (core first, then domain packages)
	cd packages/maass-form-core && make install
	cd packages/maass-forms-hilbert && make install
	cd packages/maass-forms-klein && make install

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
