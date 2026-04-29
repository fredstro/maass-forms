# SageMath is needed in path unless we use docker
ifeq (docker,$(findstring docker, $(MAKECMDGOALS)))
  DOCKER = 1
endif
ifeq (,$(DOCKER)$(shell which sage))
  $(error "No sage in $(PATH). Please install sage first or add the relevant path.")
endif
# Command line arguments:
# Set to 0 to use local source instead of GIT
REMOTE_SRC=1
# Set to desired git branch
GIT_BRANCH=main
# Port for notebook
NBPORT=8888
# Arguments for tox
TOX_ARGS=doctest,coverage,pycodestyle,relint,codespell,flake8
# Bind notebook example directory to a local directory
LOCAL_EXAMPLES=0

ifeq (0, $(LOCAL_EXAMPLES))
  BINDMOUNT=
else
  BINDMOUNT=-v $(LOCAL_EXAMPLES):/home/sage/maass_forms_klein/examples
endif
$(info $$LOCAL_EXAMPLES is [${LOCAL_EXAMPLES}])

$(info $$BINDMOUNT is [${BINDMOUNT}])
ifeq (0,$(REMOTE_SRC))
  TAG=local
  GIT_BRANCH=
else
  TAG=$(GIT_BRANCH)
endif

EXAMPLES_ARGS:=$(GIT_BRANCH) $(NBPORT)

build:
	sage -python setup.py build_ext --inplace $* 2>&1

sdist:
	sage -python setup.py sdist

sinstall: clean sdist
	sage -pip install --no-build-isolation dist/*.gz

install:
	sage -python setup.py install

test:
	sage -t src/*

examples:
	sage -n jupyter --no-browser --ip='0.0.0.0' --port=$(NBPORT)\
                            --notebook-dir=examples\
                            --NotebookApp.custom_display_url=http://127.0.0.1:$(NBPORT)\
                            --NotebookApp.use_redirect_file=False\
                            --NotebookApp.browser=x-www-browser

tox:
	sage -pip install tox meson
	sage --python -m tox src -c tox.ini -e $(TOX_ARGS)

docker:
	docker build --platform=linux/x86_64 --build-arg GIT_BRANCH=$(GIT_BRANCH) --build-arg REMOTE_SRC=$(REMOTE_SRC) -t maass_forms_klein-$(TAG) .

docker-rebuild:
	docker build --platform=linux/x86_64 --build-arg GIT_BRANCH=$(GIT_BRANCH) --build-arg REMOTE_SRC=$(REMOTE_SRC) --no-cache -t maass_forms_klein-$(TAG) .

docker-test: docker
	docker run --platform=linux/x86_64 -it -e GIT_BRANCH=$(GIT_BRANCH) --init maass_forms_klein-$(TAG) test

docker-examples: docker
	docker run --platform=linux/x86_64 -p $(NBPORT):$(NBPORT) -it -e GIT_BRANCH=$(GIT_BRANCH) -e NBPORT=$(NBPORT) $(BINDMOUNT) --init maass_forms_klein-$(TAG) examples $(EXAMPLES_ARGS)

docker-tox: docker
	docker run --platform=linux/x86_64 -it -e GIT_BRANCH=$(GIT_BRANCH) -e TOX_ARGS=$(TOX_ARGS) --init maass_forms_klein-$(TAG) tox

docker-shell: docker
	docker run --platform=linux/x86_64 -it -e GIT_BRANCH=$(GIT_BRANCH) --init maass_forms_klein-$(TAG) shell

docker-sage: docker
	docker run --platform=linux/x86_64 -it -e GIT_BRANCH=$(GIT_BRANCH) --init maass_forms_klein-$(TAG) run


clean:
	rm -rf src/maass_forms_klein/*/*.c
	rm -rf src/maass_forms_klein/*/*.so
	rm -rf src/maass_forms_klein/*/*.cpp
	rm -rf maass_forms_klein.egg-info
	rm -rf build
	rm -rf dist
	rm -rf Library
.PHONY: build sdist sinstall install test examples clean shell tox docker docker-examples docker-examples docker-rebuild \
		docker-shell docker-test docker-tox