#!/bin/bash

if [ -z "$NBPORT" ]; then
  NBPORT=8888
fi
if [ -z "$TOX_ARGS" ]; then
  TOX_ARGS="doctest,coverage,pycodestyle,relint,codespell"
fi
if [ -n "$GIT_BRANCH" ]; then
  echo "PULL from $GIT_BRANCH"
  git pull -ff origin $GIT_BRANCH;
fi
make clean
sage setup.py install
case $1 in
    test)
      echo "Docker container running Sage doctests for the maass_forms_hilbert package."
      sage -t src
      ;;
    tox)
      echo "Docker container running tox with $TOX_ARGS"
      sage -pip install tox
      sage -python -m tox -c tox.ini src -e $TOX_ARGS
      ;;
    examples)
      echo "Docker container with Jupyter Notebook interface to run example notebooks."
      echo "NOTE: The Jupyter Notebook server is only accessible using the URL from outside the container."
      sage -n jupyter --no-browser --ip='*' --port=$NBPORT\
                            --notebook-dir=examples\
                            --Application.log_level=10\
                            --NotebookApp.custom_display_url=http://127.0.0.1:$NBPORT\
                            --NotebookApp.use_redirect_file=True\
                            --NotebookApp.allow_remote_access=True\
                            --NotebookApp.allow_root=True\
                            --NotebookApp.browser=x-www-browser
      ;;
    run)
      sage
      ;;
    shell)
      /bin/bash
      ;;
    *)
      echo "Unknown command $1. Exiting..."
      ;;
esac
