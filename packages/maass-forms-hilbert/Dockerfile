ARG REMOTE_SRC=1
ARG GIT_BRANCH=develop

FROM  --platform=linux/x86_64 sagemath/sagemath:10.2 as base
RUN sudo apt-get update ds
RUN sudo apt-get -y install git

FROM base as use-git-1
RUN git clone https://github.com/fredstro/maass-forms-hilbert.git
WORKDIR "maass-forms-hilbert"
RUN git config pull.rebase false && git checkout ${GIT_BRANCH}

FROM base as use-git-0
ARG GIT_BRANCH=''
COPY --chown=sage . maass-forms-hilbert
WORKDIR "maass-forms-hilbert"

FROM use-git-${REMOTE_SRC} AS final
RUN sudo apt-get -y install make
COPY entrypoint.sh /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]
CMD ["run"]