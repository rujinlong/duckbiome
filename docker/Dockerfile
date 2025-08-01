FROM jinlongru/duckbiome_base:latest
COPY ./docker/renv.lock.prod renv.lock
RUN --mount=type=cache,id=renv-cache,target=/root/.cache/R/renv R -e 'renv::restore()'
COPY ./docker/duckbiome_*.tar.gz /app.tar.gz
RUN R -e 'remotes::install_local("/app.tar.gz",upgrade="never")'
RUN rm /app.tar.gz
COPY . /home/rstudio/
RUN R -e 'BiocManager::install(c("TreeSummarizedExperiment", "ANCOMBC", "ALDEx2", "DESeq2", "mia", "miaViz"))'
EXPOSE 8087
# RUN useradd --create-home --shell /bin/bash rstudio
# RUN echo "rstudio:rstudio" | chpasswd
# USER rstudio
# CMD R -e "options('shiny.port'=80,shiny.host='0.0.0.0');library(duckbiome);duckbiome::run_app()"
