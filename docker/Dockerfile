FROM jinlongru/duckbiome_base:latest
# COPY ./docker/renv.lock.prod renv.lock
# RUN R -e 'options(renv.config.pak.enabled = FALSE);renv::restore()'
COPY ./docker/duckbiome_*.tar.gz /app.tar.gz
RUN R -e 'remotes::install_local("/app.tar.gz",upgrade="never")'
RUN rm /app.tar.gz
EXPOSE 8087
# RUN useradd --create-home --shell /bin/bash rstudio
# RUN echo "rstudio:rstudio" | chpasswd
# USER rstudio
# CMD R -e "options('shiny.port'=80,shiny.host='0.0.0.0');library(duckbiome);duckbiome::run_app()"
