FROM rocker/shiny:latest

LABEL maintainer="[Анисимов Тимур ИКС-431] <[tanisimov2006@gmail.com]>"
LABEL description="Hex Memory Game - Shiny Application"
LABEL version="1.0"

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('shiny', 'glue'), repos='https://cloud.r-project.org/')"

COPY global.R /srv/shiny-server/
COPY server.R /srv/shiny-server/
COPY ui.R /srv/shiny-server/
COPY modules/ /srv/shiny-server/modules/
COPY functions/ /srv/shiny-server/functions/
COPY www/ /srv/shiny-server/www/

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]