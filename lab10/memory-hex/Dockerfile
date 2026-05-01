FROM rocker/shiny:latest

LABEL maintainer="[Анисимов Тимур] <[tanisimov2006@gmail.com]>"
LABEL version="1.0"
LABEL description="Hex Memory Game - Shiny Application"

RUN apt-get update && apt-get install -y \
    --no-install-recommends \
    libssl-dev \
    libcurl4-gnutls-dev \
    curl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ENV _R_SHLIB_STRIP_=true

COPY Rprofile.site /etc/R


RUN install2.r --error --skipinstalled \
    shiny \
    glue

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -y language-pack-ru && \
    locale-gen ru_RU.UTF-8 && \
    dpkg-reconfigure locales

ENV LANGUAGE=ru_RU.UTF-8
ENV LANG=ru_RU.UTF-8
ENV LC_ALL=ru_RU.UTF-8

COPY global.R /srv/shiny-server/
COPY server.R /srv/shiny-server/
COPY ui.R /srv/shiny-server/
COPY modules/ /srv/shiny-server/modules/
COPY functions/ /srv/shiny-server/functions/
COPY www/ /srv/shiny-server/www/


COPY shiny-server.conf /etc/shiny-server/shiny-server.conf


WORKDIR /srv/shiny-server

RUN addgroup --system shiny && \
    adduser --system --ingroup shiny shiny && \
    chown -R shiny:shiny /srv/shiny-server

EXPOSE 3838

USER shiny

CMD ["/usr/bin/shiny-server"]