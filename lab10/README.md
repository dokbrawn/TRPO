# Лабораторная работа №10  
## Docker — создание образов и запуск контейнеров

**Выполнил:** [Анисимов Тимур]  
**Группа:** [ИКС-431]  
**Дата:** 01 мая 2026 г.


### Цель лабораторной работы

Изучить основы работы с Docker: создание Docker-образов с помощью Dockerfile, запуск контейнеров, проброс портов, а также управление образами и контейнерами.


## Часть 1. Работа с Node.js приложением (psweb)

### 1.1 Клонирование проекта


git clone https://github.com/nigelpoulton/psweb.git
cd psweb

**Что делает команда:**  
`git clone` скачивает исходный код веб-приложения с GitHub в локальную папку `psweb`.

**Скриншот 1:** Клонирование репозитория psweb  
![Клонирование psweb](screenshots/git-clone.png)

### 1.2 Сборка Docker-образа

docker build -t example:latest .

**Что делает команда:**  
Создаёт Docker-образ на основе инструкций из файла `Dockerfile`, находящегося в текущей директории.

**Расшифровка флагов:**
- `-t example:latest` — присваивает образу имя `example` и тег `latest`
- `.` — указывает контекст сборки (Docker отправляет все файлы текущей папки в демон Docker)

**Скриншот 2:** Сборка Docker-образа `example:latest`  
![Сборка example:latest](screenshots/docker-build.png)

### 1.3 Просмотр списка образов

docker images

**Что делает команда:**  
Выводит список всех Docker-образов, хранящихся на компьютере.

**Скриншот 3:** Список Docker-образов  
![docker images](screenshots/docker-images.png)

### 1.4 Запуск контейнера

docker run -d --name web -p 8081:8080 example:latest


**Что делает команда:**  
Запускает новый контейнер из созданного образа `example:latest`.

**Расшифровка флагов:**
- `-d` (`--detach`) — запускает контейнер в фоновом режиме (не занимает терминал)
- `--name web` — присваивает контейнеру имя `web`
- `-p 8081:8080` (`--publish`) — публикует порт: открывает порт 8081 на хосте и связывает его с портом 8080 внутри контейнера
- `example:latest` — образ, из которого запускается контейнер

**Скриншот 4:** Запуск контейнера `web`  
![Запуск контейнера web](screenshots/docker-run.png)

### 1.5 Просмотр запущенных контейнеров

docker ps

**Что делает команда:**  
Показывает список всех **активных** (работающих) контейнеров.

**Скриншот 5:** Список активных контейнеров  
![docker ps](screenshots/docker-ps.png)

### 1.6 Проверка работы приложения

Открываем в браузере: **http://localhost:8081**

**Скриншот 6:** Работающее Node.js веб-приложение  
![Node.js приложение](screenshots/web-app.png)

### 1.7 Остановка и удаление ресурсов

docker rm -f web && docker rmi example:latest

**Что делают команды:**
- `docker rm -f web` — принудительно останавливает и удаляет контейнер `web`
- `docker rmi example:latest` — удаляет Docker-образ `example:latest`

**Скриншот 7:** Очистка ресурсов Node.js части  
![Очистка example](screenshots/cleanup.png)


## Часть 2. R Shiny приложение (Hex Memory Game)

### 2.1 Клонирование проекта

cd ..
git clone https://github.com/dreamRs/memory-hex.git
cd memory-hex

**Что делает команда:**  
Скачивает проект интерактивной игры на R Shiny.

**Скриншот 8:** Клонирование репозитория memory-hex  
![Клонирование memory-hex](screenshots/git-clone-memory-hex.png)

### 2.2 Просмотр Dockerfile

**Скриншот 9:** Содержимое файла Dockerfile  
![Dockerfile для Shiny приложения](screenshots/dockerfile.png)

### 2.3 Сборка образа Shiny-приложения

docker build -t ggweb:latest .

**Что делает команда:**  
Собирает Docker-образ `ggweb:latest` на основе Dockerfile проекта memory-hex (включая установку R-пакетов).

**Скриншот 10:** Сборка образа `ggweb:latest`  
![Сборка ggweb:latest](screenshots/docker-build-ggweb.png)

### 2.4 Запуск контейнера Shiny-приложения

docker run -d --name web2 -p 8082:3838 ggweb:latest

**Что делает команда:**  
Запускает контейнер из образа `ggweb:latest` в фоновом режиме.

**Расшифровка флагов:**
- `-d` — фоновый режим
- `--name web2` — имя контейнера
- `-p 8082:3838` — пробрасывает порт 8082 хоста на порт 3838 внутри контейнера (стандартный порт Shiny-сервера)

**Скриншот 11:** Запуск контейнера `web2`  
![Запуск web2](screenshots/docker-run-web2.png)

### 2.5 Просмотр запущенных контейнеров

docker ps

**Что делает команда:**  
Отображает список работающих контейнеров (включая `web2`).

**Скриншот 12:** Активные контейнеры (с web2)  
![docker ps с web2](screenshots/docker-ps-web2.png)

### 2.6 Проверка работы Shiny-приложения

Открываем в браузере: **http://localhost:8082**

**Скриншот 13:** Работающая игра Hex Memory Game  
![Hex Memory Game](screenshots/shiny-app.png)

### 2.7 Очистка ресурсов Shiny-приложения

docker rm -f web2 && docker rmi ggweb:latest

**Что делают команды:**
- `docker rm -f web2` — принудительно удаляет контейнер `web2`
- `docker rmi ggweb:latest` — удаляет образ `ggweb:latest`

**Скриншот 14:** Очистка ресурсов Shiny-приложения  
![Очистка ggweb](screenshots/cleanup-web2.png)


## Итоги и выводы

В рамках лабораторной работы были успешно:
- Клонированы два проекта с GitHub
- Создано два Docker-образа (`example:latest` и `ggweb:latest`)
- Запущены два контейнера (`web` и `web2`) с пробросом портов
- Выполнена очистка всех созданных ресурсов

### Основные изученные команды Docker:

| Команда                              | Назначение                                              | Основные флаги                     |
|--------------------------------------|---------------------------------------------------------|------------------------------------|
| `git clone <url>`                    | Скачивание проекта с GitHub                             | —                                  |
| `docker build -t name .`             | Сборка Docker-образа из Dockerfile                      | `-t`                               |
| `docker run -d --name NAME -p ...`   | Запуск контейнера в фоне с публикацией портов           | `-d`, `--name`, `-p`               |
| `docker ps`                          | Просмотр запущенных контейнеров                         | —                                  |
| `docker images`                      | Просмотр имеющихся образов                              | —                                  |
| `docker rm -f <name>`                | Принудительное удаление контейнера                      | `-f` (force)                       |
| `docker rmi <image>`                 | Удаление Docker-образа                                  | —                                  |

### Возникшие сложности и решения:
- Конфликт портов (8080 был занят) → использовали порт **8081**
- Долгая сборка образа `ggweb:latest` — нормальное явление из-за установки R-пакетов


## Структура папки lab10

```
lab10/
├── README.md
├── psweb/
├── memory-hex/
└── screenshots/
    ├── git-clone.png
    ├── docker-build.png
    ├── docker-images.png
    ├── docker-run.png
    ├── docker-ps.png
    ├── web-app.png
    ├── cleanup.png
    ├── git-clone-memory-hex.png
    ├── dockerfile.png
    ├── docker-build-ggweb.png
    ├── docker-run-web2.png
    ├── docker-ps-web2.png
    ├── shiny-app.png
    └── cleanup-web2.png
