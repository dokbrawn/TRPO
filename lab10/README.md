```markdown
# Лабораторная работа №10  
## Docker — создание образов и запуск контейнеров

**Выполнил:** [Ваше ФИО]  
**Группа:** [Ваша группа]  
**Дата:** 01 мая 2026 г.

---

### Цель лабораторной работы

Изучить основные возможности Docker:
- создание Docker-образов с помощью `Dockerfile`;
- запуск контейнеров с пробросом портов;
- управление образами и контейнерами;
- работу с готовыми приложениями (Node.js и R Shiny).

---

## Часть 1. Работа с простым Node.js приложением (psweb)

### 1.1 Клонирование проекта из GitHub

```bash
git clone https://github.com/nigelpoulton/psweb.git
cd psweb
```

**Пояснение:**  
Команда `git clone` загружает исходный код веб-приложения, написанного на Node.js и Express. Это готовый пример из книги Nigel Poulton "Docker Deep Dive".

**Скриншот 1:** Клонирование репозитория psweb  
![Клонирование репозитория psweb](screenshots/git-clone.png)

### 1.2 Сборка Docker-образа

```bash
docker build -t example:latest .
```

**Подробная расшифровка команды:**
- `docker build` — команда для сборки образа из `Dockerfile`
- `-t example:latest` — **tag** (метка). Задаёт имя образа (`example`) и версию (`latest`)
- `.` — контекст сборки (Docker отправляет все файлы текущей директории на Docker daemon)

**Что происходит во время сборки:**
Docker последовательно выполняет все инструкции из `Dockerfile` (FROM, RUN, COPY, WORKDIR, CMD и т.д.), создавая слои образа. Каждый слой кэшируется, что ускоряет последующие сборки.

**Скриншот 2:** Процесс сборки Docker-образа `example:latest`  
![Сборка образа example:latest](screenshots/docker-build.png)

### 1.3 Просмотр созданных образов

```bash
docker images
```

**Скриншот 3:** Список всех Docker-образов (виден `example:latest`)  
![Просмотр Docker images](screenshots/docker-images.png)

---

## Часть 2. Запуск и тестирование контейнера

### 2.1 Запуск контейнера в фоновом режиме

```bash
docker run -d --name web -p 8081:8080 example:latest
```

**Подробная расшифровка флагов:**
- `-d` (`--detach`) — запускает контейнер в фоновом режиме (detached mode)
- `--name web` — присваивает контейнеру человекочитаемое имя `web`
- `-p 8081:8080` (`--publish`) — **публикация порта**. Пробрасывает порт 8081 хоста на порт 8080 внутри контейнера.  
  *Почему 8081?* Порт 8080 уже был занят другим контейнером (`bookings-gateway`).
- `example:latest` — имя образа, из которого создаётся контейнер

**Скриншот 4:** Запуск контейнера `web`  
![Запуск контейнера web](screenshots/docker-run.png)

### 2.2 Просмотр запущенных контейнеров

```bash
docker ps
```

**Пояснение:**  
Команда показывает только **активные** (работающие) контейнеры. Вывод включает ID, имя, образ, статус, порты и время работы.

**Скриншот 5:** Список активных контейнеров  
![docker ps — активные контейнеры](screenshots/docker-ps.png)

### 2.3 Проверка работы приложения в браузере

Переходим по адресу: **http://localhost:8081**

**Скриншот 6:** Работающее Node.js веб-приложение  
![Веб-приложение в браузере](screenshots/web-app.png)

### 2.4 Остановка и удаление ресурсов

```bash
docker stop web          # остановка контейнера
docker rm web            # удаление остановленного контейнера
docker rmi example:latest # удаление образа
```

**Скриншот 7:** Очистка ресурсов (остановка и удаление)  
![Очистка ресурсов — пример](screenshots/cleanup.png)

---

## Часть 3. Развёртывание R Shiny приложения (Hex Memory Game)

### 3.1 Клонирование проекта

```bash
cd ..
git clone https://github.com/dreamRs/memory-hex.git
cd memory-hex
```

**Скриншот 8:** Клонирование репозитория memory-hex  
![Клонирование memory-hex](screenshots/git-clone-memory-hex.png)

### 3.2 Просмотр Dockerfile

**Скриншот 9:** Содержимое файла `Dockerfile` для Shiny-приложения  
![Dockerfile для R Shiny](screenshots/dockerfile.png)

### 3.3 Сборка образа Shiny-приложения

```bash
docker build -t ggweb:latest .
```

**Особенности сборки:**
- Базовый образ: `rocker/shiny:latest`
- Установка системных зависимостей (`libcurl4-openssl-dev`, `libssl-dev`)
- Установка R-пакетов (`shiny`, `glue`)
- Копирование всех файлов приложения в директорию `/srv/shiny-server/`

**Скриншот 10:** Сборка образа `ggweb:latest`  
![Сборка ggweb:latest](screenshots/docker-build-ggweb.png)

### 3.4 Запуск контейнера с Shiny-приложением

```bash
docker run -d --name web2 -p 8082:3838 ggweb:latest
```

**Расшифровка:**
- `-p 8082:3838` — внешний порт хоста **8082** → внутренний порт Shiny-сервера **3838**

**Скриншот 11:** Запуск контейнера `web2`  
![Запуск web2](screenshots/docker-run-web2.png)

### 3.5 Проверка запущенных контейнеров

```bash
docker ps
```

**Скриншот 12:** Активные контейнеры (виден `web2`)  
![docker ps с web2](screenshots/docker-ps-web2.png)

### 3.6 Проверка работы Shiny-приложения

Открываем в браузере: **http://localhost:8082**

**Скриншот 13:** Работающая игра «Hex Memory Game»  
![Hex Memory Game — Shiny приложение](screenshots/shiny-app.png)

### 3.7 Финальная очистка ресурсов

```bash
docker stop web2
docker rm web2
docker rmi ggweb:latest
```

**Скриншот 14:** Очистка ресурсов Shiny-приложения  
![Очистка ggweb](screenshots/cleanup-web2.png)

---

## Итоги и выводы

В ходе лабораторной работы были успешно выполнены следующие задачи:

### Созданные Docker-образы:
- `example:latest` — Node.js веб-приложение
- `ggweb:latest` — R Shiny приложение (игра Memory Hex)

### Запущенные контейнеры:
- `web` — порт **8081** (Node.js)
- `web2` — порт **8082** (R Shiny)

### Основные команды Docker, изученные в работе:

| Команда                          | Назначение                                      | Ключевые флаги                  |
|----------------------------------|--------------------------------------------------|---------------------------------|
| `docker build -t name .`         | Сборка образа из Dockerfile                     | `-t`                            |
| `docker run -d --name ...`       | Запуск контейнера в фоне                        | `-d`, `--name`, `-p`            |
| `docker ps`                      | Просмотр запущенных контейнеров                 | —                               |
| `docker images`                  | Просмотр всех образов                           | —                               |
| `docker stop <name>`             | Остановка контейнера                            | —                               |
| `docker rm <name>`               | Удаление контейнера                             | —                               |
| `docker rmi <image>`             | Удаление образа                                 | —                               |

### Проблемы, с которыми пришлось столкнуться:

1. **Конфликт портов** — порт 8080 был занят. Решение: использование порта 8081.
2. **Длительная сборка** образа `ggweb:latest` — вызвана установкой R-пакетов и зависимостей. Это нормальное поведение.
3. Необходимость правильного проброса портов (`-p host:container`).

### Вывод

Docker значительно упрощает процесс развёртывания приложений, обеспечивая единообразие среды выполнения. Мы на практике освоили полный цикл: от клонирования кода до сборки образа, запуска контейнера и очистки ресурсов.

---

## Структура лабораторной папки

```
lab10/
├── README.md
├── psweb/
│   └── Dockerfile
├── memory-hex/
│   └── Dockerfile
└── screenshots/
    ├── 01-git-clone.png
    ├── 02-docker-build.png
    ├── 03-docker-images.png
    ├── 04-docker-run.png
    ├── 05-docker-ps.png
    ├── 06-web-app.png
    ├── 07-cleanup.png
    ├── 08-git-clone-memory-hex.png
    ├── 09-dockerfile.png
    ├── 10-docker-build-ggweb.png
    ├── 11-docker-run-web2.png
    ├── 12-docker-ps-web2.png
    ├── 13-shiny-app.png
    └── 14-cleanup-web2.png
```
