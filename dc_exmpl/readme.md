# Лабораторная работа №11
## Docker Compose: клиент‑серверный пример `dc_exampl` + BLAS Level 3 (`cblaslevel3`)

**Выполнил:** Анисимов Тимур  
**Группа:** ИКС-431  
**Дата:** 06 мая 2026 г.

---

## Цель работы

1. Освоить базовую работу с Docker Compose на примере клиент‑серверного приложения `dc_exampl`.
2. Подготовить и описать структуру проекта с отдельными контейнерами `server` и `client`.
3. Зафиксировать команды сборки/запуска и ожидаемые результаты.
4. Отразить структуру и новые файлы части `lab2(cblaslevel3)` с тестами BLAS Level 3.

---

## Часть 1. Пример `dc_exampl` (Docker Compose)

### 1.1 Структура проекта

Пример организуется в виде:

```text
dc_exampl/
├── docker-compose.yml
├── server/
│   ├── server.py
│   ├── index.html
│   └── Dockerfile
└── client/
    ├── client.py
    └── Dockerfile
```

### 1.2 Содержимое `server/server.py`

Сервер поднимает простой HTTP‑обработчик на порту `1234`:

```python
#!/usr/bin/env python3
import http.server
import socketserver

handler = http.server.SimpleHTTPRequestHandler

with socketserver.TCPServer(("", 1234), handler) as httpd:
    httpd.serve_forever()
```

### 1.3 Содержимое `server/index.html`

```html
Love it!
```

### 1.4 Содержимое `server/Dockerfile`

```dockerfile
FROM python:latest
WORKDIR /server/
COPY server.py index.html /server/
CMD ["python", "./server.py"]
```

### 1.5 Содержимое `client/client.py`

Клиент обращается к серверу и печатает ответ:

```python
#!/usr/bin/env python3
import urllib.request

fp = urllib.request.urlopen("http://localhost:1234/")
encodedContent = fp.read()
decodedContent = encodedContent.decode("utf8")
print(decodedContent)
fp.close()
```

### 1.6 Содержимое `client/Dockerfile`

```dockerfile
FROM python:latest
WORKDIR /client/
COPY client.py /client/
CMD ["python", "./client.py"]
```

### 1.7 Содержимое `docker-compose.yml`

Ключевые требования:
- объявить имя проекта и сервисы `server`/`client`;
- для `server` — проброс `1234:1234` и запуск `python ./server.py`;
- для `client` — запуск `python ./client.py`, сеть `host`, зависимость от `server`.

Пример:

```yaml
name: dc_exampl

services:
  server:
    build: ./server
    command: python ./server.py
    ports:
      - "1234:1234"

  client:
    build: ./client
    command: python ./client.py
    network_mode: host
    depends_on:
      - server
```

### 1.8 Команды запуска

```bash
docker compose build
docker compose up
```

Проверка результата в браузере:

```text
http://localhost:1234/
```

Ожидаемый результат: страница с текстом `Love it!`, а клиент в логах печатает HTML‑ответ.

---

## Часть 2. `lab2(cblaslevel3)`: добавленные и актуальные файлы

В каталоге лабораторной по BLAS Level 3 сейчас присутствуют реализация, инфраструктура тестирования и набор новых тестовых файлов.

### 2.1 Основные файлы сборки и запуска

- `Dockerfile` — контейнерная сборка проекта.
- `docker-compose.yml` — оркестрация контейнера для запуска сборки/тестов.
- `CMakeLists.txt` — конфигурация сборки CMake.
- `readme.md` — текущий отчет.

### 2.2 Заголовочные файлы (`include/`)

- `include/blas_types.h`
- `include/syrk.h`
- `include/test_macros.h`
- `include/test_utils.h`

### 2.3 Исходные файлы (`src/`)

Реализация и исполняемые модули:

- `src/syrk.c`
- `src/main.c`
- `src/benchmark_syrk.c`

Новые/актуальные тесты BLAS Level 3:

- `src/test_syrk.c`
- `src/test_syrk_herk.c`
- `src/test_syr2k_her2k.c`
- `src/test_gemm.c`
- `src/test_symm_hemm.c`
- `src/test_trmm_trsm.c`

### 2.4 Вспомогательные и результатные файлы

- `lib/mock_openblas.c` — заглушка/мок для сценариев сравнения.
- `results/log.txt` — текстовый лог запусков.
- `results/results.db` — база результатов тестирования.
- `artifacts/test_syrk.log`, `artifacts/benchmark_syrk.log` — журналы.
- `artifacts/test_syrk.svg`, `artifacts/benchmark_syrk.svg` — скриншоты/графические артефакты.

---

## Часть 3. Проверка `lab2(cblaslevel3)`

Пример локального запуска:

```bash
cmake -S "lab2(cblaslevel3)" -B "lab2(cblaslevel3)/build"
cmake --build "lab2(cblaslevel3)/build" -j
ctest --test-dir "lab2(cblaslevel3)/build" --output-on-failure
```

Запуск бенчмарка:

```bash
"lab2(cblaslevel3)/build/benchmark_syrk"
```

---
