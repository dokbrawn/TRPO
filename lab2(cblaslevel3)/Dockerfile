FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    libopenblas-dev \
    libopenblas-openmp-dev \
    sqlite3 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY . .

# создаём папки (если volume пустой)
RUN mkdir -p /app/build /app/results

CMD bash -c "\
echo '=====================================' && \
echo 'SYRK BLAS - Сборка и тестирование' && \
echo '=====================================' && \

echo 'Сборка...' && \
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release && \
cmake --build build -j$(nproc) && \

echo 'Инициализация БД...' && \
sqlite3 /app/results/results.db \"CREATE TABLE IF NOT EXISTS results (\
id INTEGER PRIMARY KEY AUTOINCREMENT,\
test_name TEXT,\
output TEXT,\
created_at DATETIME DEFAULT CURRENT_TIMESTAMP);\" && \

echo 'Запуск тестов...' && \
FUNC_OUT=\$(./build/test_syrk) && \
BENCH_OUT=\$(./build/benchmark_syrk 5 512 512) && \

echo 'Сохранение логов...' && \
echo '===== FUNCTIONAL TEST =====' > /app/results/log.txt && \
echo \"\$FUNC_OUT\" >> /app/results/log.txt && \
echo '' >> /app/results/log.txt && \
echo '===== BENCHMARK =====' >> /app/results/log.txt && \
echo \"\$BENCH_OUT\" >> /app/results/log.txt && \

sqlite3 /app/results/results.db \"INSERT INTO results (test_name, output) VALUES ('functional', '\$FUNC_OUT');\" && \
sqlite3 /app/results/results.db \"INSERT INTO results (test_name, output) VALUES ('benchmark', '\$BENCH_OUT');\" && \

echo '=====================================' && \
echo 'build/ и results/ сохранены на хосте' && \
echo '=====================================' \
"