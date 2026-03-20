# Makefile — spring_grid
# Supports Linux (sdl2-config) and macOS/Homebrew (Apple Silicon + Intel)
# Usage: make / make debug / make clean

CC     = gcc
TARGET = spring_grid
SRC    = spring_grid.c

# ── locate sdl2-config ────────────────────────────────────────────────────────
SDL2_CONFIG := $(shell \
    command -v sdl2-config 2>/dev/null || \
    ls /opt/homebrew/bin/sdl2-config 2>/dev/null || \
    ls /usr/local/bin/sdl2-config 2>/dev/null)

ifeq ($(SDL2_CONFIG),)
    $(error sdl2-config not found. Install SDL2: brew install sdl2  OR  apt install libsdl2-dev)
endif

SDL2_CFLAGS := $(shell $(SDL2_CONFIG) --cflags)
SDL2_LIBS   := $(shell $(SDL2_CONFIG) --libs)

# ── locate fftw3 ─────────────────────────────────────────────────────────────
FFTW_CFLAGS := $(shell pkg-config --cflags fftw3 2>/dev/null || \
    ([ -f /opt/homebrew/include/fftw3.h ] && echo -I/opt/homebrew/include) || \
    ([ -f /usr/local/include/fftw3.h   ] && echo -I/usr/local/include))

FFTW_LIBS   := $(shell pkg-config --libs fftw3 2>/dev/null || \
    ([ -f /opt/homebrew/lib/libfftw3.a ] && echo -L/opt/homebrew/lib -lfftw3) || \
    echo -lfftw3)

# ── flags ─────────────────────────────────────────────────────────────────────
CFLAGS = -std=c17 -Wall -Wextra -Wpedantic $(SDL2_CFLAGS) $(FFTW_CFLAGS)
LIBS   = $(SDL2_LIBS) $(FFTW_LIBS) -lm

RFLAGS = -O2
DFLAGS = -O0 -g3 -fsanitize=address,undefined

# ── targets ───────────────────────────────────────────────────────────────────
.PHONY: all debug clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(RFLAGS) $< -o $@ $(LIBS)

debug: $(SRC)
	$(CC) $(CFLAGS) $(DFLAGS) $< -o $(TARGET)_debug $(LIBS)

clean:
	rm -f $(TARGET) $(TARGET)_debug
