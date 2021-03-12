.PHONY: build

all: build run draw_spectro

build:
	g++ ./src/main.cpp -I./src/include ./src/include/fft.cpp -o ./build/main

run:
	./build/main.exe $(file)

draw_spectro:
	./venv/Scripts/python.exe ./src/draw_spectrogram.py
