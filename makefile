.PHONY: build

all: build run draw_spectro

build:
	g++ ./src/compute_spectrogram.cpp -I./src/include ./src/include/fft.cpp -o ./build/compute_spectrogram

run:
	./build/compute_spectrogram.exe $(file)

draw_spectro:
	./venv/Scripts/python.exe ./src/draw_spectrogram.py
