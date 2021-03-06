#include <iostream>
#include <assert.h>
#include <string>
#include <chrono>

#include "fft.h"
#include "AudioFile.h"

#define PI 3.14159265

int main(int argc, char** argv) 
{
	assert(argc==2 && "Please provide a WAV file path to analyze.");

	std::string audio_file_path(argv[1]);
	AudioFile<float> a;
	bool loadedOK = a.load(audio_file_path);
	assert (loadedOK);
	Fft::vec_real sample = a.samples[0]; //Analyze first channel

	size_t F_s = a.getSampleRate();
	double delta_t = 1./F_s;
	size_t N = sample.size();
	size_t window_size = pow(2, 10);
	size_t window_step = window_size/2;
	double delta_f = double(F_s)/window_size;

	auto t_start = std::chrono::high_resolution_clock::now();
	Fft::mat_complex spectrogram = Fft::stft_fft(sample, window_size, window_step, true);
	auto t_end = std::chrono::high_resolution_clock::now();
	uint64_t execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

	std::vector<std::string> infos = {
		"delta_t:" + std::to_string(double(window_step)/F_s),
		"delta_f:" + std::to_string(delta_f),
		"path:" + audio_file_path,
		"execution time (ms):" + std::to_string(execution_time),
	};
	Fft::save_spectrogram(spectrogram, infos ,"./spectrograms/spectrogram.csv");

	std::cout << "\t- Spectrogram successfully computed." << std::endl;
    return 0;
}