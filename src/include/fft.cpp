#include "fft.h"

#include <assert.h>
#include <fstream>

#define PI 3.14159265

static size_t reverseBits(size_t val, size_t width) {
	size_t result = 0;
	for (size_t i = 0; i < width; i++, val >>= 1)
		result = (result << 1) | (val & 1U);
	return result;
}

void Fft::transformRadix2(vec_complex &signal) {
	// Length variables
	size_t n = signal.size();
	size_t levels = 0;
	for (size_t temp = n; temp > 1U; temp >>= 1)
		levels++;
	assert(static_cast<size_t>(1U) << levels == n);

	// Trigonometric table
	vec_complex exp_table(n / 2);
	for (size_t i = 0; i < n / 2; i++)
		exp_table[i] = std::polar(1.0, -2. * PI * i / n);
	
	// Bit-reversed addressing permutation
	for (size_t i = 0; i < n; i++) {
		size_t j = reverseBits(i, levels);
		if (j > i)
			std::swap(signal[i], signal[j]);
	}
	
	// Cooley-Tukey decimation-in-time radix-2 FFT
	for (size_t size = 2; size <= n; size *= 2) {
		size_t half_size = size / 2;
		size_t tablestep = n / size;
		for (size_t i = 0; i < n; i += size) {
			for (size_t j = i, k = 0; j < i + half_size; j++, k += tablestep) {
				std::complex<float> temp = signal[j + half_size] * exp_table[k];
				signal[j + half_size] = signal[j] - temp;
				signal[j] += temp;
			}
		}
		if (size == n)  // Prevent overflow in 'size *= 2'
			break;
	}
}

Fft::vec_complex Fft::fft(const Fft::vec_complex &signal, size_t vec_begin, size_t vec_size) {
	Fft::vec_complex::const_iterator first = signal.begin() + vec_begin;
	Fft::vec_complex::const_iterator last = signal.begin() + vec_begin + vec_size;
	Fft::vec_complex subsignal(first, last);
	Fft::transformRadix2(subsignal);
	return subsignal;
}

Fft::vec_complex Fft::fft(const Fft::vec_real &signal, size_t vec_begin, size_t vec_size) {
	Fft::vec_real::const_iterator first = signal.begin() + vec_begin;
	Fft::vec_real::const_iterator last = signal.begin() + vec_begin + vec_size;
	Fft::vec_complex subsignal(first, last);
	Fft::transformRadix2(subsignal);
	return subsignal;
}

Fft::vec_complex Fft::dft(const Fft::vec_complex &signal, size_t vec_begin, size_t vec_size) {
	Fft::vec_complex output;
	size_t n = signal.size();
	if (vec_size == 0) vec_size = n;
	size_t vec_end = vec_begin + vec_size;
	assert(n >= vec_end);

	for (size_t k = 0; k < vec_size/2.; k++) {
		std::complex<float> sum = 0.;
		for (size_t t = 0; t < vec_size; t++) {
			float angle = 2 * PI * t * k / n;
			sum += signal[t+vec_begin] * std::exp(std::complex<float>(0, -angle));
		}
		output.push_back(sum);
	}
	return output;
}

Fft::mat_complex Fft::stft_dft(Fft::vec_complex &vec, size_t window_size, size_t window_step){
	Fft::mat_complex spectrogram = Fft::mat_complex();
	for (size_t begin = 0; begin < vec.size() - window_size; begin+=window_step)
	{
		spectrogram.push_back(Fft::dft(vec, begin, window_size));
	}
	return spectrogram;
}

Fft::mat_complex Fft::stft_fft(Fft::vec_complex &vec, size_t window_size, size_t window_step){
	Fft::mat_complex spectrogram = Fft::mat_complex();
	for (size_t begin = 0; begin < vec.size() - window_size; begin+=window_step)
	{
		spectrogram.push_back(Fft::fft(vec, begin, window_size));
	}
	return spectrogram;
}

Fft::mat_complex Fft::stft_fft(Fft::vec_real &vec, size_t window_size, size_t window_step){
	Fft::mat_complex spectrogram = Fft::mat_complex();
	for (size_t begin = 0; begin < vec.size() - window_size; begin+=window_step)
	{
		spectrogram.push_back(Fft::fft(vec, begin, window_size));
	}
	return spectrogram;
}

void Fft::save_spectrogram(mat_complex& spectro, std::string name){
	std::ofstream myfile;
	myfile.open(name);
	for (auto& row : spectro) { 
		for (auto& elem : row) { 
        	myfile << abs(elem) << ',';
		}
		myfile << '\n';
    }
	myfile.close();
}