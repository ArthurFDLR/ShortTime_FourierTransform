#ifndef FFT_H
#define FFT_H

#include <complex>
#include <vector>
#include <string>

namespace Fft {
	
	typedef std::vector<std::complex<float>> vec_complex;
	typedef std::vector<float> vec_real;
	typedef std::vector<vec_complex> mat_complex;

	void transformRadix2(vec_complex &vec);
	vec_complex dft(const vec_complex &signal, size_t vec_begin = 0, size_t vec_size = 0);
	vec_complex fft(const vec_complex &signal, size_t vec_begin = 0, size_t vec_size = 0);
	vec_complex fft(const vec_real &signal, size_t vec_begin = 0, size_t vec_size = 0);
	mat_complex stft_dft(vec_complex &vec, size_t window_size, size_t window_step);
	mat_complex stft_fft(vec_complex &vec, size_t window_size, size_t window_step);
	mat_complex stft_fft(vec_real &vec, size_t window_size, size_t window_step);

	void save_spectrogram(mat_complex& spectro, std::vector<std::string> infos, std::string name = "spectrogram.csv");
}

#endif