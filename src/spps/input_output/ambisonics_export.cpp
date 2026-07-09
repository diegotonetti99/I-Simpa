/* ----------------------------------------------------------------------
* I-SIMPA (https://i-simpa.univ-gustave-eiffel.fr). This file is part of I-SIMPA.
*
* I-SIMPA is a GUI for 3D numerical sound propagation modelling dedicated
* to scientific acoustic simulations.
* Copyright (C) UMRAE, CEREMA, Univ Gustave Eiffel - Judicael Picaut, Nicolas Fortin
*
* I-SIMPA is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 3 of the License, or
* (at your option) any later version.
*
* I-SIMPA is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software Foundation,
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA or
* see <http://ww.gnu.org/licenses/>
*
* For more information, please consult: <https://i-simpa.univ-gustave-eiffel.fr> or
* send an email to contact@noise-planet.org
*
* ----------------------------------------------------------------------*/

#include "ambisonics_export.h"
#include <fstream>
#include <cmath>
#include <cstdint>
#include <random>
#include <algorithm>

namespace AmbisonicsExport
{
	static const int SAMPLE_RATE = 44100;
	static const double W_SCALE = 0.70710678118654752440; // 1/sqrt(2), traditional B-format W scaling
	static const int NB_FILTER_SECTIONS = 3; // cascaded identical biquads, sharpens the bandpass skirts

	namespace
	{
		// Direct form I biquad, one section = 12dB/octave skirts around the bandpass center
		struct Biquad
		{
			double b0, b1, b2, a1, a2;
			double x1, x2, y1, y2;
			Biquad() : b0(1), b1(0), b2(0), a1(0), a2(0), x1(0), x2(0), y1(0), y2(0) {}

			inline float Process(float in)
			{
				double out = b0 * in + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
				x2 = x1; x1 = in;
				y2 = y1; y1 = out;
				return (float)out;
			}
		};

		// RBJ Audio EQ Cookbook constant-skirt-gain bandpass
		void DesignBandpass(double centerFreq, double bandwidthHz, double sampleRate, Biquad& bq)
		{
			double q = bandwidthHz > 0. ? std::max(0.1, centerFreq / bandwidthHz) : 10.0;
			double w0 = 2.0 * M_PI * centerFreq / sampleRate;
			double alpha = sin(w0) / (2.0 * q);
			double cosw0 = cos(w0);
			double b0 = alpha, b1 = 0.0, b2 = -alpha;
			double a0 = 1.0 + alpha, a1 = -2.0 * cosw0, a2 = 1.0 - alpha;
			bq.b0 = b0 / a0; bq.b1 = b1 / a0; bq.b2 = b2 / a0;
			bq.a1 = a1 / a0; bq.a2 = a2 / a0;
		}

		// One continuous band-limited noise carrier covering the whole IR duration
		std::vector<float> GenerateBandNoise(double centerFreq, double bandwidthHz, size_t nbSamples, std::mt19937& rng)
		{
			std::normal_distribution<double> dist(0.0, 1.0);
			std::vector<float> noise(nbSamples);
			for(size_t i = 0; i < nbSamples; i++)
				noise[i] = (float)dist(rng);

			Biquad reference;
			DesignBandpass(centerFreq, bandwidthHz, SAMPLE_RATE, reference);
			for(int section = 0; section < NB_FILTER_SECTIONS; section++)
			{
				Biquad stage = reference; // fresh filter state per section, same coefficients
				for(size_t i = 0; i < nbSamples; i++)
					noise[i] = stage.Process(noise[i]);
			}
			return noise;
		}

		bool WriteWavFloat4(const CoreString& path, int sampleRate, const std::vector<float>& W,
			const std::vector<float>& X, const std::vector<float>& Y, const std::vector<float>& Z)
		{
			std::ofstream ofs(path.c_str(), std::ios::binary);
			if(!ofs.is_open())
				return false;

			const size_t nbFrames = W.size();
			const uint16_t numChannels = 4;
			const uint16_t bitsPerSample = 32;
			const uint32_t byteRate = (uint32_t)sampleRate * numChannels * (bitsPerSample / 8);
			const uint16_t blockAlign = numChannels * (bitsPerSample / 8);
			const uint32_t dataSize = (uint32_t)(nbFrames * blockAlign);
			const uint32_t fmtChunkSize = 16;
			const uint32_t riffSize = 4 + (8 + fmtChunkSize) + (8 + dataSize);
			const uint16_t audioFormat = 3; // IEEE float

			ofs.write("RIFF", 4);
			ofs.write((const char*)&riffSize, 4);
			ofs.write("WAVE", 4);
			ofs.write("fmt ", 4);
			ofs.write((const char*)&fmtChunkSize, 4);
			ofs.write((const char*)&audioFormat, 2);
			ofs.write((const char*)&numChannels, 2);
			uint32_t sr = (uint32_t)sampleRate;
			ofs.write((const char*)&sr, 4);
			ofs.write((const char*)&byteRate, 4);
			ofs.write((const char*)&blockAlign, 2);
			ofs.write((const char*)&bitsPerSample, 2);
			ofs.write("data", 4);
			ofs.write((const char*)&dataSize, 4);

			for(size_t i = 0; i < nbFrames; i++)
			{
				float frame[4] = { W[i], X[i], Y[i], Z[i] };
				ofs.write((const char*)frame, sizeof(frame));
			}
			return ofs.good();
		}
	}

	bool WriteBFormatImpulseResponse(const CoreString& wavFilePath, const std::vector<t_BandContribution>& bands, int nbTimeStep, double timeStep, unsigned int randomSeed)
	{
		if(bands.empty() || nbTimeStep <= 0)
			return false;

		// Per-bin sample boundaries computed from cumulative rounding so bins never drift out of sync
		std::vector<size_t> binStart(nbTimeStep + 1);
		for(int i = 0; i <= nbTimeStep; i++)
			binStart[i] = (size_t)llround(i * timeStep * SAMPLE_RATE);
		const size_t nbSamples = binStart[nbTimeStep];
		if(nbSamples == 0)
			return false;

		std::vector<float> W(nbSamples, 0.f), X(nbSamples, 0.f), Y(nbSamples, 0.f), Z(nbSamples, 0.f);
		std::mt19937 rng(randomSeed != 0 ? randomSeed : 1u);

		// Sort bands by center frequency so edges can be derived from neighboring centers
		std::vector<size_t> order(bands.size());
		for(size_t i = 0; i < order.size(); i++)
			order[i] = i;
		std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return bands[a].centerFreq < bands[b].centerFreq; });

		for(size_t oi = 0; oi < order.size(); oi++)
		{
			const t_BandContribution& band = bands[order[oi]];
			const double logCenter = log(band.centerFreq);
			const double logPrev = oi > 0 ? log(bands[order[oi - 1]].centerFreq)
				: (oi + 1 < order.size() ? logCenter - (log(bands[order[oi + 1]].centerFreq) - logCenter) : logCenter - log(2.0));
			const double logNext = oi + 1 < order.size() ? log(bands[order[oi + 1]].centerFreq)
				: logCenter + (logCenter - logPrev);
			const double lowEdge = exp((logPrev + logCenter) / 2.0);
			const double highEdge = exp((logCenter + logNext) / 2.0);
			const double bandwidth = highEdge - lowEdge;

			std::vector<float> carrier = GenerateBandNoise(band.centerFreq, bandwidth, nbSamples, rng);

			for(int bin = 0; bin < nbTimeStep; bin++)
			{
				const size_t s0 = binStart[bin], s1 = binStart[bin + 1];
				if(s1 <= s0)
					continue;
				const double energy = band.energyPerBin[bin];
				if(energy <= 0. || !std::isfinite(energy))
					continue;

				double sumSq = 0.;
				for(size_t n = s0; n < s1; n++)
					sumSq += (double)carrier[n] * carrier[n];
				if(sumSq <= 0.)
					continue;
				// Scale this bin's burst so its energy exactly matches the tracked value
				const double gain = sqrt(energy / sumSq);

				const veci_t& dirVec = band.directionPerBin[bin];
				const double dirLen = dirVec.length();
				double ux = 0., uy = 0., uz = 0.;
				if(dirLen > 0.)
				{
					// Tracked intensity points away from the source (direction of particle travel);
					// B-format X/Y/Z encode direction *to* the perceived source, so negate.
					ux = -dirVec.x / dirLen; uy = -dirVec.y / dirLen; uz = -dirVec.z / dirLen;
				}

				for(size_t n = s0; n < s1; n++)
				{
					const double c = gain * carrier[n];
					W[n] += (float)(c * W_SCALE);
					X[n] += (float)(c * ux);
					Y[n] += (float)(c * uy);
					Z[n] += (float)(c * uz);
				}
			}
		}

		// Normalize to avoid an overly hot signal, keeping headroom below full scale
		float peak = 0.f;
		for(size_t n = 0; n < nbSamples; n++)
		{
			peak = std::max(peak, std::max(std::max(fabsf(W[n]), fabsf(X[n])), std::max(fabsf(Y[n]), fabsf(Z[n]))));
		}
		if(peak > 0.f)
		{
			// Raw physical energy units are tiny, so always rescale to target peak (not just clamp down),
			// otherwise the output is correctly shaped but silent.
			const float target = 0.99f;
			const float scale = target / peak;
			for(size_t n = 0; n < nbSamples; n++)
			{
				W[n] *= scale; X[n] *= scale; Y[n] *= scale; Z[n] *= scale;
			}
		}

		return WriteWavFloat4(wavFilePath, SAMPLE_RATE, W, X, Y, Z);
	}
}
