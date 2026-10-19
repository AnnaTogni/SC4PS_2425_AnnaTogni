#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TLine.h>
#include <TMarker.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TH1.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace {

constexpr double kPi = 3.14159265358979323846;

double SignalValue(double t, const std::vector<double>& freqs) {
  double value = 0.0;
  for (double f : freqs) {
    value += std::sin(2.0 * kPi * f * t);
  }
  return value;
}

double AliasFrequency(double f, double fs) {
  if (fs <= 0.0) {
    return 0.0;
  }

  double wrapped = std::fmod(std::fabs(f), fs);
  if (wrapped < 0.0) {
    wrapped += fs;
  }
  if (wrapped > fs / 2.0) {
    wrapped = fs - wrapped;
  }
  return wrapped;
}

bool IsPowerOfTwo(int n) {
  return n > 0 && (n & (n - 1)) == 0;
}

int NextPowerOfTwo(int n) {
  int power = 1;
  while (power < n) {
    power <<= 1;
  }
  return power;
}

std::vector<std::complex<double>> ComputeNaiveDFT(const std::vector<double>& samples) {
  const int n = static_cast<int>(samples.size());
  std::vector<std::complex<double>> spectrum(n, {0.0, 0.0});

  for (int k = 0; k < n; ++k) {
    std::complex<double> sum(0.0, 0.0);
    for (int m = 0; m < n; ++m) {
      const double angle = -2.0 * kPi * static_cast<double>(k) * static_cast<double>(m) / static_cast<double>(n);
      sum += samples[m] * std::complex<double>(std::cos(angle), std::sin(angle));
    }
    spectrum[k] = sum;
  }

  return spectrum;
}

std::vector<double> ComputeInverseDFT(const std::vector<std::complex<double>>& spectrum) {
  const int n = static_cast<int>(spectrum.size());
  std::vector<double> reconstructed(n, 0.0);

  for (int m = 0; m < n; ++m) {
    std::complex<double> sum(0.0, 0.0);
    for (int k = 0; k < n; ++k) {
      const double angle = 2.0 * kPi * static_cast<double>(k) * static_cast<double>(m) / static_cast<double>(n);
      sum += spectrum[k] * std::complex<double>(std::cos(angle), std::sin(angle));
    }
    reconstructed[m] = sum.real() / static_cast<double>(n);
  }

  return reconstructed;
}

std::vector<std::complex<double>> ComputeFFT(const std::vector<double>& samples) {
  const int n = static_cast<int>(samples.size());
  std::vector<std::complex<double>> spectrum(n);
  for (int i = 0; i < n; ++i) {
    spectrum[i] = std::complex<double>(samples[i], 0.0);
  }

  int j = 0;
  for (int i = 1; i < n; ++i) {
    int bit = n >> 1;
    while (j & bit) {
      j ^= bit;
      bit >>= 1;
    }
    j ^= bit;
    if (i < j) {
      std::swap(spectrum[i], spectrum[j]);
    }
  }

  for (int len = 2; len <= n; len <<= 1) {
    const double angle = -2.0 * kPi / static_cast<double>(len);
    const std::complex<double> wLen(std::cos(angle), std::sin(angle));
    for (int i = 0; i < n; i += len) {
      std::complex<double> w(1.0, 0.0);
      for (int k = 0; k < len / 2; ++k) {
        const std::complex<double> u = spectrum[i + k];
        const std::complex<double> v = spectrum[i + k + len / 2] * w;
        spectrum[i + k] = u + v;
        spectrum[i + k + len / 2] = u - v;
        w *= wLen;
      }
    }
  }

  return spectrum;
}

std::vector<double> ComputeInverseFFT(const std::vector<std::complex<double>>& spectrum) {
  const int n = static_cast<int>(spectrum.size());
  std::vector<std::complex<double>> values = spectrum;
  std::vector<double> reconstructed(n, 0.0);

  for (std::complex<double>& value : values) {
    value = std::conj(value);
  }

  int j = 0;
  for (int i = 1; i < n; ++i) {
    int bit = n >> 1;
    while (j & bit) {
      j ^= bit;
      bit >>= 1;
    }
    j ^= bit;
    if (i < j) {
      std::swap(values[i], values[j]);
    }
  }

  for (int len = 2; len <= n; len <<= 1) {
    const double angle = -2.0 * kPi / static_cast<double>(len);
    const std::complex<double> wLen(std::cos(angle), std::sin(angle));
    for (int i = 0; i < n; i += len) {
      std::complex<double> w(1.0, 0.0);
      for (int k = 0; k < len / 2; ++k) {
        const std::complex<double> u = values[i + k];
        const std::complex<double> v = values[i + k + len / 2] * w;
        values[i + k] = u + v;
        values[i + k + len / 2] = u - v;
        w *= wLen;
      }
    }
  }

  for (int i = 0; i < n; ++i) {
    reconstructed[i] = std::conj(values[i]).real() / static_cast<double>(n);
  }

  return reconstructed;
}

std::vector<double> BuildAmplitudeSpectrum(const std::vector<std::complex<double>>& spectrum) {
  const int n = static_cast<int>(spectrum.size());
  const int nPositive = n / 2 + 1;
  std::vector<double> amplitudes(nPositive, 0.0);

  for (int k = 0; k < nPositive; ++k) {
    double scale = 2.0 / static_cast<double>(n);
    if (k == 0 || (n % 2 == 0 && k == n / 2)) {
      scale = 1.0 / static_cast<double>(n);
    }
    amplitudes[k] = scale * std::abs(spectrum[k]);
  }

  return amplitudes;
}

TGraph* MakeGraph(const std::vector<double>& x,
                  const std::vector<double>& y,
                  int color,
                  int width,
                  int style = 1,
                  int markerStyle = 0,
                  double markerSize = 0.0) {
  auto* graph = new TGraph(static_cast<int>(x.size()), x.data(), y.data());
  graph->SetLineColor(color);
  graph->SetLineWidth(width);
  graph->SetLineStyle(style);
  if (markerStyle > 0) {
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(markerSize);
    graph->SetMarkerColor(color);
  }
  return graph;
}

double FindMaxAbs(const std::vector<double>& values) {
  double maxAbs = 0.0;
  for (double value : values) {
    maxAbs = std::max(maxAbs, std::fabs(value));
  }
  return maxAbs;
}

double FindMax(const std::vector<double>& values) {
  double maxValue = 0.0;
  for (double value : values) {
    maxValue = std::max(maxValue, value);
  }
  return maxValue;
}

}  // namespace

void sampling_dft_fft_alias_macro() {
  gStyle->SetOptStat(0);

  int nFreq = 0;
  std::cout << "Number of sinusoidal frequencies in the input continuous signal: ";
  std::cin >> nFreq;
  if (!std::cin || nFreq <= 0) {
    std::cerr << "Invalid number of frequencies." << std::endl;
    return;
  }

  std::vector<double> inputFreqs(nFreq, 0.0);
  for (int i = 0; i < nFreq; ++i) {
    std::cout << "Frequency " << i + 1 << " [Hz]: ";
    std::cin >> inputFreqs[i];
    if (!std::cin || inputFreqs[i] < 0.0) {
      std::cerr << "Invalid frequency value." << std::endl;
      return;
    }
  }

  double samplingFrequency = 0.0;
  std::cout << "Sampling frequency fs [Hz]: ";
  std::cin >> samplingFrequency;
  if (!std::cin || samplingFrequency <= 0.0) {
    std::cerr << "Sampling frequency must be positive." << std::endl;
    return;
  }

  double duration = 0.0;
  std::cout << "Total signal duration T [s]: ";
  std::cin >> duration;
  if (!std::cin || duration <= 0.0) {
    std::cerr << "Duration must be positive." << std::endl;
    return;
  }

  int nSamples = static_cast<int>(std::lround(duration * samplingFrequency));
  if (nSamples < 8) {
    std::cerr << "The product T * fs must give at least 8 samples." << std::endl;
    return;
  }
  if (!IsPowerOfTwo(nSamples)) {
    const int requestedSamples = nSamples;
    nSamples = NextPowerOfTwo(nSamples);
    std::cout << "\nRequested N = " << requestedSamples
              << " is not a power of two." << std::endl;
    std::cout << "For the radix-2 FFT, the macro will use N = " << nSamples
              << " samples instead." << std::endl;
  }

  const double dt = 1.0 / samplingFrequency;
  const double nyquistFrequency = samplingFrequency / 2.0;
  const double effectiveDuration = static_cast<double>(nSamples) * dt;

  std::cout << "\n===== Sampling summary =====" << std::endl;
  std::cout << "Sampling frequency fs = " << samplingFrequency << " Hz" << std::endl;
  std::cout << "Nyquist frequency f_N = fs / 2 = " << nyquistFrequency << " Hz" << std::endl;
  std::cout << "Time step dt = " << dt << " s" << std::endl;
  std::cout << "Number of samples N = " << nSamples << std::endl;
  std::cout << "Effective sampled duration = " << effectiveDuration << " s" << std::endl;

  std::cout << "\n===== Input frequencies and aliasing =====" << std::endl;
  for (double f : inputFreqs) {
    const double alias = AliasFrequency(f, samplingFrequency);
    const bool aliased = f > nyquistFrequency;
    std::cout << std::fixed << std::setprecision(6)
              << "Input frequency " << f << " Hz";
    if (aliased) {
      std::cout << "  -> above Nyquist, aliases to " << alias << " Hz";
    } else {
      std::cout << "  -> below Nyquist, observed at " << alias << " Hz";
    }
    std::cout << std::endl;
  }

  std::vector<double> timeContinuous;
  std::vector<double> signalContinuous;
  const int nContinuousPoints = std::max(2000, 20 * nSamples);
  timeContinuous.reserve(nContinuousPoints);
  signalContinuous.reserve(nContinuousPoints);
  for (int i = 0; i < nContinuousPoints; ++i) {
    const double t = effectiveDuration * static_cast<double>(i) / static_cast<double>(nContinuousPoints - 1);
    timeContinuous.push_back(t);
    signalContinuous.push_back(SignalValue(t, inputFreqs));
  }

  std::vector<double> timeSamples(nSamples, 0.0);
  std::vector<double> sampledSignal(nSamples, 0.0);
  for (int n = 0; n < nSamples; ++n) {
    timeSamples[n] = n * dt;
    sampledSignal[n] = SignalValue(timeSamples[n], inputFreqs);
  }

  const auto dftStart = std::chrono::high_resolution_clock::now();
  const std::vector<std::complex<double>> dftSpectrum = ComputeNaiveDFT(sampledSignal);
  const std::vector<double> reconstructedFromDFT = ComputeInverseDFT(dftSpectrum);
  const auto dftStop = std::chrono::high_resolution_clock::now();

  const auto fftStart = std::chrono::high_resolution_clock::now();
  const std::vector<std::complex<double>> fftSpectrum = ComputeFFT(sampledSignal);
  const std::vector<double> reconstructedFromFFT = ComputeInverseFFT(fftSpectrum);
  const auto fftStop = std::chrono::high_resolution_clock::now();

  const std::chrono::duration<double, std::milli> dftTime = dftStop - dftStart;
  const std::chrono::duration<double, std::milli> fftTime = fftStop - fftStart;

  std::cout << "\n===== Timing =====" << std::endl;
  std::cout << "Naive DFT + inverse DFT time: " << dftTime.count() << " ms" << std::endl;
  std::cout << "FFT + inverse FFT time:       " << fftTime.count() << " ms" << std::endl;

  const std::vector<double> dftAmplitude = BuildAmplitudeSpectrum(dftSpectrum);
  const std::vector<double> fftAmplitude = BuildAmplitudeSpectrum(fftSpectrum);

  std::vector<double> freqAxis(dftAmplitude.size(), 0.0);
  for (std::size_t k = 0; k < freqAxis.size(); ++k) {
    freqAxis[k] = static_cast<double>(k) * samplingFrequency / static_cast<double>(nSamples);
  }

  std::vector<double> theoreticalFreqs = inputFreqs;
  std::vector<double> theoreticalAmps(inputFreqs.size(), 1.0);

  std::vector<double> aliasedFreqs(inputFreqs.size(), 0.0);
  std::vector<double> aliasedAmps(inputFreqs.size(), 1.0);
  for (std::size_t i = 0; i < inputFreqs.size(); ++i) {
    aliasedFreqs[i] = AliasFrequency(inputFreqs[i], samplingFrequency);
  }

  auto* canvas = new TCanvas("c_sampling_dft_fft",
                             "Sampling, DFT, FFT, and Aliasing",
                             1500,
                             1200);
  canvas->Divide(2, 3);

  const double timeMax = 1.2 * std::max({FindMaxAbs(signalContinuous),
                                         FindMaxAbs(sampledSignal),
                                         FindMaxAbs(reconstructedFromDFT),
                                         FindMaxAbs(reconstructedFromFFT),
                                         1.0});
  const double freqMax = 1.25 * std::max({FindMax(dftAmplitude), FindMax(fftAmplitude), 1.0});

  canvas->cd(1);
  auto* mgInputTime = new TMultiGraph();
  auto* gContinuous = MakeGraph(timeContinuous, signalContinuous, kBlue + 1, 2);
  auto* gSamples = MakeGraph(timeSamples, sampledSignal, kRed + 1, 2, 1, 20, 0.8);
  mgInputTime->Add(gContinuous, "L");
  mgInputTime->Add(gSamples, "LP");
  mgInputTime->Draw("A");
  mgInputTime->SetTitle("Input signal in time domain;time [s];signal amplitude");
  mgInputTime->GetYaxis()->SetRangeUser(-timeMax, timeMax);
  auto* legInput = new TLegend(0.56, 0.72, 0.88, 0.88);
  legInput->AddEntry(gContinuous, "continuous input", "l");
  legInput->AddEntry(gSamples, "sampled points", "lp");
  legInput->Draw();

  canvas->cd(2);
  auto* frameTheory = gPad->DrawFrame(0.0, 0.0, std::max(*std::max_element(inputFreqs.begin(), inputFreqs.end()), samplingFrequency) * 1.05, 1.35);
  frameTheory->SetTitle("Input and aliased frequencies;frequency [Hz];relative amplitude");
  for (double f : inputFreqs) {
    auto* line = new TLine(f, 0.0, f, 1.0);
    line->SetLineColor(kBlue + 1);
    line->SetLineWidth(2);
    line->Draw();
  }
  for (double fAlias : aliasedFreqs) {
    auto* line = new TLine(fAlias, 0.0, fAlias, 0.85);
    line->SetLineColor(kRed + 1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();
  }
  auto* nyqLine = new TLine(nyquistFrequency, 0.0, nyquistFrequency, 1.2);
  nyqLine->SetLineColor(kGreen + 2);
  nyqLine->SetLineStyle(3);
  nyqLine->SetLineWidth(2);
  nyqLine->Draw();
  auto* gOriginalFreqs = MakeGraph(theoreticalFreqs, theoreticalAmps, kBlue + 1, 2, 1, 20, 1.0);
  auto* gAliasedFreqs = MakeGraph(aliasedFreqs, aliasedAmps, kRed + 1, 2, 1, 24, 1.1);
  gOriginalFreqs->Draw("P");
  gAliasedFreqs->Draw("P");
  auto* legTheory = new TLegend(0.56, 0.68, 0.88, 0.88);
  legTheory->AddEntry(gOriginalFreqs, "input frequencies", "p");
  legTheory->AddEntry(gAliasedFreqs, "observed after sampling", "p");
  legTheory->AddEntry(nyqLine, "Nyquist frequency", "l");
  legTheory->Draw();

  canvas->cd(3);
  auto* mgDFTTime = new TMultiGraph();
  auto* gSampledAgain = MakeGraph(timeSamples, sampledSignal, kGray + 2, 1, 1, 20, 0.75);
  auto* gReconstructedDFT = MakeGraph(timeSamples, reconstructedFromDFT, kMagenta + 1, 2);
  mgDFTTime->Add(gSampledAgain, "P");
  mgDFTTime->Add(gReconstructedDFT, "L");
  mgDFTTime->Draw("A");
  mgDFTTime->SetTitle("Time domain from explicit DFT + inverse DFT;time [s];signal amplitude");
  mgDFTTime->GetYaxis()->SetRangeUser(-timeMax, timeMax);
  auto* legDFTTime = new TLegend(0.49, 0.72, 0.88, 0.88);
  legDFTTime->AddEntry(gSampledAgain, "sampled signal", "p");
  legDFTTime->AddEntry(gReconstructedDFT, "reconstructed from DFT", "l");
  legDFTTime->Draw();

  canvas->cd(4);
  auto* gDFTFreq = MakeGraph(freqAxis, dftAmplitude, kMagenta + 1, 2);
  gDFTFreq->SetTitle("Frequency domain: explicit DFT;frequency [Hz];amplitude");
  gDFTFreq->Draw("AL");
  gDFTFreq->GetXaxis()->SetLimits(0.0, nyquistFrequency);
  gDFTFreq->GetYaxis()->SetRangeUser(0.0, freqMax);

  canvas->cd(5);
  auto* mgFFTTime = new TMultiGraph();
  auto* gSampledFFT = MakeGraph(timeSamples, sampledSignal, kGray + 2, 1, 1, 20, 0.75);
  auto* gReconstructedFFT = MakeGraph(timeSamples, reconstructedFromFFT, kOrange + 7, 2);
  mgFFTTime->Add(gSampledFFT, "P");
  mgFFTTime->Add(gReconstructedFFT, "L");
  mgFFTTime->Draw("A");
  mgFFTTime->SetTitle("Time domain from FFT + inverse FFT;time [s];signal amplitude");
  mgFFTTime->GetYaxis()->SetRangeUser(-timeMax, timeMax);
  auto* legFFTTime = new TLegend(0.49, 0.72, 0.88, 0.88);
  legFFTTime->AddEntry(gSampledFFT, "sampled signal", "p");
  legFFTTime->AddEntry(gReconstructedFFT, "reconstructed from FFT", "l");
  legFFTTime->Draw();

  canvas->cd(6);
  auto* gFFTFreq = MakeGraph(freqAxis, fftAmplitude, kOrange + 7, 2);
  gFFTFreq->SetTitle("Frequency domain: FFT;frequency [Hz];amplitude");
  gFFTFreq->Draw("AL");
  gFFTFreq->GetXaxis()->SetLimits(0.0, nyquistFrequency);
  gFFTFreq->GetYaxis()->SetRangeUser(0.0, freqMax);

  canvas->Update();

  std::cout << "\nCanvas created with 2 columns x 3 rows." << std::endl;
  std::cout << "Left column: time-domain input, DFT reconstruction, FFT reconstruction." << std::endl;
  std::cout << "Right column: input/aliased frequencies, DFT spectrum, FFT spectrum." << std::endl;
}
