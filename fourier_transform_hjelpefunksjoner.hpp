#pragma once

// Dette er en samling hjelpefunksjoner for prosjektoppgaven om foruiertransformasjonen og lyd i MM2

#include <matplot/matplot.h>
#include <AudioFile.h>
#include <fftw3.h>


struct Audio {
    AudioFile<double> file;
    std::vector<std::vector<double>> channels_fourier_transform;

    void load(std::string filename)
    {
        file.load(filename);
    }

    void save(std::string filename)
    {
        file.save(filename);
    }
};


std::vector<std::vector<double>> fourier_transform(AudioFile<double>& audiofile)
{
    std::vector<std::vector<double>> transformed(audiofile.getNumChannels());
    for (int i = 0; i < audiofile.getNumChannels(); ++i) {
        fftw_plan plan;
        std::vector<double>& raw_data = audiofile.samples[i];
        std::vector<double> fft_transform(raw_data.size());

        plan = fftw_plan_r2r_1d(raw_data.size(), raw_data.data(), fft_transform.data(), FFTW_REDFT00, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        transformed[i] = std::move(fft_transform);
    }

    return transformed;
}

void fourier_transform(Audio& audio_struct)
{
    audio_struct.channels_fourier_transform = fourier_transform(audio_struct.file);
}

void add_frequency(AudioFile<double>& audiofile, std::vector<std::vector<double>>& fourier_transformed, double Hz, double amplitude)
{
    uint64_t samplesize = fourier_transformed[0].size();
    uint64_t index = 2 * Hz * samplesize / audiofile.getSampleRate();

    for (int i = 0; i < fourier_transformed.size(); ++i)
        fourier_transformed[i][index] = amplitude;
}

void add_frequency(Audio& audio_struct, double Hz, double amplitude)
{
    add_frequency(audio_struct.file, audio_struct.channels_fourier_transform, Hz, amplitude);
}



void remove_frequency_range(AudioFile<double> audiofile, std::vector<std::vector<double>>& fourier_transformed, double Hz_low, double Hz_high)
{
    uint64_t samplesize = fourier_transformed[0].size();
    uint64_t index_low = 2 * Hz_low * samplesize / audiofile.getSampleRate();
    uint64_t index_high = 2 * Hz_high * samplesize / audiofile.getSampleRate();

    for (int i = 0; i < fourier_transformed.size(); ++i) {
        for (uint64_t j = index_low; j < index_high && j < samplesize; ++j)
            fourier_transformed[i][j] = 0;
    }
}

void remove_frequency_range(Audio& audio_struct, double Hz_low, double Hz_high)
{
    remove_frequency_range(audio_struct.file, audio_struct.channels_fourier_transform, Hz_low, Hz_high);
}

void remove_frequency_over(AudioFile<double>& audiofile, std::vector<std::vector<double>>& fourier_transformed, double Hz)
{
    uint64_t samplesize = fourier_transformed[0].size();
    uint64_t index = 2 * Hz * samplesize / audiofile.getSampleRate();

    for (int i = 0; i < fourier_transformed.size(); ++i) {
        for (uint64_t j = index; j < samplesize; ++j)
            fourier_transformed[i][j] = 0;
    }
}

void remove_frequency_over(Audio& audio_struct, double Hz)
{
    remove_frequency_over(audio_struct.file, audio_struct.channels_fourier_transform, Hz);
}


void remove_frequency_under(AudioFile<double>& audiofile, std::vector<std::vector<double>>& fourier_transformed, double Hz)
{
    uint64_t samplesize = fourier_transformed[0].size();
    uint64_t index = 2 * Hz * samplesize / audiofile.getSampleRate();

    for (int i = 0; i < fourier_transformed.size(); ++i) {
        for (uint64_t j = 0; j < index; ++j)
            fourier_transformed[i][j] = 0;
    }
}

void remove_frequency_under(Audio& audio_struct, double Hz)
{
    remove_frequency_under(audio_struct.file, audio_struct.channels_fourier_transform, Hz);
}

void inverse_fourier_transform(AudioFile<double>& audiofile, std::vector<std::vector<double>>& fourier_transformed)
{
    for (int i = 0; i < audiofile.getNumChannels(); ++i) {
        fftw_plan plan;
        std::vector<double>& fourier_data = fourier_transformed[i];
        std::vector<double>& out_data = audiofile.samples[i];

        plan = fftw_plan_r2r_1d(fourier_data.size(), fourier_data.data(), out_data.data(), FFTW_REDFT00, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
    }

    for (int i = 0; i < audiofile.getNumChannels(); ++i) {
        uint64_t samplesize = audiofile.samples[0].size();
        for (uint64_t j = 0; j < samplesize; ++j) {
            audiofile.samples[i][j] /= samplesize;
        }
    }
}

void inverse_fourier_transform(Audio& audio_struct)
{
    inverse_fourier_transform(audio_struct.file, audio_struct.channels_fourier_transform);
}

void display_frequencies(AudioFile<double>& audiofile, std::vector<std::vector<double>>& fourier_transformed)
{
    double Hz_max = static_cast<double>(audiofile.getSampleRate()) / 2;
    std::vector<double> x = matplot::linspace(0, Hz_max, fourier_transformed[0].size());
    /*
    std::vector<double>display_fft(channel_fourier_transform[0].size() / 100);
    for (uint64_t i = 0; i < channel_fourier_transform[0].size() / 100; ++i) {
        display_fft[i] = channel_fourier_transform[0][i * 100];
    }

    std::vector<double> x = matplot::linspace(0, channel_fourier_transform[0].size() * 100 - 1, channel_fourier_transform[0].size());*/

    matplot::cla();
    matplot::hold(matplot::on);
    for (std::vector<double>& channel : fourier_transformed) {
        matplot::plot(x, channel);
    }
    matplot::show();
}

void display_frequencies(Audio& audio_struct)
{
    display_frequencies(audio_struct.file, audio_struct.channels_fourier_transform);
}

/*
void display_frequencies(AudioFile<double>& audiofile, std::vector<std::vector<double>>& fourier_transformed, double Hz_threshold)
{
    matplot::cla();
    matplot::hold(matplot::on);
    double Hz_max = audiofile.getSampleRate() / 2;
    for (std::vector<double>& channel : fourier_transformed) {
        std::vector<uint64_t> indexes;
        for (uint64_t i = 0; i < channel.size(); ++i) {
            if (abs(channel[i]) == threshold) {
                indexes.push_back(i);
            }
        }//index = 2 * Hz * samplesize / audiofile.getSampleRate()

        std::vector<double> x(indexes.size());
        std::vector<double> displaying_freq(indexes.size());
        for (uint64_t i : indexes) {
            x.push_back(i);
            displaying_freq.push_back(i * audiofile.getSampleRate() / (channel.size() * 2));
        }
        //std::vector<double> x = matplot::transform(indexes, [audiofile, channel](double index) {return index * audiofile.getSampleRate() / (channel.size() * 2); });
        //std::vector<double> displaying_freq = matplot::transform(x, [&audiofile, &channel](double index) {return index * audiofile.getSampleRate() / (channel.size() * 2); });
        matplot::scatter(x, displaying_freq);
    }


    matplot::show();
}*/
