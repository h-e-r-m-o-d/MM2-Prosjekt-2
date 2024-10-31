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

        assert(raw_data.size() <= INT_MAX);

        plan = fftw_plan_r2r_1d(static_cast<int>(raw_data.size()), raw_data.data(), fft_transform.data(), FFTW_REDFT00, FFTW_ESTIMATE);
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
    uint64_t index = static_cast<uint64_t>(2 * Hz * samplesize / audiofile.getSampleRate());

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
    uint64_t index_low = static_cast<uint64_t>(2 * Hz_low * samplesize / audiofile.getSampleRate());
    uint64_t index_high = static_cast<uint64_t>(2 * Hz_high * samplesize / audiofile.getSampleRate());

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
    uint64_t index = static_cast<uint64_t>(2 * Hz * samplesize / audiofile.getSampleRate());

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
    uint64_t index = static_cast<uint64_t>(2 * Hz * samplesize / audiofile.getSampleRate());

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

        assert(fourier_data.size() <= INT_MAX);

        plan = fftw_plan_r2r_1d(static_cast<int>(fourier_data.size()), fourier_data.data(), out_data.data(), FFTW_REDFT00, FFTW_ESTIMATE);
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
    /*for (std::vector<double>& channel : fourier_transformed) {
        matplot::plot(x, channel);
    }*/
    matplot::xlabel("Hz");

    for (size_t i = 0; i < fourier_transformed.size(); ++i) {
        matplot::plot(x, fourier_transformed[i])->display_name("Channel " + std::to_string(i));
    }
        
    matplot::show();
}

void display_frequencies(Audio& audio_struct)
{
    display_frequencies(audio_struct.file, audio_struct.channels_fourier_transform);
}

void display_sound(AudioFile<double>& audiofile)
{
    double time = static_cast<double>(audiofile.getLengthInSeconds()) / 2;
    std::vector<double> x = matplot::linspace(0, time, audiofile.samples[0].size());
    /*
    std::vector<double>display_fft(channel_fourier_transform[0].size() / 100);
    for (uint64_t i = 0; i < channel_fourier_transform[0].size() / 100; ++i) {
        display_fft[i] = channel_fourier_transform[0][i * 100];
    }

    std::vector<double> x = matplot::linspace(0, channel_fourier_transform[0].size() * 100 - 1, channel_fourier_transform[0].size());*/

    matplot::cla();
    matplot::hold(matplot::on);
    for (std::vector<double>& channel : audiofile.samples) {
        matplot::plot(x, channel);
    }

    matplot::show();
}

void display_sound(Audio& audio_struct)
{
    display_sound(audio_struct.file);
}


void display_frequencies(AudioFile<double>& audiofile, std::vector<std::vector<double>>& fourier_transformed, double percent_threshold)
{
    double highest = 0;

    for (std::vector<double>& channel : fourier_transformed) {
        for (double i : channel) {
            if (abs(i) > highest) {
                highest = abs(i);
            }
        }
    }

    matplot::cla();
    matplot::hold(matplot::on);
    matplot::xlabel("Hz");
    double Hz_max = static_cast<double>(audiofile.getSampleRate()) / 2;

    int i = 0;
    for (std::vector<double>& channel : fourier_transformed) {
        std::vector<uint64_t> indexes;
        for (uint64_t i = 0; i < channel.size(); ++i) {
            if (abs(channel[i]) >= highest*percent_threshold/100) {
                indexes.push_back(i);
            }
        }//index = 2 * Hz * samplesize / audiofile.getSampleRate()

        std::vector<double> x(indexes.size());
        std::vector<double> displaying_freq(indexes.size());

        uint64_t last_i = 0;
        for (uint64_t i : indexes) {
            if (i - last_i > 1) {
                x.push_back(static_cast<double>((last_i + 1) * audiofile.getSampleRate() / (channel.size() * 2)));
                displaying_freq.push_back(0);
                x.push_back(static_cast<double>((i - 1) * audiofile.getSampleRate() / (channel.size() * 2)));
                displaying_freq.push_back(0);
            }
            x.push_back(static_cast<double>(i * audiofile.getSampleRate() / (channel.size() * 2)));
            displaying_freq.push_back(channel[i]);
            //i * audiofile.getSampleRate() / (channel.size() * 2)
            last_i = i;
        }

        if (last_i != indexes[indexes.size() - 1]) {
            x.push_back(static_cast<double>((last_i + 1) * audiofile.getSampleRate() / (channel.size() * 2)));
            displaying_freq.push_back(0);
        }

        //std::vector<double> x = matplot::transform(indexes, [audiofile, channel](double index) {return index * audiofile.getSampleRate() / (channel.size() * 2); });
        //std::vector<double> displaying_freq = matplot::transform(x, [&audiofile, &channel](double index) {return index * audiofile.getSampleRate() / (channel.size() * 2); });
        matplot::plot(x, displaying_freq)->display_name("Channel " + std::to_string(i));
        ++i;
    }


    matplot::show();
}

void display_frequencies(Audio& audio_struct, double percent_threshold)
{
    display_frequencies(audio_struct.file, audio_struct.channels_fourier_transform, percent_threshold);
}

void display_frequencies_cutoff(AudioFile<double>& audiofile, std::vector<std::vector<double>>& fourier_transformed, double Hz_cutoff)
{
    double Hz_max = static_cast<double>(audiofile.getSampleRate()) / 2;
    uint64_t samplesize = audiofile.samples[0].size();
    uint64_t index_max = 2 * static_cast<uint64_t>(Hz_cutoff * samplesize / audiofile.getSampleRate());

    
    std::vector<double> x = matplot::linspace(0, Hz_cutoff, index_max+1);
    
    /*
    std::vector<double>display_fft(channel_fourier_transform[0].size() / 100);
    for (uint64_t i = 0; i < channel_fourier_transform[0].size() / 100; ++i) {
        display_fft[i] = channel_fourier_transform[0][i * 100];
    }

    std::vector<double> x = matplot::linspace(0, channel_fourier_transform[0].size() * 100 - 1, channel_fourier_transform[0].size());*/

    matplot::cla();
    matplot::hold(matplot::on);
    matplot::xlabel("Hz");
    int i = 0;
    for (std::vector<double>& channel : fourier_transformed) {
        std::vector<double> display_val(channel.begin(), channel.begin() + index_max);
        matplot::plot(x, display_val)->display_name("Channel " + std::to_string(i));
        ++i;
    }
    matplot::legend();
    matplot::show();
}

void display_frequencies_cutoff(Audio& audio_struct, double Hz_cutoff)
{
    display_frequencies_cutoff(audio_struct.file, audio_struct.channels_fourier_transform, Hz_cutoff);
}




// ingen av disse komprimerer egentlig :)
void komprimere(AudioFile<double>& audiofile, std::vector<std::vector<double>>& fourier_transformed, double percent)
{
    const uint64_t samplesize = audiofile.samples[0].size();
    const uint64_t start_index = percent / 100 * samplesize;

    for (std::vector<double>& channel : fourier_transformed) {
        for (uint64_t i = start_index; i < channel.size() - 1; i += 2) {
            channel[i] += channel[i + 1];
            channel[i] /= 2;
            channel[i + 1] = 0;
        }
    }
}
void komprimere(Audio& audio_struct, double percent)
{
    komprimere(audio_struct.file ,audio_struct.channels_fourier_transform, percent);
}


//denne bare setter ned samplingraten
void komprimere2(AudioFile<double>& audiofile, std::vector<std::vector<double>>& fourier_transformed, double percent)
{
    const uint64_t samplesize = audiofile.samples[0].size();
    const uint64_t end_index = percent / 100 * samplesize;

    std::vector<std::vector<double>> new_vec;
    new_vec.reserve(fourier_transformed.size());

    for (std::vector<double>& channel : fourier_transformed) {
        std::vector<double> temp;
        temp.reserve(end_index);
        for (uint64_t i = 0; i < end_index; i += 1) {
            temp.push_back(channel[i]);
            //channel[i] += channel[i + 1];
            //channel[i] /= 2;
            //channel[i] = 0;
            //channel[i + 1] = 0;
        }
        new_vec.push_back(temp);
        channel.resize(end_index);
    }
    
    audiofile.samples = std::move(new_vec);
    audiofile.setSampleRate(audiofile.getSampleRate() * percent / 100);
}

void komprimere2(Audio& audio_struct, double percent)
{
    komprimere2(audio_struct.file, audio_struct.channels_fourier_transform, percent);
}