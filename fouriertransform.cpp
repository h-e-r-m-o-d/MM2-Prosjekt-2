
#include <iostream>

#include <matplot/matplot.h>
#include <fftw3.h>

#include <AudioFile.h>
#include <fstream>

constexpr double PI = 3.14159265359;
namespace plt = matplot;
 
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

void add_frequency(AudioFile<double> audiofile, std::vector<std::vector<double>>& fourier_transformed, double Hz, double amplitude)
{
    uint64_t samplesize = fourier_transformed[0].size();
    uint64_t index = 2 * Hz * samplesize / audiofile.getSampleRate();

    for (int i = 0; i < fourier_transformed.size(); ++i)
        fourier_transformed[i][index] = amplitude;
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

void remove_frequency_over(AudioFile<double>& audiofile, std::vector<std::vector<double>>& fourier_transformed, double Hz)
{
    uint64_t samplesize = fourier_transformed[0].size();
    uint64_t index = 2* Hz * samplesize / audiofile.getSampleRate();

    for (int i = 0; i < fourier_transformed.size(); ++i) {
        for (uint64_t j = index; j < samplesize; ++j)
            fourier_transformed[i][j] = 0;
    }
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

void display_frequencies(AudioFile<double>& audiofile, std::vector<std::vector<double>>& fourier_transformed)
{
    double Hz_max = audiofile.getSampleRate() / 2;
    std::vector<double> x = plt::linspace(0, Hz_max, fourier_transformed[0].size());
    /*
    std::vector<double>display_fft(channel_fourier_transform[0].size() / 100);
    for (uint64_t i = 0; i < channel_fourier_transform[0].size() / 100; ++i) {
        display_fft[i] = channel_fourier_transform[0][i * 100];
    }

    std::vector<double> x = plt::linspace(0, channel_fourier_transform[0].size() * 100 - 1, channel_fourier_transform[0].size());*/

    plt::cla();
    plt::hold(plt::on);
    for (std::vector<double>& channel : fourier_transformed) {
        plt::area(x, channel);
    }
    plt::show();
}

void display_frequencies(AudioFile<double>& audiofile, std::vector<std::vector<double>>& fourier_transformed, double threshold)
{
    plt::cla();
    plt::hold(plt::on);
    double Hz_max = audiofile.getSampleRate() / 2;
    for (std::vector<double>& channel : fourier_transformed) {
        std::vector<uint64_t> indexes;
        for (uint64_t i = 0; i < channel.size(); ++i) {
            if (abs(channel[i]) == threshold) {
                indexes.push_back(i);
            }
        }//index = 2 * Hz * samplesize / audiofile.getSampleRate()
            std::vector<double> x = plt::transform(x, [audiofile, channel](double index) {return index * audiofile.getSampleRate() / (channel.size() * 2); });
            std::vector<double> displaying_freq = plt::transform(x, [&audiofile, &channel](double index) {return index*audiofile.getSampleRate() / (channel.size() * 2); });
            plt::scatter(x, displaying_freq);
    }


    plt::show();
}



int main()
{
    AudioFile<double> audiofile;
    audiofile.load("cell.wav");
    std::cout << audiofile.getLengthInSeconds() << '\n';
    std::vector<std::vector<double>> channel_fourier_transform = fourier_transform(audiofile);

    //remove_frequency_over(audiofile, channel_fourier_transform, 12000);
    add_frequency(audiofile, channel_fourier_transform, 12000, 500);
    add_frequency(audiofile, channel_fourier_transform, 12100, 500);
    add_frequency(audiofile, channel_fourier_transform, 12200, 500);
    add_frequency(audiofile, channel_fourier_transform, 500, 500);
    add_frequency(audiofile, channel_fourier_transform, 510, 500);
    add_frequency(audiofile, channel_fourier_transform, 520, 500);

    display_frequencies(audiofile, channel_fourier_transform,200);

    inverse_fourier_transform(audiofile, channel_fourier_transform);

    audiofile.save("cell_edit.wav");


    audiofile.load("cell_edit.wav");
    std::cout << audiofile.getLengthInSeconds() << '\n';
    channel_fourier_transform = fourier_transform(audiofile);

    //remove_frequency_over(audiofile, channel_fourier_transform, 12000);
    remove_frequency_over(audiofile, channel_fourier_transform, 11500);

    //display_frequencies(audiofile, channel_fourier_transform);

    inverse_fourier_transform(audiofile, channel_fourier_transform);

    audiofile.save("cell_fix.wav");

    //std::vector<double> audio = audiofile.samples[0];
    /*
    std::vector<double>& raw_data = audiofile.samples[0];

    std::vector<double> fft_transform(raw_data.size());
    fftw_plan plan;

    plan = fftw_plan_r2r_1d(raw_data.size(), raw_data.data(), fft_transform.data(), FFTW_REDFT00, FFTW_ESTIMATE);
    fftw_execute(plan); // signal to spectrum
    fftw_destroy_plan(plan);

    //mess up audio

    fft_transform[50000] = 10000;

    std::vector<double> edit_raw(raw_data.size());
    plan = fftw_plan_r2r_1d(raw_data.size(), fft_transform.data(), edit_raw.data(), FFTW_REDFT00, FFTW_ESTIMATE);
    fftw_execute(plan); // signal to spectrum
    fftw_destroy_plan(plan);

    fftw_cleanup();

    for (uint64_t i = 0; i < fft_transform.size(); ++i) {
        edit_raw[i] /= fft_transform.size() * 2;
        //display_fft[i] = raw_data[i*100];
    }
   
    
    std::vector<double>display_fft(fft_transform.size() / 100);
    */
    /*std::vector<double>display_fft(audio.size() / 100);
    for (uint64_t i = 0; i < audio.size()/100; ++i) {
        display_fft[i] = audio[i*100];
    }

    std::vector<double> x = plt::linspace(0, (audio.size()-1), audio.size());

    //plt::plot(x, display_fft); *
    plt::area(x, display_fft);
    plt::show();

    //std::vector<double> x = plt::linspace(0, (fft_transform.size() - 1), fft_transform.size());
    //plt::area(x, display_fft);
    //plt::show();*/
    
    return 0;
}

