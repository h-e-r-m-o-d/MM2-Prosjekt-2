
#include <iostream>

#include <matplot/matplot.h>
#include <fftw3.h>

#include <AudioFile.h>

#include "fourier_transform_hjelpefunksjoner.hpp"


int main()
{
    Audio lyd;
    lyd.load("fourier.wav");
    //lyd.file.printSummary();
    fourier_transform(lyd);

    //remove_frequency_over(audiofile, channel_fourier_transform, 12000);
    //add_frequency(lyd, 12000, 500);
    //add_frequency(lyd, 12100, 500);
    //add_frequency(lyd, 12200, 500);
    //add_frequency(lyd, 500, 500);
    //add_frequency(lyd, 510, 500);
    //add_frequency(lyd, 520, 500);

    display_frequencies_cutoff(lyd, 15000);
    
    //remove_frequency_range(lyd, 150, 270);
    remove_frequency_range(lyd, 6600,7100);

    display_frequencies_cutoff(lyd,15000);

    inverse_fourier_transform(lyd);
    

    lyd.save("lyd_edit.wav");

    //-------------------------------------

    //lyd.load("lyd_edit.wav");

    //fourier_transform(lyd);

    //remove_frequency_over(audiofile, channel_fourier_transform, 12000);
    //remove_frequency_over(lyd, 11500);

    //display_frequencies(audiofile, channel_fourier_transform);

    //inverse_fourier_transform(lyd);

    //lyd.save("lyd_fix.wav");

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

