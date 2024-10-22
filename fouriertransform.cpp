
#include <iostream>

#include <matplot/matplot.h>
#include <fftw3.h>

#include <AudioFile.h>
#include <fstream>

constexpr double PI = 3.14159265359;

double f(double x) // func 0-8
{
    return cos(2 * PI * 3 * x) * exp(-PI * (x - 4) * (x - 4));
}


namespace plt = matplot;
int main()
{
    AudioFile<double> audiofile;
    //audiofile.load("cell.wav");
    std::vector<double> x_ = plt::linspace(0, 8, 1000);
    std::vector<double> raw_data = plt::transform(x_, [](double x) {return f(x); }); //audiofile.samples[0];

    std::vector<double> fft_transform(raw_data.size());
    fftw_plan plan;

    plan = fftw_plan_r2r_1d(raw_data.size(), raw_data.data(), fft_transform.data(), FFTW_REDFT00, FFTW_ESTIMATE);
    fftw_execute(plan); // signal to spectrum
    fftw_destroy_plan(plan);

    fftw_cleanup();
    
    /*std::vector<double>display_fft(fft_transform.size() / 100);

    for (int i = 0; i < fft_transform.size()/100; ++i) {
        display_fft[i] = abs(fft_transform[i*100]);
    }
    
    std::vector<double> x = plt::linspace(0, (fft_transform.size()-1), display_fft.size());

    //plt::plot(x, display_fft);*/
    plt::area(x_, raw_data);
    plt::show();

    std::vector<double> x = plt::linspace(0, (fft_transform.size() - 1), fft_transform.size());
    plt::area(x, fft_transform);
    plt::show();

    return 0;
}

