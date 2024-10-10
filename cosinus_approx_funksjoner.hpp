#pragma once

// Dette er en samling hjelpefunksjoner for prosjektoppgaven om varmelikningen i MM2

#include <matplot/matplot.h>
#include <math.h>

// for multithreading
#include <chrono>
#include <thread>

struct Cosbølge {
    double amplitude;
    uint32_t frekvens;
    double offset;
};
constexpr double PI = 3.14159265358979323846;

double integral(double func(double), double start, double slutt, double delta)
{
    double sum = 0;
    for (double i = start; i < slutt; i += delta) {
        sum += func(i);
    }
    return sum;
}

double Cn(const std::vector<double>& verdier, double start, double Len, int n)
{
    // antar at verdiene er jevnt fordelt over et tidsrom L, fra start til slutt
    // regner ut 2/Len * integral(verdier*cos(PI*n*x)*dx)
    // ledd 0 (frekvens 0) er kun gjennomsnitt av funksjonen

    double sum = 0;
    double x = 0;
    double dx = Len / verdier.size();

    for (double val : verdier) {
        sum += val * cos(PI * n * x / Len) * dx;
        x += dx;
    }
    sum /= Len;

    if (n != 0) { // C_0 / 2
        sum *= 2;
    }


    return sum;
}

std::vector<Cosbølge> cosinus_approximasjon(const std::vector<double>& punkter, double start, double Len, int antall_cos)
{
    // "frekvens" betyr multippel av basisfrekvensen PI/Len. Hele cosinusfunksjonen må offsettes med startverdien, fordi de nye cosinusfunksjonene 
    // starter i x = 0
    std::vector<Cosbølge> approksimasjon;
    approksimasjon.reserve(antall_cos);

    for (int n = 0; n < antall_cos; ++n) {
        Cosbølge temp;
        temp.amplitude = Cn(punkter, start, Len, n);
        temp.frekvens = n;
        temp.offset = start;
        approksimasjon.push_back(temp);
    }

    return approksimasjon;
}

std::vector<double> parse_approximation(const std::vector<Cosbølge>& approksimasjon, const std::vector<double>& x)
{
    std::vector<double> vals(x.size());
    double Len = x.back() - x.front();

    for (const Cosbølge& i : approksimasjon) {

        double amplitude = i.amplitude;
        uint32_t frekvens = i.frekvens;
        double offset = i.offset;

        for (int j = 0; j < x.size(); ++j) {
            vals[j] += amplitude * cos(PI * (x[j] - offset) * frekvens / Len);
        }
    }

    return vals;
}

void plot_decomposed(const matplot::axes_handle plotn, const std::vector<Cosbølge>& approksimasjon, const std::vector<double>& x, int n = 10, double threshold = 0.0)
{
    //plotter hver enkelt multippel av basisfrekvensen for seg selv.

    std::vector<double> vals(x.size());
    double Len = x.back() - x.front();

    for (int i = 0; i < approksimasjon.size() && i < n; ++i) {

        double amplitude = approksimasjon[i].amplitude;
        uint32_t frekvens = approksimasjon[i].frekvens;
        double offset = approksimasjon[i].offset;

        if (abs(amplitude) < threshold) {
            continue;
        }

        for (int j = 0; j < x.size(); ++j) {
            vals[j] = amplitude * cos(PI * (x[j] - offset) * frekvens / Len);
        }
        if (plotn != 0) {
            matplot::plot(plotn, x, vals)->line_width(1).display_name(std::to_string(frekvens));
        }
        else {
            matplot::plot(x, vals)->line_width(1).display_name(std::to_string(frekvens));
        }
    }
}
void plot_decomposed(const std::vector<Cosbølge>& approksimasjon, const std::vector<double>& x, int n = 10, double threshold = 0.0)
{
    plot_decomposed(0, approksimasjon, x, n, threshold);
}

std::vector<double> parse_varmeligning_t(const std::vector<Cosbølge>& approksimasjon, const std::vector<double>& x, double t, double heat_conductivity)
{
    // hvordan varmefordelingen blir over "staven" etter "t" tid.
    std::vector<double> vals(x.size());
    double Len = x.back() - x.front();

    for (const Cosbølge& i : approksimasjon) {

        double amplitude = i.amplitude;
        uint32_t frekvens = i.frekvens;
        double offset = i.offset;

        double omega = PI * frekvens / Len;

        for (int j = 0; j < x.size(); ++j) {
            vals[j] += amplitude * cos(omega * (x[j] - offset)) * std::exp(-heat_conductivity * omega * omega * t);
        }
    }

    return vals;
}

double parse_varmeligning_x_t(const std::vector<Cosbølge>& approksimasjon, double Len, double x, double t, double heat_conductivity)
{
    // hva varmen blir i et punkt etter t.
    double varme_punkt = 0;
    for (const Cosbølge& i : approksimasjon) {

        double amplitude = i.amplitude;
        uint32_t frekvens = i.frekvens;
        double offset = i.offset;

        double omega = PI * frekvens / Len;

        varme_punkt += amplitude * cos(omega * (x - offset)) * std::exp(-heat_conductivity * omega * omega * t);
    }
    return varme_punkt;
}

void cosinus_approximasjon_range(std::vector<Cosbølge>* approksimasjon_temp, const std::vector<double>* punkter, double start, double Len, int cos_start, int cos_slutt, char* finished)
{
    approksimasjon_temp->clear();
    approksimasjon_temp->reserve(cos_slutt - cos_start);

    for (int n = cos_start; n < cos_slutt; ++n) {
        Cosbølge temp;
        temp.amplitude = Cn(*punkter, start, Len, n);
        temp.frekvens = n;
        temp.offset = start;
        approksimasjon_temp->push_back(temp);
    }
    *finished = true;
}


std::vector<Cosbølge> cosinus_approximasjon_threaded(const std::vector<double>& punkter, double start, double Len, int antall_cos)
{
    // dette er en mutithreaded versjon av funksjonen "cosinus_approximasjon". Se den for å lettere forstå hva funksjonene gjør

    if (antall_cos <= 1000) {
        return cosinus_approximasjon(punkter, start, Len, antall_cos);
    }

    std::vector<Cosbølge> approksimasjon;
    approksimasjon.reserve(antall_cos);

    int current_cos = 0;

    int sz = std::thread::hardware_concurrency();
    std::vector<std::thread> threadpool(sz);
    std::vector<std::vector<Cosbølge>> thread_vector_return(sz);
    std::vector<char> finished(sz);

    for (int i = 0; i < sz && current_cos < antall_cos; ++i) {
        if (current_cos + 1000 < antall_cos) {
            threadpool[i] = std::thread{ cosinus_approximasjon_range, &(thread_vector_return[i]), &punkter, start, Len, current_cos,  current_cos + 1000, &(finished[i]) };
            current_cos += 1000;
        }
        else {
            threadpool[i] = std::thread{ cosinus_approximasjon_range, &(thread_vector_return[i]), &punkter, start, Len, current_cos,  antall_cos, &(finished[i]) };
            current_cos = antall_cos;
        }
    }

    while (current_cos < antall_cos) {

        for (int i = 0; i < sz; ++i) {
            if (finished[i]) {
                if (threadpool[i].joinable()) {
                    finished[i] = false;
                    threadpool[i].join();
                    for (int j = 0; j < thread_vector_return[i].size(); ++j) {
                        approksimasjon.push_back(thread_vector_return[i][j]);
                    }
                    if (current_cos + 1000 < antall_cos) {
                        threadpool[i] = std::thread{ cosinus_approximasjon_range, &(thread_vector_return[i]), &punkter, start, Len, current_cos,  current_cos + 1000, &(finished[i]) };
                        current_cos += 1000;
                    }
                    else {
                        threadpool[i] = std::thread{ cosinus_approximasjon_range, &(thread_vector_return[i]), &punkter, start, Len, current_cos,  antall_cos, &(finished[i]) };
                        current_cos = antall_cos;
                    }
                }
            }
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }

    for (int i = 0; i < sz; ++i) {
        if (threadpool[i].joinable()) {
            threadpool[i].join();
            for (int j = 0; j < thread_vector_return[i].size(); ++j) {
                approksimasjon.push_back(thread_vector_return[i][j]);
            }
        }
        
    }
    
    return approksimasjon;
}