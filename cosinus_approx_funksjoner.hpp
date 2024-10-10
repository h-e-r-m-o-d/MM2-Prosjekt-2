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
};

struct Cos_approkimasjon {
    std::vector<Cosbølge> cosbølger;
    double start;
    double Len;
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

double Cn(const std::vector<double>& verdier, double start, double Len, uint32_t n)
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

Cos_approkimasjon cosinus_approximasjon(const std::vector<double>& punkter, double start, double Len, uint32_t antall_cos)
{
    // "frekvens" betyr multippel av basisfrekvensen PI/Len. Hele cosinusfunksjonen må offsettes med startverdien, fordi de nye cosinusfunksjonene 
    // starter i x = 0

    std::vector<Cosbølge> cosbølger;
    cosbølger.reserve(antall_cos);

    for (int n = 0; n < antall_cos; ++n) {
        Cosbølge temp;
        temp.amplitude = Cn(punkter, start, Len, n);
        temp.frekvens = n;
        cosbølger.push_back(temp);
    }

    Cos_approkimasjon approksimasjon;
    approksimasjon.cosbølger = std::move(cosbølger);
    approksimasjon.Len = Len;
    approksimasjon.start = start;

    return approksimasjon;
}

std::vector<double> parse_approximation(const Cos_approkimasjon& approksimasjon, const std::vector<double>& x)
{
    std::vector<double> vals(x.size());

    for (const Cosbølge& i : approksimasjon.cosbølger) {

        double amplitude = i.amplitude;
        uint32_t frekvens = i.frekvens;

        for (int j = 0; j < x.size(); ++j) {
            vals[j] += amplitude * cos(PI * (x[j] - approksimasjon.start) * frekvens / approksimasjon.Len);
        }
    }

    return vals;
}

void plot_decomposed(const matplot::axes_handle plotn, const Cos_approkimasjon& approksimasjon, const std::vector<double>& x, uint32_t n, double threshold = 0.0)
{
    //plotter hver enkelt multippel av basisfrekvensen for seg selv.

    std::vector<double> vals(x.size());

    for (int i = 0; i < approksimasjon.cosbølger.size() && i < n; ++i) {

        double amplitude = approksimasjon.cosbølger[i].amplitude;
        uint32_t frekvens = approksimasjon.cosbølger[i].frekvens;

        if (abs(amplitude) < threshold) {
            continue;
        }

        for (int j = 0; j < x.size(); ++j) {
            vals[j] = amplitude * cos(PI * (x[j] - approksimasjon.start) * frekvens / approksimasjon.Len);
        }

        if (plotn != nullptr) {
            matplot::plot(plotn, x, vals)->line_width(1).display_name(std::to_string(frekvens));
        }
        else {
            matplot::plot(x, vals)->line_width(1).display_name(std::to_string(frekvens));
        }
    }
}
void plot_decomposed(const Cos_approkimasjon& approksimasjon, const std::vector<double>& x, uint32_t n = 10, double threshold = 0.0)
{
    plot_decomposed(nullptr, approksimasjon, x, n, threshold);
}

std::vector<double> parse_varmeligning_t(const Cos_approkimasjon& approksimasjon, const std::vector<double>& x, double t, double heat_conductivity)
{
    // hvordan varmefordelingen blir over "staven" etter "t" tid.
    std::vector<double> vals(x.size());

    for (const Cosbølge& i : approksimasjon.cosbølger) {

        double amplitude = i.amplitude;
        uint32_t frekvens = i.frekvens;

        double omega = PI * frekvens / approksimasjon.Len;

        for (int j = 0; j < x.size(); ++j) {
            vals[j] += amplitude * cos(omega * (x[j] - approksimasjon.start)) * std::exp(-heat_conductivity * omega * omega * t);
        }
    }

    return vals;
}

double parse_varmeligning_x_t(const Cos_approkimasjon& approksimasjon, double x, double t, double heat_conductivity)
{
    // hva varmen blir i et punkt etter t.
    double varme_punkt = 0;
    for (const Cosbølge& i : approksimasjon.cosbølger) {

        double amplitude = i.amplitude;
        uint32_t frekvens = i.frekvens;

        double omega = PI * frekvens / approksimasjon.Len;

        varme_punkt += amplitude * cos(omega * (x - approksimasjon.start)) * std::exp(-heat_conductivity * omega * omega * t);
    }
    return varme_punkt;
}

void cosinus_approximasjon_range(std::vector<Cosbølge>* approksimasjon_temp, const std::vector<double>* punkter, double start, double Len, uint32_t cos_start, uint32_t cos_slutt, char* finished)
{
    // denne er brukt til multithreading. se "cosinus approximasjon" for å lette forstå hva funksjonene gjør

    approksimasjon_temp->clear();
    approksimasjon_temp->reserve(cos_slutt - cos_start);

    for (int n = cos_start; n < cos_slutt; ++n) {
        Cosbølge temp;
        temp.amplitude = Cn(*punkter, start, Len, n);
        temp.frekvens = n;
        approksimasjon_temp->push_back(temp);
    }
    *finished = true;
}


Cos_approkimasjon cosinus_approximasjon_threaded(const std::vector<double>& punkter, double start, double Len, uint32_t antall_cos)
{
    // dette er en mutithreaded versjon av funksjonen "cosinus_approximasjon". Se den for å lettere forstå hva funksjonene gjør

    if (antall_cos <= 1000) {
        return cosinus_approximasjon(punkter, start, Len, antall_cos);
    }

    std::vector<Cosbølge> cosbølger;
    cosbølger.reserve(antall_cos);

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
                        cosbølger.push_back(thread_vector_return[i][j]);
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
                cosbølger.push_back(thread_vector_return[i][j]);
            }
        }
        
    }
    
    Cos_approkimasjon approksimasjon;
    approksimasjon.cosbølger = std::move(cosbølger);
    approksimasjon.Len = Len;
    approksimasjon.start = start;

    return approksimasjon;
}