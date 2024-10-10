
#include <iostream>
#include "cosinus_approx_funksjoner.hpp"
#include <matplot/matplot.h>

void figur1();
void figur2();
void figur3();
void figur4();
void figur5();
void figur6og7();

namespace plt = matplot;

int main()
{
    /*figur1();
    plt::cla();
    figur2();
    plt::cla();
    figur3();
    plt::cla();
    figur4();
    plt::cla();
    figur5();
    plt::cla();*/
    figur6og7();

    return 0;
}

double f1(double x)
{
    return x * x;
}
double f2(double x)
{
    return 4;
}
double f3(double x)
{
    return 2 * cos(x) + 5 * sin(4 * x) + 5 * abs(x) / (abs(x) + 1);
}
double f5(double x)
{
    return sin(2 * PI * x);
}
double f6(double x)
{
    if (x < 0)
        return -1;
    if (x > 0)
        return 1;
    return 0;
}


void figur1() {

    double x_start = 0;
    double x_slutt = 5;
    size_t x_num = 1000;

    double Len = x_slutt - x_start;

    std::vector<Cosbølge> cosbølger;
    std::vector<double> x_utregning = plt::linspace(x_start, x_slutt, x_num);
    std::vector<double> y_utregning = plt::transform(x_utregning, [](double x) {return f1(x); });

    Cos_approkimasjon approksimajon = cosinus_approximasjon(y_utregning, x_start, Len, 10);


    std::vector<double> x = plt::linspace(x_start, x_slutt, 500);
    std::vector<double> y = plt::transform(x, [](double x) {return f1(x); });
    std::vector<double> y_approx = parse_approximation(approksimajon, x);


    /* plt::hold(plt::on);
     plt::title("Graf mot tilnærming");
     plt::plot(x, y)->display_name("f(x) = x*x");
     plt::plot(x, y_approx)->display_name("Cosinustilnærming");
     plt::legend();
     plt::show();

     plt::cla();*/

    plt::hold(plt::on);
    plt::title("Dekomponering av tilnærming");
    plot_decomposed(approksimajon, x, 10);
    plt::legend()->title("Antall halve perioder");
    plt::show();

    /*plt::cla();

    plt::title("Amplituder i tilnærming");
    std::vector<int> freq;
    std::vector<double> freq_amplitudes;
    freq.reserve(approksimasjon.size());
    freq_amplitudes.reserve(approksimasjon.size());

    for (int i = 0; i < approksimasjon.size() && i < 10; ++i) {
        Cosbølge& j = approksimasjon[i];
        freq.push_back(j.frekvens);
        freq_amplitudes.push_back(j.amplitude);
    }

    plt::bar(freq, freq_amplitudes);
    plt::xlabel("n");
    plt::ylabel("Amplitude");
    plt::legend();
    plt::show();*/
}

void figur2()
{
    double x_start = 0;
    double x_slutt = 5;
    size_t x_num = 100;

    double Len = x_slutt - x_start;

    Cos_approkimasjon approksimasjon;
    std::vector<double> x_utregning = plt::linspace(x_start, x_slutt, x_num);
    std::vector<double> y_utregning = plt::transform(x_utregning, [](double x) {return f2(x); });

    approksimasjon = cosinus_approximasjon(y_utregning, x_start, Len, 10);

    std::vector<double> x = plt::linspace(x_start, x_slutt, 1000);
    std::vector<double> y = plt::transform(x, [](double x) {return f2(x); });
    std::vector<double> y_approx = parse_approximation(approksimasjon, x);


    plt::hold(plt::on);
    plt::title("Graf mot tilnærming");
    plt::plot(x, y)->display_name("f(x) = 4");
    plt::plot(x, y_approx)->display_name("Cosinustilnærming");
    plt::legend();
    plt::show();
}

void figur3()
{
    double x_start = -10;
    double x_slutt = 10;
    size_t x_num = 100000;

    double Len = x_slutt - x_start;

    Cos_approkimasjon approksimasjon;
    std::vector<double> x_utregning = plt::linspace(x_start, x_slutt, x_num);
    std::vector<double> y_utregning = plt::transform(x_utregning, [](double x) {return f3(x); });



    approksimasjon = cosinus_approximasjon(y_utregning, x_start, Len, 25);

    std::vector<double> x = plt::linspace(x_start, x_slutt, 1000);
    std::vector<double> y = plt::transform(x, [](double x) {return f3(x); });
    std::vector<double> y_approx = parse_approximation(approksimasjon, x);

    plt::hold(plt::on);
    plt::title("Graf mot tilnærming");
    plt::plot(x, y)->display_name("f(x) = 2*cos(x) + 5*sin(4*x) + 5*abs(x)/(abs(x)+1)");
    plt::plot(x, y_approx)->display_name("Cosinustilnærming");
    plt::legend();
    plt::show();
}

void figur4()
{
    double x_start = -10;
    double x_slutt = 10;
    size_t x_num = 100000;

    double Len = x_slutt - x_start;

    Cos_approkimasjon approksimasjon;
    std::vector<double> x_utregning = plt::linspace(x_start, x_slutt, x_num);
    std::vector<double> y_utregning = plt::transform(x_utregning, [](double x) {return f3(x); });


    approksimasjon = cosinus_approximasjon(y_utregning, x_start, Len, 60);

    plt::title("Amplituder i tilnærming");
    std::vector<int> freq;
    std::vector<double> freq_amplitudes;
    freq.reserve(approksimasjon.cosbølger.size());
    freq_amplitudes.reserve(approksimasjon.cosbølger.size());

    for (int i = 0; i < approksimasjon.cosbølger.size() && i < 60; ++i) {
        Cosbølge& j = approksimasjon.cosbølger[i];
        freq.push_back(j.frekvens);
        freq_amplitudes.push_back(j.amplitude);
    }

    plt::bar(freq, freq_amplitudes);
    plt::xlabel("n");
    plt::ylabel("Amplitude");
    plt::legend();
    plt::show();
}

void figur5()
{
    double x_start = 0;
    double x_slutt = 5;
    size_t x_num = 11;

    double Len = x_slutt - x_start;

    Cos_approkimasjon approksimasjon;
    std::vector<double> x_utregning = plt::linspace(x_start, x_slutt, x_num);
    std::vector<double> y_utregning = plt::transform(x_utregning, [](double x) {return f5(x); });


    approksimasjon = cosinus_approximasjon(y_utregning, x_start, Len, 12);

    std::vector<double> x = plt::linspace(x_start, x_slutt, 1000);
    std::vector<double> y = plt::transform(x, [](double x) {return f5(x); });
    std::vector<double> y_approx = parse_approximation(approksimasjon, x);

    plt::hold(plt::on);
    plt::title("Graf mot tilnærming");
    plt::plot(x, y)->display_name("f(x) = sin(4 * PI * x)");
    plt::plot(x, y_approx)->display_name("Cosinustilnærming");
    plt::scatter(x_utregning, y_utregning)->display_name("Sampling");
    plt::legend();
    plt::show();
}

void figur6og7()
{
    double x_start = -10;
    double x_slutt = 10;
    size_t x_num = 3'600'000;

    double Len = x_slutt - x_start;

    Cos_approkimasjon approksimasjon;
    std::vector<double> x_utregning = plt::linspace(x_start, x_slutt, x_num);
    std::vector<double> y_utregning = plt::transform(x_utregning, [](double x) {return f6(x); });



    approksimasjon = cosinus_approximasjon_threaded(y_utregning, x_start, Len, 36'000);

    std::vector<double> x = plt::linspace(-10, 10, 5000);
    std::vector<double> y = plt::transform(x, [](double x) {return f6(x); });
    std::vector<double> y_approx = parse_approximation(approksimasjon, x);

    plt::hold(plt::on);
    plt::plot(x, y)->display_name("Step-funksjonen");
    plt::plot(x, y_approx)->display_name("Cosinustilnærming");
    plt::legend();
    plt::show();

    plt::cla();

    auto [X, t] = plt::meshgrid(plt::linspace(x_start, x_slutt, 200), plt::linspace(0, 50, 25));
    auto heat = plt::transform(X, t, [&approksimasjon](double x, double t) {
        return parse_varmeligning_x_t(approksimasjon, x, t, 0.5);
        });

    plt::surf(X, t, heat);
    plt::title("");
    plt::xlabel("x");
    plt::ylabel("Tid");
    plt::zlabel("Temperatur");
    plt::legend(false);
    plt::show();
}