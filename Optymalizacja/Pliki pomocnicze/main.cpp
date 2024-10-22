/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include"opt_alg.h"
#include <random>

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		lab1();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	system("pause");
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;
	int Nmax = 10000;
	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
	solution opt;
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Wahadlo
	Nmax = 1000;
	epsilon = 1e-2;
	lb = 0;
	ub = 5;
	double teta_opt = 1;
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(opt.x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
	ofstream Sout("symulacja_lab0.csv");
	Sout << hcat(Y[0], Y[1]);
	Sout.close();
	Y[0].~matrix();
	Y[1].~matrix();
}

void lab1()
{
    double epsilon = 0.0001;
    double gamma = 0.0001;
    int Nmax = 100;
    double x0;
    double d = 0.5;
    double alpha = 1.5;
    double alpha2 = 3;
    double alpha3 = 5.5;

    std::ofstream Sout_expansion("results_expansion3.csv");
    std::ofstream Sout_fibonacci("results_fibonacci3.csv");
    std::ofstream Sout_lagrange("results_lagrange3.csv");

    Sout_expansion << "x0; a; b; f_calls" << std::endl;
    Sout_fibonacci << "x; y; f_calls" << std::endl;
    Sout_lagrange << "x; y; f_calls" << std::endl;

	std::mt19937 gen(42);
	std::uniform_real_distribution<double> unif(-100.0, 100.0);

    for (int j = 0; j < 100; ++j)
    {
        x0 = unif(gen);

        solution exp = expansion(ff1T, x0, d, alpha3, Nmax);
        matrix interval = exp.x;
        double a = interval(0, 0);
        double b = interval(1, 0);
		
        Sout_expansion << x0 << "; " << a << "; " << b << "; " << solution::f_calls << std::endl;
        solution::clear_calls();

        solution fibonacci = fib(ff1T, a, b, epsilon);
        Sout_fibonacci << fibonacci.x(0) << "; " << ff1T(fibonacci.x) << solution::f_calls << std::endl;
        solution::clear_calls();
		
        solution lagrange = lag(ff1T, a, b, epsilon, gamma, Nmax);
		Sout_lagrange << lagrange.x(0) << "; ";

		if (lagrange.x == 666) {
			Sout_lagrange << 666 << "; ";
		}
		else
			Sout_lagrange << ff1T(lagrange.x);

		Sout_lagrange << solution::f_calls << std::endl;
        solution::clear_calls();
    }

    Sout_expansion.close();
    Sout_fibonacci.close();
    Sout_lagrange.close();
}

void lab2()
{

}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
