#include"opt_alg.h"
#include <cmath>

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

solution expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution interval;
		interval.x = matrix(2, 1);
		int i = 0;
		solution X0 = x0;
		solution X1 = X0.x + d;

		if (X1.fit_fun(ff) == X0.fit_fun(ff))
		{
			interval.x(0, 0) = X0.x(0);
			interval.x(1, 0) = X1.x(0);
			return interval;
		}

		if (X1.fit_fun(ff) > X0.fit_fun(ff))
		{
			d = -d;
			X1 = X0.x + d;
			if (X1.fit_fun(ff) >= X0.fit_fun(ff))
			{
				interval.x(0, 0) = X1.x(0);
				interval.x(1, 0) = X0.x(0) - d;
				return interval;
			}
		}

		solution xPre;

		while (true)
		{
			if (i >= Nmax) throw "Maximum number of function calls exceeded";

			i++;
			xPre = X0;
			X0 = X1;
			X1.x(0) = X0.x(0) + pow(alpha, i) * d;

			if (X0.fit_fun(ff) <= X1.fit_fun(ff))
				break;
		}

		if (d > 0)
		{
			interval.x(0, 0) = xPre.x(0);
			interval.x(1, 0) = X1.x(0);
		}
		else
		{
			interval.x(0, 0) = X1.x(0);
			interval.x(1, 0) = xPre.x(0);
		}

		return interval;
	}
	catch (string ex_info)
	{
		throw ("solution expansion(...):\n" + ex_info);
	}
}
vector<double> generateFibonacci(int n) {
	vector<double> fib(n + 1);
	fib[0] = 1;
	fib[1] = 1;
	for (int i = 2; i <= n; ++i) {
		fib[i] = fib[i - 1] + fib[i - 2];
	}
	return fib;
}
solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		int k = 0; 
		double fibKMinus2 = 1; 
		double fibKMinus1 = 1; 
		double fibK = fibKMinus2 + fibKMinus1;

		while (fibK <= (b - a) / epsilon)
		{
			fibKMinus2 = fibKMinus1;
			fibKMinus1 = fibK;
			fibK = fibKMinus1 + fibKMinus2;
			k++;
		}
		vector<double> fib_sequence = generateFibonacci(k);
		solution A0 = a;
		solution B0 = b;
		solution C0 = b - (fibKMinus1 / fibK) * (b - a);
		solution D0 = a + b - C0.x(0);
		solution A1, B1, C1, D1;

		for (int i = 0; i <= k - 3; i++)
		{
			if (C0.fit_fun(ff) < D0.fit_fun(ff))
			{
				A1 = A0; 
				B1 = D0; 
			}
			else
			{
				B1 = B0; 
				A1 = C0; 
			}
			
			C1 = B1.x(0) - (fib_sequence[k - i - 2] / fib_sequence[k - i - 1]) * (B1.x(0) - A1.x(0));
			D1 = A1.x + B1.x - C1.x;

			
		}

		Xopt = C1; 
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		solution A0 = a;
		solution B0 = b;
		solution C0 = 0.5 * (A0.x + B0.x);
		solution D0;
		solution D1, A1, B1, C1;
		solution l, m;
		int i = 0;

		while (true)
		{
			if (i > Nmax)
				throw "Maximum number of function calls exceeded";

			l = A0.fit_fun(ff) * (pow(B0.x, 2) - pow(C0.x, 2)) +
				B0.fit_fun(ff) * (pow(C0.x, 2) - pow(A0.x, 2)) +
				C0.fit_fun(ff) * (pow(A0.x, 2) - pow(B0.x, 2));
			m = A0.fit_fun(ff) * (B0.x - C0.x) +
				B0.fit_fun(ff) * (C0.x - A0.x) +
				C0.fit_fun(ff) * (A0.x - B0.x);

			if (m.x <= 0) 
			{
				Xopt.x(0) = 666;
				Xopt.y = 666;
				return Xopt;
			}
				//throw "Error: m is less than or equal to zero, cannot proceed";
				

			D0 = 0.5 * l.x / m.x;

			if (A0.x < D0.x && D0.x < C0.x)
			{
				if (D0.fit_fun(ff) < C0.fit_fun(ff))
				{
					A1 = A0;
					B1 = C0;
					C1 = D0;
				}
				else
				{
					A1 = D0;
					C1 = C0;
					B1 = B0;
				}
			}

			else if (C0.x < D0.x && D0.x < B0.x)
			{
				if (D0.fit_fun(ff) < C0.fit_fun(ff))
				{
					A1 = C0;
					C1 = D0;
					B1 = B0;
				}
				else
				{
					A1 = A0;
					C1 = C0;
					B1 = D0;
				}
			}
			else
			{
				//throw "Error: D(i) is out of bounds";
				Xopt.x(0) = 666;
				Xopt.y(0) = 666;
				return Xopt;
			}

			if ((B0.x - A0.x < epsilon) || (fabs(D0.x(0) - D1.x(0)) < gamma))
				break;
			D1.x(0) = D0.x(0);
			A0 = A1;
			B0 = B1;
			C0 = C1;
			D1 = D0;
			i++;
		}

		Xopt = D0;
		Xopt.flag = 0;
		return Xopt;
	}

	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
