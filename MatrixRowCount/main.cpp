#include <iostream>
#include "mpi.h"
#include <chrono>
#include <random>
#include <string>
using namespace std;

void printMatrix(int size,double* matr)
{
	cout << "Generated matrix:" << endl << endl;
	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; ++j)
			cout << matr[size*i + j] << " | ";
		cout << endl << endl;
	}
}

int main(int argc, char* argv[])
{
	int rank, proc;
	double time = 0.0;
	double* matrix = nullptr;
	double* ans = nullptr;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &proc);

	const int size = stoi(string(argv[1]));
	bool showMatr = false;
	if (argc == 3 && string(argv[2]) == "show")
		showMatr = true;

	if (rank == 0)
	{
		// создаём генератор случайных чисел с seed равным количеству времени с начала эпохи
		default_random_engine generator(chrono::system_clock::now().time_since_epoch().count());
		// создаём равномерное распределение случайной величины типа double в диапазоне 
		//   [-100, 100]
		const uniform_real_distribution <double> distribution(-100, 100);

		matrix = new double[size*size];

		cout.precision(5);
		for (int i = 0; i < size*size; ++i)
			matrix[i] = distribution(generator);

		if (showMatr)
		{
			cout << endl;
			printMatrix(size, matrix);
		}
		/*
		time = MPI_Wtime();
		
		double* ansIter = new double[size];

		for (int i = 0; i < size; ++i)
		{
			double s = 0;
			for (int j = 0; j < size; ++j)
			{
				s += matrix[i*size + j];
			}
			ansIter[i] = s;
		}

		if (showMatr)
		{
			cout << endl << "Sum of elems in each row:" << endl;
			for (int i = 0; i < size; ++i)
			{
				cout << ansIter[i] << " | ";
			}
			cout << endl;
		}
		cout << endl << "Iterative algo time: " << MPI_Wtime() - time << "sec" << endl << endl;
		*/
		time = MPI_Wtime();
	}

	const int nRows = size / proc;

	auto mas = new double[size*nRows];

	MPI_Scatter(matrix, nRows*size, MPI_DOUBLE, mas, nRows*size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	auto res = new double[nRows];

	for (int i = 0; i < nRows; ++i)
	{
		double s = 0.0;
		for (int j = 0; j < size; ++j)
		{
			s += mas[i*size + j];
		}
		res[i] = s;
	}

	if (rank == 0)
		ans = new double[size];

	MPI_Gather(res, nRows, MPI_DOUBLE, ans, nRows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(rank == 0)
	{
		if (showMatr)
		{
			cout << endl << "Sum of elems in each row:" << endl;
			for (int i = 0; i < size; ++i)
			{
				cout << ans[i] << " | ";
			}
			cout << endl;
		}
		time = MPI_Wtime() - time;
		cout << endl << "Parallel algo time: " << time << "sec" << endl;
	}
	MPI_Finalize();
	return 0;
}
