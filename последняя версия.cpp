#include <iostream>
#include <cstdio>
#include <cmath>
#include <chrono>
#include <time.h>
#include <fstream>
using namespace std;

void DGEMM_BLAS_1(double *X, double *Y, double *R, int size, double& time) {
	int i, j, k;
	auto start = chrono::steady_clock::now();
	for (i = 0; i < size; i++)
		for (j = 0; j < size; j++) {

			for (k = 0; k < size; k++)
				R[i * size + j] += X[i * size + k] * Y[k * size + j];
		}
	auto end = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    time = elapsed_seconds.count();
}

void DGEMM_BLAS_2(double** X, double** Y, double** R, int N, double& time)
{	int i,j,k;
    auto start = chrono::steady_clock::now();
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
        {
            for (k = 0; k < N; k++)
                R[i][j] += X[i][k] * Y[k][j];
        }
	auto end = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    time = elapsed_seconds.count();
}

void DGEMM_opt_1_1(double* X, double* Y, double *R, int N, double& time)
{   int i, j, k;
    auto start = chrono::steady_clock::now();
    for (i = 0; i < N; i++) 
	for (k = 0; k < N; k++) 
		for (j = 0; j < N; j++) 
			R[i * N + j] += X[i * N + k] * Y[k * N + j];
    auto end = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    time = elapsed_seconds.count();
}

void DGEMM_opt_1_2(double** X, double** Y, double **R, int N, double& time)
{
    auto start = chrono::steady_clock::now();
    for (int i=0;i<N;i++)
        for (int k=0;k<N;k++)
            for (int j=0;j<N;j++)
                R[i][j]+=X[i][k]*Y[k][j];
    auto end = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    time = elapsed_seconds.count();
}


void DGEMM_opt_2(int size, double *X, double *Y, double *R, int size_block, double& time)
{
    int i, j, k, ik, jk, kk;
    for(i = 0; i < size; i++)
    {
        for(j = 0; j < size; j++)
        {
            R[i * size + j] = 0;
        }
    }
    auto start = chrono::steady_clock::now();
    for(ik = 0; ik < size; ik+= size_block)
        for(kk = 0; kk < size; kk+= size_block)
            for(jk = 0; jk < size; jk+= size_block)
                for(i = 0; i < size_block; i++ )
                    for(k = 0; k < size_block; k++ )
#pragma simd
                        for(j = 0; j < size_block; j++ )
                            R[(ik + i) * size + (jk + j)] += X[(ik + i) * size + (kk + k)] * Y[(kk + k) * size + (jk + j)];
    auto end = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    time = elapsed_seconds.count();
}

void DGEMM_opt_3(int n, double *A, double *B, double *C, double& time)
{
    int i, j, k;
    auto start = chrono::steady_clock::now();
    for(i = 0; i < n; i++ )
        for(j = 0; j < n; j++ )
            C[i * n + j] = 0;
    for(i = 0; i < n; i++ )
        for(k = 0; k < n; k++ )
// делает векторизацию кода
#pragma simd
            for(j = 0; j < n; j++ )
                C[i * n + j] += A[i * n + k] * B[k * n + j];
    auto end = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    time = elapsed_seconds.count();
}


int main(int argc, char *argv[])
{

    srand(time(nullptr));

    double *X, *Y, *R;
    double **A, **B, **C;
    int size = 500;
    double someTime1, sumTime1, averageTime1, someTime2, sumTime2, averageTime2;
    int choice = 2;
    double speedBLAS[4];

    int i, j, k, l, p, w;
	
	ofstream fout;
	fout.open("result.csv");
	
    switch (choice)
    {
    case 0:
        for(k = 0; k < 4; k++)
        {
            sumTime1 = 0;
            averageTime1 = 0;
			sumTime2 = 0;
            averageTime2 = 0;
            for(p = 0; p < 20; p++)
            {
               X = new double[size * size];
	Y = new double[size * size];
	R = new double[size * size];
	for ( i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			X[i * size + j] = rand() % 10 + 0.000001 * (rand() % 100000);
			Y[i * size + j] = rand() % 10 + 0.000001 * (rand() % 100000);
			R[i * size + j] = 0;
		}
	}
                A = new double*[size];
                B = new double*[size];
                C = new double*[size];
                for(i = 0; i < size; i++)
                {
                    A[i] = new double[size];
                    B[i] = new double[size];
                    C[i] = new double[size];
                    for(j = 0; j < size; j++)
                    {
                        A[i][j] = X[i*size+j];
                        B[i][j] = Y[i*size+j];
                        C[i][j] = 0;
                    }
                }
                DGEMM_BLAS_1(X, Y, R, size, someTime1);
                cout << "№" << p + 1 << "DGEMM_BLAS_1" << size << " " << someTime1 << endl;
                sumTime1 += someTime1;
				DGEMM_BLAS_2(A, B, C, size, someTime2);
				cout << "№" << p + 2 << "DGEMM_BLAS_2" << size << " " << someTime2 << endl;
				sumTime1 += someTime2;
                 delete[]X;
                 delete[]Y;
                 delete[]R;
				for(i = 0; i < size; i++)
                {
                    delete[]A[i];
                    delete[]B[i];
                    delete[]C[i];
					
                }
                delete[]A;
                delete[]B;
                delete[]C;
            }
            averageTime1 = sumTime1/10;
            cout << "averageTime1 = " << averageTime1 << endl;
			fout << "DGEMM_BLAS_1;" << size << ";" << averageTime1 << "\n";
			averageTime2 = sumTime2/10;
            cout << "averageTime2 = " << averageTime2 << endl;
			fout << "DGEMM_BLAS_2;" << size << ";" << averageTime2 << "\n";
            size += 500 + 500*k;
        }
        break;

    case 1:
        for(k = 0; k < 4; k++)
        {
            sumTime1 = 0;
            averageTime1 = 0;
			sumTime2 = 0;
            averageTime2 = 0;
	
             for(p = 0; p < 20; p++)
            {
               X = new double[size * size];
			   Y = new double[size * size];
			   R = new double[size * size];
	
	for ( i = 0; i < size; i++) {
		for ( j = 0; j < size; j++) {
			X[i * size + j] = rand() % 10 + 0.000001 * (rand() % 100000);
			Y[i * size + j] = rand() % 10 + 0.000001 * (rand() % 100000);
			R[i * size + j] = 0;
		}
	}
                A = new double*[size];
                B = new double*[size];
                C = new double*[size];
                for(i = 0; i < size; i++)
                {
                    A[i] = new double[size];
                    B[i] = new double[size];
                    C[i] = new double[size];
                    for(j = 0; j < size; j++)
                    {
                        A[i][j] = X[i*size+j];
                        B[i][j] = Y[i*size+j];
                        C[i][j] = 0;
                    }
                }
                DGEMM_opt_1_1(X, Y, R, size, someTime1);
                cout << "№" << p + 1 << "DGEMM_opt_1_1 " << size << " " << someTime1 << endl;
                sumTime1 += someTime1;
				DGEMM_opt_1_2(A, B, C, size, someTime2);
                cout << "№" << p + 2 << "DGEMM_opt_1_1 " << size << " " << someTime2 << endl;
                sumTime2 += someTime2;
				 delete[]X;
                 delete[]Y;
                 delete[]R;
                for(i = 0; i < size; i++)
                {
                    delete[]A[i];
                    delete[]B[i];
                    delete[]C[i];
                }
                delete[]A;
                delete[]B;
                delete[]C;
            }
            averageTime1 = sumTime1/10;
			averageTime2 = sumTime2/10;
            cout << "averageTime1 = " << averageTime1 << endl;
			cout << "averageTime2 = " << averageTime2 << endl;
             cout << "Speed: x" << averageTime/speedBLAS[k] << endl;
			fout << "DGEMM_opt_1_1;" << size << ";" << averageTime1 << "\n";
			fout << "DGEMM_opt_1_2;" << size << ";" << averageTime2 << "\n";
            size += 500 + 500*k;
        }
        break;

    case 2:
        for(k = 0; k < 4; k++)
        {
            int blockSize[5] = { 4, 20, 50, 100, 250 };
            double blockTime[5], someTime1;  
            int id = 0;
            for(l = 0; l < 5; l++)
            {
                sumTime1 = 0;
                averageTime1 = 0;

                for(p = 0; p < 10; p++)
                {
                    X = new double[size*size];
                    Y = new double[size*size];
                    R = new double[size*size];
                    for(i = 0; i < size; i++)
                    {
                        for(j = 0; j < size; j++)
                        {
                            X[i*size + j] = rand()%101 - 50;
                            Y[i*size + j] = rand()%101 - 50;
                            R[i*size + j] = 0;
                        }
                    }
                    DGEMM_opt_2(size, X, Y, R, blockSize[l], someTime1);
                    cout << "№" << p + 1 << "DGEMM_opt_2 " << size << " " << someTime1 << " Block size: " << blockSize[l] << endl;
                    sumTime1 += someTime1;
                    delete[]X;
                    delete[]Y;
                    delete[]R;
                }
                averageTime1 = sumTime1/10;
                blockTime[l] = averageTime1;
                cout << "averageTime = " << averageTime1 << endl;
				fout << "DGEMM_opt_2;" << size << ";" << blockSize[l] << ";" << averageTime1 << "\n";
                // cout << "Speed: x" << averageTime/speedBLAS[k] << endl;               
            }
            for(l = 0; l < 5; l++)
            {
                if(blockTime[l] < blockTime[id])
                    id = l;
            }
            cout << "Оптимальный Block size для матрицы размерности " << size << ": " << blockSize[id] << endl << endl;
			size += 500 + 500*k;
        }
        break;

    case 3:
        for(k = 0; k < 4; k++)
        {
            sumTime1 = 0;
            averageTime1 = 0;

            for(p = 0; p < 10; p++)
            {
                X= new double[size*size];
                Y = new double[size*size];
                R = new double[size*size];
                for(i = 0; i < size; i++)
                {
                    for(j = 0; j < size; j++)
                    {
                        X[i*size + j] = rand()%101 - 50;
                        Y[i*size + j] = rand()%101 - 50;
                        R[i*size + j] = 0;
                    }
                }
                DGEMM_opt_3(size, X, Y, R, someTime1);
                cout << "№" << p + 1 << "DGEMM_opt_3 " << size << " " << someTime1 << endl;
                sumTime1 += someTime1;
                delete[]X;
                delete[]Y;
                delete[]R;
            }
            averageTime1 = sumTime1/10;
            cout << "averageTime = " << averageTime1 << endl;
            // cout << "Speed: x" << averageTime/speedBLAS[k] << endl;
			fout << "DGEMM_opt_3;" << size << ";" << averageTime1 << "\n";
            size += 500 + 500*k;
        }
        break;
    }

	fout.close();

    return 0;
}
