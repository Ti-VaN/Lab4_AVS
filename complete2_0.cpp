#include <iostream>
#include <cstdio>
#include <cmath>
#include <chrono>
#include <time.h>
#include <fstream>
using namespace std;

void DGEMM_BLAS(double** X, double** Y, double** R, int N, double& time)
{
    auto start = chrono::steady_clock::now();
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < N; k++)
                R[i][j] += X[i][k] * Y[k][j];
        }
    auto end = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    time = elapsed_seconds.count();
}

void DGEMM_opt_1(double** X, double** Y, double **R, int N, double& time1)
{
    auto start = chrono::steady_clock::now();
    for (int i=0;i<N;i++)
        for (int j=0;j<N;j++)
            for (int k=0;k<N;k++)
                R[i][k]+=X[i][j]*Y[j][k];
    auto end = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    time1 = elapsed_seconds.count();
}

void DGEMM_opt_2(int size, double *A, double *B, double *C, int size_block, double& time1)
{
    int i, j, k, ik, jk, kk;
    for(j = 0; j < size; j++)
    {
        for(i = 0; i < size; i++)
        {
            C[j * size + i] = 0;
        }
    }
    auto start = chrono::steady_clock::now();
    for(jk = 0; jk < size; jk+= size_block)
        for(kk = 0; kk < size; kk+= size_block)
            for(ik = 0; ik < size; ik+= size_block)
                for(j = 0; j < size_block; j++ )
                    for(k = 0; k < size_block; k++ )
#pragma simd
                        for(i = 0; i < size_block; i++ )
                            C[(jk + j) * size + (ik + i)] += A[(jk + j) * size + (kk + k)] * B[(kk + k) * size + (ik + i)];
    auto end = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    time1 = elapsed_seconds.count();
}

void DGEMM_opt_3(int n, double *A, double *B, double *C, double& time1)
{
    int i, j, k;
    auto start = chrono::steady_clock::now();
    for(j = 0; j < n; j++ )
        for(i = 0; i < n; i++ )
            C[j * n + i] = 0;
    for(j = 0; j < n; j++ )
        for(k = 0; k < n; k++ )
// делает векторизацию кода
#pragma simd
            for(i = 0; i < n; i++ )
                C[j * n + i] += A[j * n + k] * B[k * n + i];
    auto end = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    time1 = elapsed_seconds.count();
}


int main(int argc, char *argv[])
{

    srand(time(nullptr));

    double **X, **Y, **R;
    double *A, *B, *C;
    int size = 1000;
    double someTime, sumTime, averageTime;
    int choice = 3;
    double speedBLAS[4];

    int i, j, k, l, p, w;
	
	ofstream fout;
	fout.open("result.csv");
	
    switch (choice)
    {
    case 0:
        for(k = 0; k < 4; k++)
        {
            sumTime = 0;
            averageTime = 0;

            for(p = 0; p < 10; p++)
            {
                X = new double*[size];
                Y = new double*[size];
                R = new double*[size];
                for(i = 0; i < size; i++)
                {
                    X[i] = new double[size];
                    Y[i] = new double[size];
                    R[i] = new double[size];
                    for(j = 0; j < size; j++)
                    {
                        X[i][j] = rand()%101 - 50;
                        Y[i][j] = rand()%101 - 50;
                        R[i][j] = 0;
                    }
                }
                DGEMM_BLAS(X, Y, R, size, someTime);
                cout << "№" << p + 1 << "DGEMM_BLAS " << size << " " << someTime << endl;
                sumTime += someTime;
                for(i = 0; i < size; i++)
                {
                    delete[]X[i];
                    delete[]Y[i];
                    delete[]R[i];
                }
                delete[]X;
                delete[]Y;
                delete[]R;
            }
            averageTime = sumTime/10;
            cout << "averageTime = " << averageTime << endl;
			fout << "DGEMM_BLAS;" << size << ";" << averageTime << "\n";
            size += 1000 + 2000*k;
        }
        break;

    case 1:
        for(k = 0; k < 4; k++)
        {
            sumTime = 0;
            averageTime = 0;

            for(p = 0; p < 10; p++)
            {
                X = new double*[size];
                Y = new double*[size];
                R = new double*[size];
                for(i = 0; i < size; i++)
                {
                    X[i] = new double[size];
                    Y[i] = new double[size];
                    R[i] = new double[size];
                    for(j = 0; j < size; j++)
                    {
                        X[i][j] = rand()%101 - 50;
                        Y[i][j] = rand()%101 - 50;
                        R[i][j] = 0;
                    }
                }
                DGEMM_opt_1(X, Y, R, size, someTime);
                cout << "№" << p + 1 << "DGEMM_opt_1 " << size << " " << someTime << endl;
                sumTime += someTime;
                for(i = 0; i < size; i++)
                {
                    delete[]X[i];
                    delete[]Y[i];
                    delete[]R[i];
                }
                delete[]X;
                delete[]Y;
                delete[]R;
            }
            averageTime = sumTime/10;
            cout << "averageTime = " << averageTime << endl;
            // cout << "Speed: x" << averageTime/speedBLAS[k] << endl;
			fout << "DGEMM_opt_1;" << size << ";" << averageTime << "\n";
            size += 1000 + 2000*k;
        }
        break;

    case 2:
        for(k = 0; k < 4; k++)
        {
            for(l = 0; l < 5; l++)
            {
                sumTime = 0;
                averageTime = 0;

                for(p = 0; p < 10; p++)
                {
                    A = new double[size*size];
                    B = new double[size*size];
                    C = new double[size*size];
                    for(i = 0; i < size; i++)
                    {
                        for(j = 0; j < size; j++)
                        {
                            A[i*size + j] = rand()%101 - 50;
                            B[i*size + j] = rand()%101 - 50;
                            C[i*size + j] = 0;
                        }
                    }
                    DGEMM_opt_2(size, A, B, C, (int)pow(2, l+3), someTime);
                    cout << "№" << p + 1 << "DGEMM_opt_2 " << size << " " << someTime << " Block size: " << pow(2, l+3) << endl;
                    sumTime += someTime;
                    delete[]A;
                    delete[]B;
                    delete[]C;
                }
                averageTime = sumTime/10;
                cout << "averageTime = " << averageTime << endl;
				fout << "DGEMM_opt_2;" << size << ";" << pow(2, l+3) << ";" << averageTime << "\n";
                // cout << "Speed: x" << averageTime/speedBLAS[k] << endl;               
            }
			size += 1000 + 2000*k;
        }
        break;

    case 3:
        for(k = 0; k < 4; k++)
        {
            sumTime = 0;
            averageTime = 0;

            for(p = 0; p < 10; p++)
            {
                A = new double[size*size];
                B = new double[size*size];
                C = new double[size*size];
                for(i = 0; i < size; i++)
                {
                    for(j = 0; j < size; j++)
                    {
                        A[i*size + j] = rand()%101 - 50;
                        B[i*size + j] = rand()%101 - 50;
                        C[i*size + j] = 0;
                    }
                }
                DGEMM_opt_3(size, A, B, C, someTime);
                cout << "№" << p + 1 << "DGEMM_opt_3 " << size << " " << someTime << endl;
                sumTime += someTime;
                delete[]A;
                delete[]B;
                delete[]C;
            }
            averageTime = sumTime/10;
            cout << "averageTime = " << averageTime << endl;
            // cout << "Speed: x" << averageTime/speedBLAS[k] << endl;
			fout << "DGEMM_opt_3;" << size << ";" << averageTime << "\n";
            size += 1000 + 2000*k;
        }
        break;
    }

	fout.close();

    return 0;
}
