#include <iostream>
#include <cstdio>
#include <cmath>
#include <string.h>
#include <chrono>
#include <time.h>
#include <unistd.h>
#include <fstream>
using namespace std;

static const char *optString = "p:";

void DGEMM_BLAS_1(double *X, double *Y, double *R, int size, double &time)
{
    int i, j, k;
    auto start = chrono::steady_clock::now();
    for (i = 0; i < size; i++)
        for (j = 0; j < size; j++)
        {

            for (k = 0; k < size; k++)
                R[i * size + j] += X[i * size + k] * Y[k * size + j];
        }
    auto end = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    time = elapsed_seconds.count();
}

void DGEMM_BLAS_2(double **X, double **Y, double **R, int N, double &time)
{
    int i, j, k;
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

void DGEMM_opt_1(double *X, double *Y, double *R, int N, double &time)
{
    int i, j, k;
    auto start = chrono::steady_clock::now();
    for (i = 0; i < N; i++)
        for (k = 0; k < N; k++)
            for (j = 0; j < N; j++)
                R[i * N + j] += X[i * N + k] * Y[k * N + j];
    auto end = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    time = elapsed_seconds.count();
}

void DGEMM_opt_2(int size, double *X, double *Y, double *R, int size_block, double &time)
{
    int i, j, k, ik, jk, kk;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            R[i * size + j] = 0;
        }
    }
    auto start = chrono::steady_clock::now();
    for (ik = 0; ik < size; ik += size_block)
        for (kk = 0; kk < size; kk += size_block)
            for (jk = 0; jk < size; jk += size_block)
                for (i = 0; i < size_block; i++)
                    for (k = 0; k < size_block; k++)
#pragma simd
                        for (j = 0; j < size_block; j++)
                            R[(ik + i) * size + (jk + j)] += X[(ik + i) * size + (kk + k)] * Y[(kk + k) * size + (jk + j)];
    auto end = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    time = elapsed_seconds.count();
}

void DGEMM_opt_3(int n, double *A, double *B, double *C, double &time)
{
    int i, j, k;
    auto start = chrono::steady_clock::now();
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            C[i * n + j] = 0;
    for (i = 0; i < n; i++)
        for (k = 0; k < n; k++)
// делает векторизацию кода
#pragma simd
            for (j = 0; j < n; j++)
                C[i * n + j] += A[i * n + k] * B[k * n + j];
    auto end = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    time = elapsed_seconds.count();
}

int main(int argc, char *argv[])
{

    srand(time(nullptr));

    int opt = 0;

    opt = getopt(argc, argv, optString);

    double *X, *Y, *R;
    double **A, **B, **C;
    int size = 500;
    double someTime1, sumTime1, averageTime1, someTime2, sumTime2, averageTime2;
    double speedBLAS[4] = {2.067, 15.5995, 141.126, 857.241};

    int i, j, k, l, p, w;

    ofstream fout;
    fout.open("result.csv");

    while (opt != -1)
    {
        switch (opt)
        {
            case 'p':
                if (strcmp(optarg, "BLAS") == 0)
                {
                    for (k = 0; k < 4; k++)
                    {
                        sumTime1 = 0;
                        averageTime1 = 0;
                        sumTime2 = 0;
                        averageTime2 = 0;
                        for (p = 0; p < 10; p++)
                        {
                            X = new double[size * size];
                            Y = new double[size * size];
                            R = new double[size * size];
                            for (i = 0; i < size; i++)
                            {
                                for (j = 0; j < size; j++)
                                {
                                    X[i * size + j] = rand() % 10 + 0.000001 * (rand() % 100000);
                                    Y[i * size + j] = rand() % 10 + 0.000001 * (rand() % 100000);
                                    R[i * size + j] = 0;
                                }
                            }
                            A = new double *[size];
                            B = new double *[size];
                            C = new double *[size];
                            for (i = 0; i < size; i++)
                            {
                                A[i] = new double[size];
                                B[i] = new double[size];
                                C[i] = new double[size];
                                for (j = 0; j < size; j++)
                                {
                                    A[i][j] = X[i * size + j];
                                    B[i][j] = Y[i * size + j];
                                    C[i][j] = 0;
                                }
                            }
                            DGEMM_BLAS_1(X, Y, R, size, someTime1);
                            cout << "№" << p + 1 << " DGEMM_BLAS_1 " << size << " " << someTime1 << endl;
                            sumTime1 += someTime1;
                            DGEMM_BLAS_2(A, B, C, size, someTime2);
                            cout << "№" << p + 1 << " DGEMM_BLAS_2 " << size << " " << someTime2 << endl;
                            sumTime2 += someTime2;
                            delete[] X;
                            delete[] Y;
                            delete[] R;
                            for (i = 0; i < size; i++)
                            {
                                delete[] A[i];
                                delete[] B[i];
                                delete[] C[i];
                            }
                            delete[] A;
                            delete[] B;
                            delete[] C;
                        }
                        averageTime1 = sumTime1 / 10;
                        cout << "averageTime1 = " << averageTime1 << endl;
                        fout << "DGEMM_BLAS_1;" << size << ";" << averageTime1 << "\n";
                        averageTime2 = sumTime2 / 10;
                        cout << "averageTime2 = " << averageTime2 << endl;
                        fout << "DGEMM_BLAS_2;" << size << ";" << averageTime2 << "\n";
                        cout << "Speed: x" << averageTime2 / averageTime1 << endl;
                        size += 500 + 500 * k;
                    }
                }
                else if (strcmp(optarg, "opt_1") == 0)
                {
                    for (k = 0; k < 4; k++)
                    {
                        sumTime1 = 0;
                        averageTime1 = 0;

                        for (p = 0; p < 10; p++)
                        {
                            X = new double[size * size];
                            Y = new double[size * size];
                            R = new double[size * size];

                            for (i = 0; i < size; i++)
                            {
                                for (j = 0; j < size; j++)
                                {
                                    X[i * size + j] = rand() % 10 + 0.000001 * (rand() % 100000);
                                    Y[i * size + j] = rand() % 10 + 0.000001 * (rand() % 100000);
                                    R[i * size + j] = 0;
                                }
                            }
                            DGEMM_opt_1(X, Y, R, size, someTime1);
                            cout << "№" << p + 1 << " DGEMM_opt_1 " << size << " " << someTime1 << endl;
                            sumTime1 += someTime1;

                            delete[] X;
                            delete[] Y;
                            delete[] R;
                        }
                        averageTime1 = sumTime1 / 10;
                        cout << "averageTime1 = " << averageTime1 << endl;
                        cout << "Speed: x" << speedBLAS[k]/averageTime1 << endl;
                        fout << "DGEMM_opt_1;" << size << ";" << averageTime1 << speedBLAS[k]/averageTime1 << ";" << "\n";
                        size += 500 + 500 * k;
                    }
                }
                else if (strcmp(optarg, "opt_2") == 0)
                {
                    for (k = 0; k < 4; k++)
                    {
                        int blockSize[5] = {4, 20, 50, 100, 250};
                        double blockTime[5], someTime1;
                        int id = 0;
                        for (l = 0; l < 5; l++)
                        {
                            sumTime1 = 0;
                            averageTime1 = 0;

                            for (p = 0; p < 10; p++)
                            {
                                X = new double[size * size];
                                Y = new double[size * size];
                                R = new double[size * size];
                                for (i = 0; i < size; i++)
                                {
                                    for (j = 0; j < size; j++)
                                    {
                                        X[i * size + j] = rand() % 10 + 0.000001 * (rand() % 100000);
                                        Y[i * size + j] = rand() % 10 + 0.000001 * (rand() % 100000);
                                        R[i * size + j] = 0;
                                    }
                                }
                                DGEMM_opt_2(size, X, Y, R, blockSize[l], someTime1);
                                cout << "№" << p + 1 << " DGEMM_opt_2 " << size << " " << someTime1 << " Block size: " << blockSize[l] << endl;
                                sumTime1 += someTime1;
                                delete[] X;
                                delete[] Y;
                                delete[] R;
                            }
                            averageTime1 = sumTime1 / 10;
                            blockTime[l] = averageTime1;
                            cout << "averageTime = " << averageTime1 << endl;
                            fout << "DGEMM_opt_2;" << size << ";" << blockSize[l] << ";" << averageTime1 << speedBLAS[k]/averageTime1 << ";" << "\n";
                            cout << "Speed: x" << speedBLAS[k]/averageTime1 << endl;
                        }
                        for (l = 0; l < 5; l++)
                        {
                            if (blockTime[l] < blockTime[id])
                                id = l;
                        }
                        cout << "Оптимальный Block size для матрицы размерности " << size << ": " << blockSize[id] << endl
                            << endl;
                        size += 500 + 500 * k;
                    }
                }
                else if (strcmp(optarg, "opt_3") == 0)
                {
                    for (k = 0; k < 4; k++)
                    {
                        sumTime1 = 0;
                        averageTime1 = 0;

                        for (p = 0; p < 10; p++)
                        {
                            X = new double[size * size];
                            Y = new double[size * size];
                            R = new double[size * size];
                            for (i = 0; i < size; i++)
                            {
                                for (j = 0; j < size; j++)
                                {
                                    X[i * size + j] = rand() % 10 + 0.000001 * (rand() % 100000);
                                    Y[i * size + j] = rand() % 10 + 0.000001 * (rand() % 100000);
                                    R[i * size + j] = 0;
                                }
                            }
                            DGEMM_opt_3(size, X, Y, R, someTime1);
                            cout << "№" << p + 1 << " DGEMM_opt_3 " << size << " " << someTime1 << endl;
                            sumTime1 += someTime1;
                            delete[] X;
                            delete[] Y;
                            delete[] R;
                        }
                        averageTime1 = sumTime1 / 10;
                        cout << "averageTime = " << averageTime1 << endl;
                        cout << "Speed: x" << speedBLAS[k]/averageTime1 << endl;
                        fout << "DGEMM_opt_3;" << size << ";" << averageTime1 << speedBLAS[k]/averageTime1 << ";" << "\n";
                        size += 500 + 500 * k;
                    }
                }
                break;

            default:
                cout << "Wrong argument. Try -p opt_1 or -p opt_2";
                break;
        }

        opt = getopt(argc, argv, optString);
    }

    fout.close();

    return 0;
}
