#include <iostream>
#include <ctime>
#include <fstream>
#include <string>
#include <chrono>
#include <windows.h>


class kholetskyDecomposition {
private:
    int n;
    double** inputArr;
    double** outputArr;
    double** exactMatrix;
    double operatingTime;

    double** matrixTranspose(double** arr) {
        double** newArray = createArr(n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                newArray[j][i] = arr[i][j];
            }
        }
        return newArray;
    }

    double** matrixMultiplication(double** A, double** B) {
        double** C = createArr(n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                for (int k = 0; k < n; k++) {
                    C[i][j] += A[i][k] * B[k][j];
                }
        return C;
    }

    double** createArr(int N) {
        double** array = new double* [N];
        for (int i = 0; i < N; ++i) {
            array[i] = new double[N];
            memset(array[i], 0, sizeof(double) * N);
        }
        return array;
    }

    void deleteArr(double** array, int N) {
        for (int i = 0; i < N; ++i) {
            delete[] array[i];
        }
        delete[] array;
    }

public:
    kholetskyDecomposition(int N) {
        n = N;
        inputArr = createArr(n);
        outputArr = createArr(n);
        exactMatrix = createArr(n);
    }

    ~kholetskyDecomposition() {
        deleteArr(inputArr, n);
        deleteArr(outputArr, n);
        deleteArr(exactMatrix, n);
    }

    double getTime() {
        return operatingTime;
    }

    void matrixGenerator(double lowerLimit = -10, double upperLimit = 10) {
        srand(time(0));

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                double v = lowerLimit + ((double)rand() / (double)(RAND_MAX)) * (upperLimit - lowerLimit);
                exactMatrix[i][j] = v;
            }
        }
        double** transposedExactMatrix = matrixTranspose(exactMatrix);
        inputArr = matrixMultiplication(exactMatrix, transposedExactMatrix);

        deleteArr(transposedExactMatrix, n);
    }
    
    void optimized_decompose() {
        auto begin = std::chrono::steady_clock::now();
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                double res = inputArr[i][j];
                for (int k = 0; k < j; ++k) {
                    res -= outputArr[i][k] * outputArr[j][k];
                }

                if (i == j) {
                    outputArr[i][j] = sqrt(res);
                } else {
                    outputArr[i][j] = res / outputArr[j][j];
                }
            }
        }
        auto end = std::chrono::steady_clock::now();
        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
        operatingTime = elapsed_ms.count();
    }
};


double algorithm(int dimension, double lowerLimit, double upperLimit) {
    kholetskyDecomposition A(dimension);
    A.matrixGenerator(lowerLimit, upperLimit);
    A.optimized_decompose();
    return A.getTime();
}

void iteration(int dimension, int maxDimention, double multiplicity = 2, double lowerLimit = -10, double upperLimit = 10,
               int precision = 5, std::string outfile = "outfile.txt") {

    std::ofstream ofile;
    ofile.open(outfile);
    std::cout.precision(precision);
    int n;
    double operatingTime;
    for ( ; dimension < maxDimention; dimension *= multiplicity) {
        n = dimension;
        operatingTime = algorithm(n, lowerLimit, upperLimit);
        ofile << n;
        ofile << "\\" << operatingTime << std::endl;

        n = dimension + (dimension / 4);
        operatingTime = algorithm(n, lowerLimit, upperLimit);
        ofile << n;
        ofile << "\\" << operatingTime << std::endl;

        n = dimension + (dimension / 2);
        operatingTime = algorithm(n, lowerLimit, upperLimit);
        ofile << n;
        ofile << "\\" << operatingTime << std::endl;

        n = dimension + 3 * (dimension / 4);
        operatingTime = algorithm(n, lowerLimit, upperLimit);
        ofile << n;
        ofile << "\\" << operatingTime << std::endl;
    }
    n = dimension;
    operatingTime = algorithm(n, lowerLimit, upperLimit);
    ofile << n;
    ofile << "\\" << operatingTime << std::endl;
    ofile.close();
}

int main() {
    for (int i = 0; i < 10; ++i) {
        iteration(128, pow(2, 12), 2, -pow(10, -1), pow(10, 1), 5, "outfile-" + std::to_string(i) + ".txt");
    }
}
