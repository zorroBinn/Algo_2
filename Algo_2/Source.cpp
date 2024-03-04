#include <iostream>
#include <fstream>
#include <conio.h>
#include <Windows.h>
#include <cmath>
#include <float.h>
#include <iomanip>
using namespace std;
#define MAX_ITERATION 10000 //Максимальное кол-во итераций для методов
#define RAND_MAX_ABS 10 //Максимальное значение (по модулю) генерируемого числа

void ArraySelecting(int size, double**& A_Matrix, double*& B, double*& Xj, double*& Xs) {
    A_Matrix = new double* [size];
    for (int i = 0; i < size; i++) {
        A_Matrix[i] = new double[size];
    }
    B = new double[size];
    Xj = new double[size];
    Xs = new double[size];
}

void RandomMatrix(int size, double** A_Matrix, double* B) {
    srand(time(0));
    double tmp = 0;
    int sgn;
    for (int i = 0; i < size; i++) {
        B[i] = -RAND_MAX_ABS + rand() % (2 * RAND_MAX_ABS + 1);
        for (int j = 0; j < size; j++) {
            if (i != j) {
                A_Matrix[i][j] = -RAND_MAX_ABS + rand() % (2 * RAND_MAX_ABS + 1);
                tmp += fabs(A_Matrix[i][j]);
            }
        }
        sgn = 0 + rand() % 2;
        A_Matrix[i][i] = tmp + (1 + rand() % RAND_MAX_ABS);
        if (sgn) A_Matrix[i][i] *= -1;
        tmp = 0;
    }
}

double DiagonalPredominance(int size, double** A_Matrix, int row) {
    double tmp = 0;
    for (int j = 0; j < size; j++) {
        if (j != row) tmp += fabs(A_Matrix[row][j]);
    }
    return fabs(A_Matrix[row][row] / tmp);
}

double DiscrepanciesNorm(int size, double** A_Matrix, double* B, double* X) {
    double* discrepancies = new double[size]; //Вектор невязок
    double tmp = 0.0;
    for (int i = 0; i < size; i++) {
        discrepancies[i] = B[i];
        for (int j = 0; j < size; j++) {
            discrepancies[i] -= A_Matrix[i][j] * X[j];
        }
    }
    for (int i = 0; i < size; i++) {
        tmp += discrepancies[i] * discrepancies[i];
    }
    return sqrt(tmp);
}

int JacobiMethod(int size, double eps, double** A_Matrix, double* B, double* Xj, int& numberOfIteration) {
    double* TempXj = new double[size];  
    do {
        for (int i = 0; i < size; i++) {
            TempXj[i] = B[i];
            for (int j = 0; j < size; j++) {
                if (i != j) TempXj[i] -= A_Matrix[i][j] * Xj[j];
            }
            TempXj[i] /= A_Matrix[i][i];
        }
        for (int i = 0; i < size; i++) {
            Xj[i] = TempXj[i];
        }
        numberOfIteration++;
        if (numberOfIteration > MAX_ITERATION) return 1;
    } while (DiscrepanciesNorm(size, A_Matrix, B, Xj) > eps);
    delete[] TempXj;
    return 0;
}

int SeidelMethod(int size, double eps, double** A_Matrix, double* B, double* Xs, int& numberOfIteration) {
    do {
        for (int i = 0; i < size; i++) {
            double tmp = 0.0;
            for (int j = 0; j < size; j++) {
                if (i != j) tmp += A_Matrix[i][j] * Xs[j];
            }
            Xs[i] = (B[i] - tmp) / A_Matrix[i][i];
        }
        numberOfIteration++;
        if (numberOfIteration > MAX_ITERATION) return 1;
    } while (DiscrepanciesNorm(size, A_Matrix, B, Xs) > eps);
    return 0;
}

void ReadingFromFile(int& size, double& eps, double& initialApproximation, double**& A_Matrix, double*& B, double*& Xj, double*& Xs) {
    char temp; //Переменная для проверки файла
    //Чтение из файла
    try {
        ifstream fin("input.txt");
        if (!fin.is_open()) { // Проверяем, удалось ли открыть файл
            throw string("Не удалось открыть файл input.txt");
        }

        fin >> size; //Чтение размерности матрицы коэффициентов

        if (size < 2) { //Если не число или не положительное - вызов исключения
            throw string("Проблема с данными в файле");
        }

        fin >> eps; //Чтение погрешности 

        if (eps <= 0) { //Если не число или не положительное - вызов исключения
            throw string("Проблема с данными в файле");
        }

        fin >> initialApproximation; //Чтение начального приближения

        ArraySelecting(size, A_Matrix, B, Xj, Xs);

        //Чтение матрицы коэффициентов
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (!(fin >> A_Matrix[i][j]) || fin.eof()) { //Если не число или достигут конец файла - вызов исключения
                    throw string("Проблема с данными в файле");
                }
            }
        }

        //Чтение столбца свободных сленов
        for (int i = 0; i < size; i++) {
            if (!(fin >> B[i])) {
                throw string("Проблема с данными в файле. В данных не число"); //Если не число - вызов исключения
            }
            if (fin.eof() && i != size - 1) { //Если достигнут конец файла до завершения запонения - вызов исключения
                throw string("Проблема с данными в файле. Размерность указана меньше действительной");
            }
            if (i == size - 1) { //Если в файле лишние данные или размерность матрицы указана меньше действительной - вызов исключения
                if (fin >> temp)
                    throw string("В файле лишние данные или размерность матрицы указана меньше действительной");
            }
        }
        fin.close();
    }
    catch (string err_message) {
        cerr << err_message << endl;
        _getch();
        exit(1);
    }
}

int ZeroOnDiagonal(int size, double**& A_Matrix, double*& B) {
    for (int i = 0; i < size; i++) {
        if (fabs(A_Matrix[i][i]) < DBL_EPSILON) return 1;
    }
}

void WritingToFile(int size, double eps, double** A_Matrix, double* B, double* Xj, double* Xs, double initialApproximation, int& numberOfIteration) {
    try {
        ofstream fout("output.txt");
        if (!fout.is_open()) { // Проверяем, удалось ли открыть файл
            throw string("Ошибка при создании файла output.txt, результаты не могут быть записаны");
        }

        fout << "Матрица коэффициентов: " << endl;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                fout << A_Matrix[i][j] << "\t";
            }
            fout << endl;
        }

        fout << endl << "Столбец свободных членов:" << endl;
        for (int i = 0; i < size; i++) {
            fout << B[i] << endl;
        }

        fout << endl << "Коэффициенты диагонального преобладания:" << endl;
        for (int i = 0; i < size; i++) {
            fout << setprecision(3) << DiagonalPredominance(size, A_Matrix, i) << " ";
        }
        fout << setprecision(5);

        fout << endl << "Заданная погрешность e: " << eps << endl;

        fout << "Начальное приближение: " << initialApproximation << endl;

        for (int i = 0; i < size; i++) {
            Xj[i] = initialApproximation;
        }

        fout << endl << "Решение методом Якоби:" << endl;
        if (JacobiMethod(size, eps, A_Matrix, B, Xj, numberOfIteration) == 1) {
            fout << "Превышен максимум итераций метода Якоби (10000), нет схождения";
        }
        else {
            for (int i = 0; i < size; i++) {
                fout << Xj[i] << endl;
            }
            fout << "Число итераций метода Якоби:: " << numberOfIteration << endl;
            fout << "Норма невязок: " << DiscrepanciesNorm(size, A_Matrix, B, Xj);
        }

        numberOfIteration = 0;
        for (int i = 0; i < size; i++) {
            Xs[i] = initialApproximation;
        }

        fout << endl << "Решение методом Зейделя:" << endl;
        if (SeidelMethod(size, eps, A_Matrix, B, Xs, numberOfIteration) == 1) {
            fout << "Превышен максимум итераций метода Зейделя (10000), нет схождения";
        }
        else {
            for (int i = 0; i < size; i++) {
                fout << Xs[i] << endl;
            }
            fout << "Число итераций: " << numberOfIteration << endl;
            fout << "Норма невязок: " << DiscrepanciesNorm(size, A_Matrix, B, Xs);
        }

        fout.close();
        cout << endl << "Программа успешно завершена. Результаты записаны в файл output.txt" << endl;
    }
    catch (string err_message) {
        cerr << err_message << endl;
    }
}

int main() {
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
    double** A_Matrix; //Двумерный массив для работы с матрицей коэффициентов
    double* B; //Массив для работы с матрицей свободных членов
    double* Xj; //Vассив для решений системы методом Якоби
    double* Xs; //Vассив для решений системы методом Зейделя
    int size; //Переменная - размерность матрицы коэффициентов
    double eps; //Погрешность
    double initialApproximation; //Начальное приближение
    int numberOfIteration = 0; //Кол-во итераций
    string IsRandom, method;

    do {
        system("cls");
        cout << "Автоматически генерировать элементы матрицы?" << endl;
        cout << "1 - Да" << endl << "2 - Считать матрицу из файла" << endl;
        cin >> IsRandom;
    } while (IsRandom != "1" && IsRandom != "2");
    
    if (IsRandom == "1") {
        do {
            cout << "Введите размерность матрицы: ";
            cin >> size;
        } while (size < 2);
        do {
            cout << "Введите погрешность e: ";
            cin >> eps;
        } while (eps <= 0);
        cout << "Введите начальное приближение: ";
        cin >> initialApproximation;
        
        ArraySelecting(size, A_Matrix, B, Xj, Xs);
        RandomMatrix(size, A_Matrix, B);
        if (ZeroOnDiagonal(size, A_Matrix, B) != 1) {
            WritingToFile(size, eps, A_Matrix, B, Xj, Xs, initialApproximation, numberOfIteration);
        }
    }
    
    if (IsRandom == "2") {
        ReadingFromFile(size, eps, initialApproximation, A_Matrix, B, Xj, Xs);
        if (ZeroOnDiagonal(size, A_Matrix, B) != 1) {
            WritingToFile(size, eps, A_Matrix, B, Xj, Xs, initialApproximation, numberOfIteration);
        }
    }

    _getch();
}