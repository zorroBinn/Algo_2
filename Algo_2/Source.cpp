#include <iostream>
#include <fstream>
#include <conio.h>
#include <Windows.h>
#include <cmath>
#include <float.h>
using namespace std;
#define MAX_ITERATION 10000 //Максимальное кол-во итераций для методов
#define RAND_MAX_ABS 10 //Максимальное значение (по модулю) генерируемого числа

void ArraySelecting(int size, double**& A_Matrix, double*& B, double*& X) {
    A_Matrix = new double* [size];
    for (int i = 0; i < size; i++) {
        A_Matrix[i] = new double[size];
    }
    B = new double[size];
    X = new double[size];
}

void RandomMatrix(int size, double** A_Matrix, double* B) {
    srand(time(0));
    int tmp = 0;
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

void DiagonalPredominance(int size, double** A_Matrix) {


}

int JacobiMethod(int size, double eps, double** A_Matrix, double* B, double* X, int& numberOfIteration) {
    double* TempX = new double[size];
    double norm;
    do {
        for (int i = 0; i < size; i++) {
            TempX[i] = B[i];
            for (int j = 0; j < size; j++) {
                if (i != j) TempX[i] -= A_Matrix[i][j] * X[j];
            }
            TempX[i] /= A_Matrix[i][i];
        }
        norm = fabs(X[0] - TempX[0]);
        for (int i = 0; i < size; i++) {
            if (fabs(X[i] - TempX[i]) > norm) norm = fabs(X[i] - TempX[i]);
            X[i] = TempX[i];
        }
        numberOfIteration++;
        if (numberOfIteration > MAX_ITERATION) {
            delete[] TempX;
            return 1;
        }
    } while (norm > eps);
    delete[] TempX;
    return 0;
}

int SeidelMethod(int size, double eps, double** A_Matrix, double* B, double* X, int& numberOfIteration) {
    double* TempX = new double[size];
    double norm;
    do {
        for (int i = 0; i < size; i++) {
            TempX[i] = B[i];
            for (int j = 0; j < size; j++) {
                if (i > j) TempX[i] -= A_Matrix[i][j] * X[j];
                else if (i < j) TempX[i] -= A_Matrix[i][j] * TempX[j];
            }
            TempX[i] /= A_Matrix[i][i];
        }
        norm = fabs(X[0] - TempX[0]);
        for (int i = 0; i < size; i++) {
            if (fabs(X[i] - TempX[i]) > norm) norm = fabs(X[i] - TempX[i]);
            X[i] = TempX[i];
        }
        numberOfIteration++;
        if (numberOfIteration > MAX_ITERATION) {
            delete[] TempX;
            return 1;
        }
    } while (norm > eps);
    delete[] TempX;
    return 0;
}

void ReadingFromFile(int& size, double& eps, double& initialApproximation, double**& A_Matrix, double*& B, double*& X) {
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

        fin >> eps; //Чтение размерности матрицы 

        if (eps <= 0) { //Если не число или не положительное - вызов исключения
            throw string("Проблема с данными в файле");
        }

        fin >> initialApproximation; //Чтение начального приближения

        ArraySelecting(size, A_Matrix, B, X);

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

void WritingToFile(int size, double& eps, double** A_Matrix, double* B, double* X) {
    try {
        ofstream fout("output.txt");
        if (!fout.is_open()) { // Проверяем, удалось ли открыть файл
            throw string("Ошибка при создании файла output.txt, результаты не могут быть записаны");
        }

        fout << "Исходная матрица коэффициентов: " << endl;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                fout << A_Matrix[i][j];
            }
            fout << endl;
        }

        fout << endl << "Исходная матрица свободных членов:" << endl;
        for (int i = 0; i < size; i++) {
            fout << B[i] << endl;
        }

        fout << endl << "Заданная точность e = " << eps << endl;

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
    double* B; //Двумерный массив для работы с матрицей свободных членов
    double* X; //Двумерный массив для решений системы
    int size; //Переменная - размерность матрицы коэффициентов
    double eps; //Точность
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
            cout << "Введите точность e: ";
            cin >> eps;
        } while (eps <= 0);
        cout << "Введите начальное приближение: ";
        cin >> initialApproximation;
        
        ArraySelecting(size, A_Matrix, B, X);
        RandomMatrix(size, A_Matrix, B);

        cout << "Сгенгерированная матрица А:" << endl;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                cout << A_Matrix[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl << "Сгенерированный столбец свободных членов:" << endl;
        for (int i = 0; i < size; i++) {
            cout << B[i] << endl;
        }

        for (int i = 0; i < size; i++) {
            X[i] = initialApproximation;
        }

        try {
            cout << endl << "Решение методом Якоби:" << endl;
            if (JacobiMethod(size, eps, A_Matrix, B, X, numberOfIteration) == 1) {
                throw string("Превышен максимум итераций метода Якоби (10000), нет схождения");
            }
            for (int i = 0; i < size; i++) {
                cout << X[i] << endl;
            }
            cout << "Число итераций: " << numberOfIteration;
        }
        catch (string err_message) {
            cerr << err_message << endl;
            _getch();
            exit(1);
        }

        numberOfIteration = 0;
        for (int i = 0; i < size; i++) {
            X[i] = initialApproximation;
        }

        try {
            cout << endl << "Решение методом Зейделя:" << endl;
            if (SeidelMethod(size, eps, A_Matrix, B, X, numberOfIteration) == 1) {
                throw string("Превышен максимум итераций метода Зейделя (10000), нет схождения");
            }
            for (int i = 0; i < size; i++) {
                cout << X[i] << endl;
            }
            cout << "Число итераций: " << numberOfIteration;
        }
        catch (string err_message) {
            cerr << err_message << endl;
            _getch();
            exit(1);
        }
    }
    
    if (IsRandom == "2") {
        
        ReadingFromFile(size, eps, initialApproximation, A_Matrix, B, X);
        
        cout << "Матрица А:" << endl;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                cout << A_Matrix[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl << "Столбец свободных членов:" << endl;
        for (int i = 0; i < size; i++) {
            cout << B[i] << endl;
        }
        cout << endl << "Точность:" << eps <<  endl;
        cout << "Начальное приближение: " << initialApproximation << endl;

        for (int i = 0; i < size; i++) {
            X[i] = initialApproximation;
        }

        try {
            cout << endl << "Решение методом Якоби:" << endl;
            if (JacobiMethod(size, eps, A_Matrix, B, X, numberOfIteration) == 1) {
                throw string("Превышен максимум итераций метода Якоби (10000), нет схождения");
            }
            for (int i = 0; i < size; i++) {
                cout << X[i] << endl;
            }
            cout << "Число итераций: " << numberOfIteration;
        }
        catch (string err_message) {
            cerr << err_message << endl;
        }

        numberOfIteration = 0;
        for (int i = 0; i < size; i++) {
            X[i] = initialApproximation;
        }

        try {
            cout << endl << "Решение методом Зейделя:" << endl;
            if (SeidelMethod(size, eps, A_Matrix, B, X, numberOfIteration) == 1) {
                throw string("Превышен максимум итераций метода Якоби (10000), нет схождения");
            }
            for (int i = 0; i < size; i++) {
                cout << X[i] << endl;
            }
            cout << "Число итераций: " << numberOfIteration;
        }
        catch (string err_message) {
            cerr << err_message << endl;
        }
    }

    _getch();
}