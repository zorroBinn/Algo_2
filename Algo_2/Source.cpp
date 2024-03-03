#include <iostream>
#include <fstream>
#include <conio.h>
#include <Windows.h>
#include <cmath>
#include <float.h>
using namespace std;
#define MAX_ITERATION 10000 //������������ ���-�� �������� ��� �������
#define RAND_MAX_ABS 10 //������������ �������� (�� ������) ������������� �����

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
    char temp; //���������� ��� �������� �����
    //������ �� �����
    try {
        ifstream fin("input.txt");
        if (!fin.is_open()) { // ���������, ������� �� ������� ����
            throw string("�� ������� ������� ���� input.txt");
        }

        fin >> size; //������ ����������� ������� �������������

        if (size < 2) { //���� �� ����� ��� �� ������������� - ����� ����������
            throw string("�������� � ������� � �����");
        }

        fin >> eps; //������ ����������� ������� 

        if (eps <= 0) { //���� �� ����� ��� �� ������������� - ����� ����������
            throw string("�������� � ������� � �����");
        }

        fin >> initialApproximation; //������ ���������� �����������

        ArraySelecting(size, A_Matrix, B, X);

        //������ ������� �������������
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (!(fin >> A_Matrix[i][j]) || fin.eof()) { //���� �� ����� ��� �������� ����� ����� - ����� ����������
                    throw string("�������� � ������� � �����");
                }
            }
        }

        //������ ������� ��������� ������
        for (int i = 0; i < size; i++) {
            if (!(fin >> B[i])) {
                throw string("�������� � ������� � �����. � ������ �� �����"); //���� �� ����� - ����� ����������
            }
            if (fin.eof() && i != size - 1) { //���� ��������� ����� ����� �� ���������� ��������� - ����� ����������
                throw string("�������� � ������� � �����. ����������� ������� ������ ��������������");
            }
            if (i == size - 1) { //���� � ����� ������ ������ ��� ����������� ������� ������� ������ �������������� - ����� ����������
                if (fin >> temp)
                    throw string("� ����� ������ ������ ��� ����������� ������� ������� ������ ��������������");
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
        if (!fout.is_open()) { // ���������, ������� �� ������� ����
            throw string("������ ��� �������� ����� output.txt, ���������� �� ����� ���� ��������");
        }

        fout << "�������� ������� �������������: " << endl;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                fout << A_Matrix[i][j];
            }
            fout << endl;
        }

        fout << endl << "�������� ������� ��������� ������:" << endl;
        for (int i = 0; i < size; i++) {
            fout << B[i] << endl;
        }

        fout << endl << "�������� �������� e = " << eps << endl;

        cout << endl << "��������� ������� ���������. ���������� �������� � ���� output.txt" << endl;
    }
    catch (string err_message) {
        cerr << err_message << endl;
    }
}

int main() {
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
    double** A_Matrix; //��������� ������ ��� ������ � �������� �������������
    double* B; //��������� ������ ��� ������ � �������� ��������� ������
    double* X; //��������� ������ ��� ������� �������
    int size; //���������� - ����������� ������� �������������
    double eps; //��������
    double initialApproximation; //��������� �����������
    int numberOfIteration = 0; //���-�� ��������
    string IsRandom, method;

    do {
        system("cls");
        cout << "������������� ������������ �������� �������?" << endl;
        cout << "1 - ��" << endl << "2 - ������� ������� �� �����" << endl;
        cin >> IsRandom;
    } while (IsRandom != "1" && IsRandom != "2");
    
    if (IsRandom == "1") {
        do {
            cout << "������� ����������� �������: ";
            cin >> size;
        } while (size < 2);
        do {
            cout << "������� �������� e: ";
            cin >> eps;
        } while (eps <= 0);
        cout << "������� ��������� �����������: ";
        cin >> initialApproximation;
        
        ArraySelecting(size, A_Matrix, B, X);
        RandomMatrix(size, A_Matrix, B);

        cout << "���������������� ������� �:" << endl;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                cout << A_Matrix[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl << "��������������� ������� ��������� ������:" << endl;
        for (int i = 0; i < size; i++) {
            cout << B[i] << endl;
        }

        for (int i = 0; i < size; i++) {
            X[i] = initialApproximation;
        }

        try {
            cout << endl << "������� ������� �����:" << endl;
            if (JacobiMethod(size, eps, A_Matrix, B, X, numberOfIteration) == 1) {
                throw string("�������� �������� �������� ������ ����� (10000), ��� ���������");
            }
            for (int i = 0; i < size; i++) {
                cout << X[i] << endl;
            }
            cout << "����� ��������: " << numberOfIteration;
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
            cout << endl << "������� ������� �������:" << endl;
            if (SeidelMethod(size, eps, A_Matrix, B, X, numberOfIteration) == 1) {
                throw string("�������� �������� �������� ������ ������� (10000), ��� ���������");
            }
            for (int i = 0; i < size; i++) {
                cout << X[i] << endl;
            }
            cout << "����� ��������: " << numberOfIteration;
        }
        catch (string err_message) {
            cerr << err_message << endl;
            _getch();
            exit(1);
        }
    }
    
    if (IsRandom == "2") {
        
        ReadingFromFile(size, eps, initialApproximation, A_Matrix, B, X);
        
        cout << "������� �:" << endl;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                cout << A_Matrix[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl << "������� ��������� ������:" << endl;
        for (int i = 0; i < size; i++) {
            cout << B[i] << endl;
        }
        cout << endl << "��������:" << eps <<  endl;
        cout << "��������� �����������: " << initialApproximation << endl;

        for (int i = 0; i < size; i++) {
            X[i] = initialApproximation;
        }

        try {
            cout << endl << "������� ������� �����:" << endl;
            if (JacobiMethod(size, eps, A_Matrix, B, X, numberOfIteration) == 1) {
                throw string("�������� �������� �������� ������ ����� (10000), ��� ���������");
            }
            for (int i = 0; i < size; i++) {
                cout << X[i] << endl;
            }
            cout << "����� ��������: " << numberOfIteration;
        }
        catch (string err_message) {
            cerr << err_message << endl;
        }

        numberOfIteration = 0;
        for (int i = 0; i < size; i++) {
            X[i] = initialApproximation;
        }

        try {
            cout << endl << "������� ������� �������:" << endl;
            if (SeidelMethod(size, eps, A_Matrix, B, X, numberOfIteration) == 1) {
                throw string("�������� �������� �������� ������ ����� (10000), ��� ���������");
            }
            for (int i = 0; i < size; i++) {
                cout << X[i] << endl;
            }
            cout << "����� ��������: " << numberOfIteration;
        }
        catch (string err_message) {
            cerr << err_message << endl;
        }
    }

    _getch();
}