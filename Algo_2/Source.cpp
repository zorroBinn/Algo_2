#include <iostream>
#include <fstream>
#include <conio.h>
#include <Windows.h>
#include <cmath>
#include <float.h>
#include <iomanip>
using namespace std;
#define MAX_ITERATION 10000 //������������ ���-�� �������� ��� �������
#define RAND_MAX_ABS 10 //������������ �������� (�� ������) ������������� �����

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
    double* discrepancies = new double[size]; //������ �������
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

        fin >> eps; //������ ����������� 

        if (eps <= 0) { //���� �� ����� ��� �� ������������� - ����� ����������
            throw string("�������� � ������� � �����");
        }

        fin >> initialApproximation; //������ ���������� �����������

        ArraySelecting(size, A_Matrix, B, Xj, Xs);

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

int ZeroOnDiagonal(int size, double**& A_Matrix, double*& B) {
    for (int i = 0; i < size; i++) {
        if (fabs(A_Matrix[i][i]) < DBL_EPSILON) return 1;
    }
}

void WritingToFile(int size, double eps, double** A_Matrix, double* B, double* Xj, double* Xs, double initialApproximation, int& numberOfIteration) {
    try {
        ofstream fout("output.txt");
        if (!fout.is_open()) { // ���������, ������� �� ������� ����
            throw string("������ ��� �������� ����� output.txt, ���������� �� ����� ���� ��������");
        }

        fout << "������� �������������: " << endl;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                fout << A_Matrix[i][j] << "\t";
            }
            fout << endl;
        }

        fout << endl << "������� ��������� ������:" << endl;
        for (int i = 0; i < size; i++) {
            fout << B[i] << endl;
        }

        fout << endl << "������������ ������������� ������������:" << endl;
        for (int i = 0; i < size; i++) {
            fout << setprecision(3) << DiagonalPredominance(size, A_Matrix, i) << " ";
        }
        fout << setprecision(5);

        fout << endl << "�������� ����������� e: " << eps << endl;

        fout << "��������� �����������: " << initialApproximation << endl;

        for (int i = 0; i < size; i++) {
            Xj[i] = initialApproximation;
        }

        fout << endl << "������� ������� �����:" << endl;
        if (JacobiMethod(size, eps, A_Matrix, B, Xj, numberOfIteration) == 1) {
            fout << "�������� �������� �������� ������ ����� (10000), ��� ���������";
        }
        else {
            for (int i = 0; i < size; i++) {
                fout << Xj[i] << endl;
            }
            fout << "����� �������� ������ �����:: " << numberOfIteration << endl;
            fout << "����� �������: " << DiscrepanciesNorm(size, A_Matrix, B, Xj);
        }

        numberOfIteration = 0;
        for (int i = 0; i < size; i++) {
            Xs[i] = initialApproximation;
        }

        fout << endl << "������� ������� �������:" << endl;
        if (SeidelMethod(size, eps, A_Matrix, B, Xs, numberOfIteration) == 1) {
            fout << "�������� �������� �������� ������ ������� (10000), ��� ���������";
        }
        else {
            for (int i = 0; i < size; i++) {
                fout << Xs[i] << endl;
            }
            fout << "����� ��������: " << numberOfIteration << endl;
            fout << "����� �������: " << DiscrepanciesNorm(size, A_Matrix, B, Xs);
        }

        fout.close();
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
    double* B; //������ ��� ������ � �������� ��������� ������
    double* Xj; //V����� ��� ������� ������� ������� �����
    double* Xs; //V����� ��� ������� ������� ������� �������
    int size; //���������� - ����������� ������� �������������
    double eps; //�����������
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
            cout << "������� ����������� e: ";
            cin >> eps;
        } while (eps <= 0);
        cout << "������� ��������� �����������: ";
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