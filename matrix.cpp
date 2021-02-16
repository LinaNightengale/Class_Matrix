#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cmath>
#include <algorithm>
using namespace std;

template <typename T>
class Matrix
{
    private:
        int m, n;
        T * matrix;

    public:
        Matrix();
        Matrix(int m, int n);
        Matrix(Matrix& matrix);
        ~Matrix();

        const T operator[](const int i);
        Matrix operator=(const Matrix & matrix_);
        Matrix operator+(const Matrix & matrix_);
        Matrix operator-(const Matrix & matrix_);
        Matrix operator*(const Matrix & matrix_);
        Matrix operator*(const T constanta);

        Matrix & input(T matrix[]);
        void output();
        void set(int i, int j, T value);
        T get(int i, int j, int n); 

        Matrix & zero();
        Matrix & identity_matrix();  
        Matrix transpose();
        void householderDecomposition(Matrix<T> & R);
        void reverse_Gaussian(Matrix<T> & matrix_a, Matrix<T> & phi, int N);

};

template <typename T>
Matrix<T>::Matrix()
{
    n = 0, m = 0;
    matrix = nullptr;
}

template <typename T>
Matrix<T>::Matrix(int m, int n)
{
   this->m = m;
   this->n = n;
   this->matrix = new T[m * n];
   
   for (int j = 0; j < m; ++j)
		for (int i = 0; i < n; ++i)
			matrix[j * n + i] = 0;
}


template <typename T>
Matrix<T>::Matrix(Matrix & matrix)
{
    m = matrix.m;
    n = matrix.n;
    
    this->matrix = new T[m * n];
	for (int i = 0; i < m * n; ++i)
        this->matrix[i] = matrix[i];
}

template <typename T>
Matrix<T>::~Matrix()
{
    delete [] matrix;
}

template <typename T>
const T Matrix<T>::operator[](const int i) 
{
    return matrix[i];
}

template <typename T>
Matrix<T> Matrix<T>::operator=(const Matrix & matrix_)
{
    if (m * n != matrix_.m * matrix_.n) 
    {
		delete [] matrix;
		matrix = new T[matrix_.m * matrix_.n];
	}

    m = matrix_.m;
    n = matrix_.n;

    for (int i = 0; i < m * n; ++i) matrix[i] = matrix_.matrix[i];

	return * this;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix & matrix_)
{
    Matrix result(m, n);

    for (int j = 0; j < m; ++j)
        for (int i = 0; i < n; ++i)
            result.matrix[j * n + i] = matrix[j * n + i] + matrix_.matrix[j * n + i];
    
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix & matrix_)
{
    Matrix result(m, n);

    for (int j = 0; j < m; ++j)
        for (int i = 0; i < n; ++i)
            result.matrix[j * n + i] = matrix[j * n + i] - matrix_.matrix[j * n + i];

    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix & matrix_) 
{
    Matrix result(m, matrix_.n);
    
    for (int j = 0; j < m; ++j)
		for (int i = 0; i < matrix_.n; ++i) 
        {
			result.matrix[j * matrix_.n + i] = 0;
			for (int k = 0; k < n; k++)
				result.matrix[j * matrix_.n + i] += matrix[j * n + k] * matrix_.matrix[k * matrix_.n + i];
		}

	return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const T constanta) // проблема тута
{
    Matrix result(m, n);

    for (int j = 0; j < m; ++j)
		for (int i = 0; i < n; ++i)
			result.matrix[j * n + i] = matrix[j * n + i] * constanta;

	return result;
}

template <typename T>
Matrix<T> Matrix<T>::transpose()
{
    Matrix result(n, m);

    for (int j = 0; j < m; ++j)
		for (int i = 0; i < n; ++i)
			result.matrix[i * m + j] = matrix[j * n + i];

	return result;
}

template <typename T>
Matrix<T> & Matrix<T>::zero() 
{
	for (int j = 0; j < m; ++j)
		for (int i = 0; i < n; ++i)
			matrix[j * n + i] = 0.0;

	return * this;
}

template <typename T>
Matrix<T> & Matrix<T>::identity_matrix() 
{
	for (int j = 0; j < m; ++j)
		for (int i = 0; i < n; ++i)
			matrix[j * n + i] = i == j ? 1.0 : 0.0;

	return * this;
}

template <typename T>
Matrix<T> & Matrix<T>::input(T matrix[]) 
{
	for (int i = 0; i < m * n; ++i) this->matrix[i] = matrix[i];
	return * this;
}

template <typename T>
void Matrix<T>::output() 
{
	for (int j = 0; j < m; ++j) 
	{
		for (int i = 0; i < n; ++i) 
			cout << setw(16) << matrix[j * n + i] << "   ";
		cout << endl;
	}
	cout << endl;
}

template <typename T>
void Matrix<T>::householderDecomposition(Matrix<T> & R) 
{
	double mag, alpha;
	Matrix<double> u(m, 1), v(m, 1);
	Matrix<double> P(m, m), I(m, m);
    P.identity_matrix();
    I.identity_matrix();
    
	R = *this;

	for (int i = 0; i < n; i++) {
		u.zero(); 
        v.zero();
		
		mag = 0.0;
		for (int j = i; j < m; j++) {
			u.matrix[j] = R.matrix[j * n + i];
			mag += u.matrix[j] * u.matrix[j];
		}
		mag = sqrt(mag);
		
		alpha = u.matrix[i] < 0 ? mag : -mag;

		mag = 0.0;
		for (int j = i; j < m; j++) {
			v.matrix[j] = j == i ? u.matrix[j] + alpha : u.matrix[j];
			mag += v.matrix[j] * v.matrix[j];
		}
		mag = sqrt(mag);

		if (mag < 0.0000000001) continue;

		for (int j = i; j < m; j++) v.matrix[j] /= mag;

		P = I - (v * v.transpose()) * 2.0;

		R = P * R;
	}
}

template <typename T>
void Matrix<T>::reverse_Gaussian(Matrix<T> & matrix_a, Matrix<T> & phi, int N)
{
    Matrix<T> matrix_x(N, N);

    T buffer_1, buffer_2;

    for (int i = N - 1; i >= 0; --i)
    {
        buffer_1 = 0.;
        for (int j = i + 1; j < N; ++j)
        {
            buffer_2 = matrix_a.matrix[i][j] * matrix_x.matrix[j];
            buffer_1 += buffer_2;
        }
        matrix_x.matrix[i] = (phi.matrix[i] - buffer_1) / matrix_a.matrix[i][i];
    }
    return matrix_x;
}
int main()
{
    return 0;
}