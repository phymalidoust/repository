#pragma once

#include <stdio.h>
#include <complex>
#include <vector>

class Matrix {
private:
	std::vector<std::vector<std::complex<double>>> M;
	unsigned int nCol, nRow;
public:
	Matrix(unsigned int row, unsigned int col);
	Matrix(const Matrix& other);
	Matrix(std::vector<std::vector<std::complex<double>>> m);

	std::vector<std::vector<std::complex<double>>> getMatrix() const;
	unsigned int getRows() const;
	unsigned int getColumns() const;

	void setValue(unsigned int row, unsigned int col, std::complex<double> value);
	std::complex<double> getValue(unsigned int row, unsigned int col) const;

// performs the matrix operators and checks
	Matrix& operator = (Matrix other);
	bool operator == (Matrix other) const;
	Matrix operator - (Matrix other) const;
	Matrix operator * (Matrix other) const;
	Matrix operator * (double other) const;
	Matrix operator * (std::complex<double> other) const;

	Matrix inverse() const;
	std::complex<double> trace() const;

	Matrix realPart() const;
	Matrix imagPart() const;

	~Matrix();
};

std::ostream& operator << (std::ostream& os, const Matrix& m);