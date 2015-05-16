/* Ellipse manipulations
    Copyright (C) 2014 Victoria Rudakova <vicrucann@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef ELLIPSE_OPERATIONS_H
#define ELLIPSE_OPERATIONS_H
#include "matrix.h"
#include <vector>

using namespace libNumerics;

template <typename T>
matrix<T> getCircleMatrix(T x, T y, T r) {
	matrix<T> circle = matrix<T>::zeros(3,3);
	circle(0,0) = 1;
	circle(1,1) = 1;
	circle(2,0) = circle(0,2) = -x;
	circle(1,2) = circle(2,1) = -y;
	circle(2,2) = -r*r+x*x+y*y;
	return circle; }

template <typename T>
void getEllipseCenter(matrix<T>& ell, T &x, T &y) 
{
	T A = ell(0,0), B = ell(0,1), D = ell(0,2), C = ell(1,1), E = ell(1,2);
	x = (B*E - C*D) / (A*C - B*B);
	y = (B*D - A*E) / (A*C - B*B);
	/*matrix<T> s22 = ell.copy(0, 1, 0, 1);
	matrix<T> s3 = ell.copy(0, 1, 2, 2);
	matrix<T> CS = -s22.inv()*s3;
	x = CS(0,0);
	y = CS(1,0); */
}

template <typename T>
matrix<T> rotateCircle(const matrix<T>& S, const matrix<T>& H) {
	return ( (H.inv()).t() ) * S * H.inv();
}

template <typename T>
void rotateVirtualImg(std::vector<matrix<T> >& S, matrix<T>& H0, matrix<T>& CH0S, matrix<T>& H0CS)
{
	int nellip = S.size(); 
	for (int i = 0; i < nellip; i++)
	{
		matrix<T> s = S[i];
		T csx, csy;
		getEllipseCenter(s, csx, csy);
		matrix<T> cs = matrix<T>::zeros(3,1); cs(0,0) = csx; cs(1,0) = csy; cs(2,0) = 1;
		matrix<T> h0cs = H0*cs;
		H0CS(0,i) = h0cs(0,0)/h0cs(2,0); H0CS(1,i) = h0cs(1,0)/h0cs(2,0);
		
		matrix<T> h0s = rotateCircle(s,H0);
		T ch0sx, ch0sy;
		getEllipseCenter(h0s, ch0sx, ch0sy);
		CH0S(0,i) = ch0sx; CH0S(1,i) = ch0sy;
	}
}

#endif
