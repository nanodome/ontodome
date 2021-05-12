/*
    NanoDome - H2020 European Project NanoDome, GA n.646121
    (www.nanodome.eu, https://github.com/nanodome/nanodome)
    e-mail: Emanuele Ghedini, emanuele.ghedini@unibo.it

    Copyright (C) 2018  Alma Mater Studiorum - Università di Bologna

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// RKF45 is a C++ library which implements the Watt and Shampine RKF45 ODE solver.
// The RKF45 ODE solver is a Runge-Kutta-Fehlberg algorithm for solving an ordinary differential equation, with automatic error estimation using rules of order 4 and 5.

#ifndef RK45_H
#define RK45_H

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

class RK45 {

public:
  float r4_abs ( float x );
  float r4_epsilon ( );
  void r4_fehl ( void f ( float t, float y[], float yp[] ), int neqn,
  float y[], float t, float h, float yp[], float f1[], float f2[], float f3[],
  float f4[], float f5[], float s[] );
  float r4_max ( float x, float y );
  float r4_min ( float x, float y );
  int r4_rkf45 ( void f ( float t, float y[], float yp[] ), int neqn,
  float y[], float yp[], float *t, float tout, float *relerr, float abserr,
  int flag );
  float r4_sign ( float x );
  double r8_abs ( double x );
  double r8_epsilon ( );
  void r8_fehl ( void f ( double t, double y[], double yp[] ), int neqn,
  double y[], double t, double h, double yp[], double f1[], double f2[], double f3[],
  double f4[], double f5[], double s[] );
  double r8_max ( double x, double y );
  double r8_min ( double x, double y );
  int r8_rkf45 ( void f ( double t, double y[], double yp[] ), int neqn,
  double y[], double yp[], double *t, double tout, double *relerr, double abserr,
  int flag );
  double r8_sign ( double x );

  void timestamp ( );

};

#endif // RK45_H
