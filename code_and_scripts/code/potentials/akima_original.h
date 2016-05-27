// CodeCogs Commercial License Agreement
// Copyright (C) 2004-2010 CodeCogs, Zyba Ltd, Broadwood, Holford, TA5 1DU, England.
//
// This software is licensed to Thomas Truskett 
// for commercial usage by version 1.2.1 of the CodeCogs Commercial Licence. You must 
// read this License (available at www.codecogs.com) before using this software.
//
// If you distribute this file it is YOUR responsibility to ensure that all 
// recipients have a valid number of commercial licenses. You must retain a
// copy of this licence in all copies you make. 
//
// This program is distributed WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the CodeCogs Commercial Licence for more details.
//---------------------------------------------------------------------------------
//! Interpolates a given set of points using Akima spline fitting.

#ifndef MATHS_INTERPOLATION_AKIMA_H
#define MATHS_INTERPOLATION_AKIMA_H

#define ABS(x) ((x) < 0) ? -(x) : (x)

namespace Maths
{

namespace Interpolation
{

//! Interpolates a given set of points using Akima spline fitting.

class Akima
{
public:

//! Class constructor

Akima(int n, double *x, double *y)
{
            m_n = n;
            m_x = new double[n];
            m_y = new double[n];
            m_z = new double[n];
            m_t = new double[n + 3];

            for (int i = 0; i < n; ++i) {
                m_x[i] = x[i];
                m_y[i] = y[i];
            }

            for (int i = 0; i < n - 1; ++i)
                m_t[i + 2] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);

            m_t[n + 1] = 2.0 * m_t[n] - m_t[n - 1];
            m_t[n + 2] = 2.0 * m_t[n + 1] - m_t[n];
            m_t[1] = 2.0 * m_t[2] - m_t[3];
            m_t[0] = 2.0 * m_t[1] - m_t[2];

            for (int i = 0; i < n; ++i)  {

                double a = ABS(m_t[i + 3] - m_t[i + 2]);
                double b = ABS(m_t[i + 1] - m_t[i]);

                if (a + b) m_z[i] = (a * m_t[i + 1] + b * m_t[i + 2]) / (a + b);
                else m_z[i] = (m_t[i + 2] + m_t[i + 1]) / 2.0;
            }
        }

//! Class destructor

~Akima()
{
            delete [] m_x;
            delete [] m_y;
            delete [] m_z;
            delete [] m_t;
        }

//! Returns an interpolated value.

double getValue(double x)
{
      if ((x <= m_x[0]) || (m_n < 3))
        return (((m_y[1] - m_y[0]) / (m_x[1] - m_x[0]))*(x - m_x[0]) + m_y[0]);      
      
      // Problems with range excess
      // Variant A
     if (x >= m_x[m_n - 1])
        return (((m_y[m_n - 1] - m_y[m_n - 2]) / (m_x[m_n - 1] - m_x[m_n - 2]))*(x - m_x[m_n - 2]) + m_y[m_n - 2]);

      // Variant B
      int i = 0;
            while ((i < (m_n - 2)) && (x > m_x[i + 1])) ++i;

            double a = x - m_x[i];
            double b = m_x[i + 1] - m_x[i];
            return m_y[i] + m_z[i] * a
            + (3.0 * m_t[i + 2] - 2.0 * m_z[i] - m_z[i + 1]) * a * a / b
            + (m_z[i] + m_z[i + 1] - 2.0 * m_t[i + 2]) * a * a * a / (b * b);
        }

private:

int m_n;

double *m_x, *m_y, *m_z, *m_t;
};


//! A static function implementing the Akima Class for one off calculations

double Akima_once(int N, double *x, double *y, double a )
{
  // This function is created to enable an Instant Calculator on CodeCogs. 
  // You probably shouldn't be using this function otherwise. 

   Maths::Interpolation:: Akima A(N, x, y);
   return A.getValue(a);
}

}

}

#undef ABS

#endif

