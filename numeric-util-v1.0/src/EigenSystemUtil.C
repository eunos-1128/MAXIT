/*
FILE:     EigenSystemUtil.C
*/
/*
VERSION:  11.200
*/
/*
DATE:     4/20/2024
*/
/*
  Comments and Questions to: sw-help@rcsb.rutgers.edu
*/
/*
COPYRIGHT 1999-2024 Rutgers - The State University of New Jersey

This software is provided WITHOUT WARRANTY OF MERCHANTABILITY OR
FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR
IMPLIED.  RUTGERS MAKE NO REPRESENTATION OR WARRANTY THAT THE
SOFTWARE WILL NOT INFRINGE ANY PATENT, COPYRIGHT OR OTHER
PROPRIETARY RIGHT.

The user of this software shall indemnify, hold harmless and defend
Rutgers, its governors, trustees, officers, employees, students,
agents and the authors against any and all claims, suits,
losses, liabilities, damages, costs, fees, and expenses including
reasonable attorneys' fees resulting from or arising out of the
use of this software.  This indemnification shall include, but is
not limited to, any and all claims alleging products liability.
*/
/*
               RCSB PDB SOFTWARE LICENSE AGREEMENT

BY CLICKING THE ACCEPTANCE BUTTON OR INSTALLING OR USING 
THIS "SOFTWARE, THE INDIVIDUAL OR ENTITY LICENSING THE  
SOFTWARE ("LICENSEE") IS CONSENTING TO BE BOUND BY AND IS 
BECOMING A PARTY TO THIS AGREEMENT.  IF LICENSEE DOES NOT 
AGREE TO ALL OF THE TERMS OF THIS AGREEMENT
THE LICENSEE MUST NOT INSTALL OR USE THE SOFTWARE.

1. LICENSE AGREEMENT

This is a license between you ("Licensee") and the Protein Data Bank (PDB) 
at Rutgers, The State University of New Jersey (hereafter referred to 
as "RUTGERS").   The software is owned by RUTGERS and protected by 
copyright laws, and some elements are protected by laws governing 
trademarks, trade dress and trade secrets, and may be protected by 
patent laws. 

2. LICENSE GRANT

RUTGERS grants you, and you hereby accept, non-exclusive, royalty-free 
perpetual license to install, use, modify, prepare derivative works, 
incorporate into other computer software, and distribute in binary 
and source code format, or any derivative work thereof, together with 
any associated media, printed materials, and on-line or electronic 
documentation (if any) provided by RUTGERS (collectively, the "SOFTWARE"), 
subject to the following terms and conditions: (i) any distribution 
of the SOFTWARE shall bind the receiver to the terms and conditions 
of this Agreement; (ii) any distribution of the SOFTWARE in modified 
form shall clearly state that the SOFTWARE has been modified from 
the version originally obtained from RUTGERS.  

2. COPYRIGHT; RETENTION OF RIGHTS.  

The above license grant is conditioned on the following: (i) you must 
reproduce all copyright notices and other proprietary notices on any 
copies of the SOFTWARE and you must not remove such notices; (ii) in 
the event you compile the SOFTWARE, you will include the copyright 
notice with the binary in such a manner as to allow it to be easily 
viewable; (iii) if you incorporate the SOFTWARE into other code, you 
must provide notice that the code contains the SOFTWARE and include 
a copy of the copyright notices and other proprietary notices.  All 
copies of the SOFTWARE shall be subject to the terms of this Agreement.  

3. NO MAINTENANCE OR SUPPORT; TREATMENT OF ENHANCEMENTS 

RUTGERS is under no obligation whatsoever to: (i) provide maintenance 
or support for the SOFTWARE; or (ii) to notify you of bug fixes, patches, 
or upgrades to the features, functionality or performance of the 
SOFTWARE ("Enhancements") (if any), whether developed by RUTGERS 
or third parties.  If, in its sole discretion, RUTGERS makes an 
Enhancement available to you and RUTGERS does not separately enter 
into a written license agreement with you relating to such bug fix, 
patch or upgrade, then it shall be deemed incorporated into the SOFTWARE 
and subject to this Agreement. You are under no obligation whatsoever 
to provide any Enhancements to RUTGERS or the public that you may 
develop over time; however, if you choose to provide your Enhancements 
to RUTGERS, or if you choose to otherwise publish or distribute your 
Enhancements, in source code form without contemporaneously requiring 
end users or RUTGERS to enter into a separate written license agreement 
for such Enhancements, then you hereby grant RUTGERS a non-exclusive,
royalty-free perpetual license to install, use, modify, prepare
derivative works, incorporate into the SOFTWARE or other computer
software, distribute, and sublicense your Enhancements or derivative
works thereof, in binary and source code form.

4. FEES.  There is no license fee for the SOFTWARE.  If Licensee
wishes to receive the SOFTWARE on media, there may be a small charge
for the media and for shipping and handling.  Licensee is
responsible for any and all taxes.

5. TERMINATION.  Without prejudice to any other rights, Licensor
may terminate this Agreement if Licensee breaches any of its terms
and conditions.  Upon termination, Licensee shall destroy all
copies of the SOFTWARE.

6. PROPRIETARY RIGHTS.  Title, ownership rights, and intellectual
property rights in the Product shall remain with RUTGERS.  Licensee 
acknowledges such ownership and intellectual property rights and will 
not take any action to jeopardize, limit or interfere in any manner 
with RUTGERS' ownership of or rights with respect to the SOFTWARE.  
The SOFTWARE is protected by copyright and other intellectual 
property laws and by international treaties.  Title and related 
rights in the content accessed through the SOFTWARE is the property 
of the applicable content owner and is protected by applicable law.  
The license granted under this Agreement gives Licensee no rights to such
content.

7. DISCLAIMER OF WARRANTY.  THE SOFTWARE IS PROVIDED FREE OF 
CHARGE, AND, THEREFORE, ON AN "AS IS" BASIS, WITHOUT WARRANTY OF 
ANY KIND, INCLUDING WITHOUT LIMITATION THE WARRANTIES THAT IT 
IS FREE OF DEFECTS, MERCHANTABLE, FIT FOR A PARTICULAR PURPOSE 
OR NON-INFRINGING.  THE ENTIRE RISK AS TO THE QUALITY AND 
PERFORMANCE OF THE SOFTWARE IS BORNE BY LICENSEE.  SHOULD THE 
SOFTWARE PROVE DEFECTIVE IN ANY RESPECT, THE LICENSEE AND NOT 
LICENSOR ASSUMES THE ENTIRE COST OF ANY SERVICE AND REPAIR.  
THIS DISCLAIMER OF WARRANTY CONSTITUTES AN ESSENTIAL PART OF 
THIS AGREEMENT.  NO USE OF THE PRODUCT IS AUTHORIZED HEREUNDER 
EXCEPT UNDER THIS DISCLAIMER.

8. LIMITATION OF LIABILITY.  TO THE MAXIMUM EXTENT PERMITTED BY
APPLICABLE LAW,  IN NO EVENT WILL LICENSOR BE LIABLE FOR ANY 
INDIRECT, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING 
OUT OF THE USE OF OR INABILITY TO USE THE SOFTWARE, INCLUDING, 
WITHOUT LIMITATION, DAMAGES FOR LOSS OF GOODWILL, WORK 
STOPPAGE, COMPUTER FAILURE OR MALFUNCTION, OR ANY AND ALL 
OTHER COMMERCIAL DAMAGES OR LOSSES, EVEN IF ADVISED OF THE
POSSIBILITY THEREOF. 
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <string>

#include "EigenSystemUtil.h"

static void jacobi(std::vector<double>& a, const int& n, std::vector<double>& d, std::vector<double>& v);
static void eigsrt(const int& n, std::vector<double>& d, std::vector<double>& v);
static void rotate(std::vector<double>& a, const int& n, const int& i, const int& j, const int& k, const int& l, double& g, double& h,
                   const double& s, const double& tau);
static void print_vector(const std::string& name, const int& dimension, const std::vector<double>& vec);

void eigen_system(const std::vector<COORD>& coords)
{
       std::list<std::vector<double> > double_coords;
       double_coords.clear();

       std::vector<double> xyz;
       for (std::vector<COORD>::const_iterator cpos = coords.begin(); cpos != coords.end(); ++cpos) {
            xyz.clear();
            xyz.push_back(cpos->x);
            xyz.push_back(cpos->y);
            xyz.push_back(cpos->z);
            double_coords.push_back(xyz);
       }
       eigen_system(3, double_coords);
}

void eigen_system(const std::list<COORD>& coords)
{
       std::list<std::vector<double> > double_coords;
       double_coords.clear();

       std::vector<double> xyz;
       for (std::list<COORD>::const_iterator cpos = coords.begin(); cpos != coords.end(); ++cpos) {
            xyz.clear();
            xyz.push_back(cpos->x);
            xyz.push_back(cpos->y);
            xyz.push_back(cpos->z);
            double_coords.push_back(xyz);
       }
       eigen_system(3, double_coords);
}

void eigen_system(const int& dimension, const std::vector<std::vector<double> >& coords)
{
       std::list<std::vector<double> > coord_list;
       coord_list.clear();
       for (std::vector<std::vector<double> >::const_iterator vpos = coords.begin(); vpos != coords.end(); ++vpos) {
            coord_list.push_back(*vpos);
       }
       eigen_system(dimension, coord_list);
}

void eigen_system(const int& dimension, const std::list<std::vector<double> >& coords)
{
       std::vector<double> Vector, Matrix, eigenValue, eigenVector;
       for (int i = 0; i < dimension; ++i) Vector.push_back(0);
       for (int i = 0; i < dimension * dimension; ++i) Matrix.push_back(0);

       for (std::list<std::vector<double> >::const_iterator vpos = coords.begin(); vpos != coords.end(); ++vpos) {
            for (int i = 0; i < dimension; ++i) {
                 Vector[i] += (*vpos)[i];
                 for (int j = i; j < dimension; ++j) {
                      Matrix[i * dimension + j] += (*vpos)[i] * (*vpos)[j];
                 }
            }
       }

       for (int i = 0; i < dimension; ++i) Vector[i] /= (double) coords.size();
       for (int i = 0; i < dimension; ++i) {
            for (int j = i; j < dimension; ++j) {
                 Matrix[i * dimension + j] /= (double) coords.size();
                 Matrix[i * dimension + j] -= Vector[i] * Vector[j];
                 if (j != i) Matrix[j * dimension + i] = Matrix[i * dimension + j];
            }
       }

       print_vector("Vector", dimension, Vector);
       print_vector("Matrix", dimension, Matrix);

       jacobi(Matrix, dimension, eigenValue, eigenVector);

       print_vector("eigenValue", dimension, eigenValue);
       print_vector("Matrix", dimension, Matrix);
       print_vector("eigenVector", dimension, eigenVector);
}

static void jacobi(std::vector<double>& a, const int& n, std::vector<double>& d, std::vector<double>& v)
{
       std::vector<double> b, z;
       d.clear();
       b.clear();
       z.clear();
       for (int ip = 0; ip < n; ++ip) {
            d.push_back(a[ip * n + ip]);
            b.push_back(a[ip * n + ip]);
            z.push_back(0);
       }

       v.clear();
       for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                 if (j == i) v.push_back(1);
                 else v.push_back(0);
            }
       }

       double t;
       for (int i = 0; i < 100; ++i) {
            double sm = 0.0;
            for (int ip = 0; ip < (n - 1); ++ip) {
                 for (int iq = ip + 1; iq < n; ++iq) sm += fabs(a[ip * n + iq]);
            }
            if (sm < XEPS) {
                 eigsrt(n, d, v);
                 return;
            }

            double tresh = 0.0;
            if (i < 4) tresh = 0.2 * sm / (double (n * n));

            for (int ip = 0; ip < (n - 1); ++ip) {
                 for (int iq = ip + 1; iq < n; ++iq) {
                      double g = 100.0 * fabs(a[ip * n + iq]);
                      if (i > 4 && (fabs(d[ip]) + g) == fabs(d[ip]) && (fabs(d[iq]) + g) == fabs(d[iq])) a[ip * n + iq] = 0.0;
                      else if (fabs(a[ip * n + iq]) > tresh) {
                           double h = d[iq] - d[ip];
                           if ((fabs(h) + g) == fabs(h)) t = a[ip * n + iq] / h;
                           else {
                                double theta = 0.5 * h / a[ip * n + iq];
                                t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                                if (theta < 0.0) t = -t;
                           }
                           double c = 1.0 / sqrt(1 + t * t);
                           double s = t * c;
                           double tau = s / (1.0 + c);
                           h = t * a[ip * n + iq];
                           z[ip] -= h;
                           z[iq] += h;
                           d[ip] -= h;
                           d[iq] += h;
                           a[ip * n + iq] = 0.0;
                           for (int j = 0; j < ip; j++) rotate(a, n, j, ip, j, iq, g, h, s, tau);
                           for (int j = ip + 1; j < iq; j++) rotate(a, n, ip, j, j, iq, g, h, s, tau);
                           for (int j = iq + 1; j < n; j++) rotate(a, n, ip, j, iq, j, g, h, s, tau);
                           for (int j = 0; j < n; j++) rotate(v, n, j, ip, j, iq, g, h, s, tau);
                      }
                 }
            }

            for (int ip = 0; ip < n; ++ip) {
                 b[ip] += z[ip];
                 d[ip] = b[ip];
                 z[ip] = 0.0;
            }
       }
       d.clear();
       v.clear();
}

static void eigsrt(const int& n, std::vector<double>& d, std::vector<double>& v)
{
       for (int i = 0; i < n; ++i) {
            double p = d[i];
            int k = i;
            for (int j = i + 1; j < n; ++j) {
                 if (d[j] < p) {
                      p = d[j];
                      k = j;
                 }
            }
            if (k != i) {
                 d[k] = d[i];
                 d[i] = p;
                 for (int j = 0; j < n; ++j) {
                      double temp = v[j * n + i];
                      v[j * n + i] = v[j * n + k];
                      v[j * n + k] = temp;
                 }
            }
       }
}

static void rotate(std::vector<double>& a, const int& n, const int& i, const int& j, const int& k, const int& l, double& g, double& h,
                   const double& s, const double& tau)
{
       g = a[i * n + j];
       h = a[k * n + l];
       a[i * n + j] = g - s * (h + g * tau);
       a[k * n + l] = h + s * (g - h * tau);
}

static void print_vector(const std::string& name, const int& dimension, const std::vector<double>& vec)
{
       printf("Print %s:", name.c_str());
       for (unsigned int i = 0; i < vec.size(); ++i) {
            if ((i % dimension) == 0) printf("\n");
            printf("%.10f ", vec[i]);
       }
       printf("\n\n");
}
