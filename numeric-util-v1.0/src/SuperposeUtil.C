/*
FILE:     SuperposeUtil.C
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

#include "SuperposeUtil.h"

static void superpose_util(double c[3][3], double *r, double *xm, double *ym);

void getSuperposeMatrix(const int& lt, const std::vector<COORD>& seqf1, const std::vector<COORD>& seqf2, double *r, double *xm, double *ym)
{
       double c[3][3];
       for (int i = 0; i < 3; i++) {
            xm[i] = 0.0;
            ym[i] = 0.0;
            for (int j = 0; j < 3; j++) c[i][j] = 0.0;
       }

       for (int i = 0; i < lt; i++) {
            xm[0] += seqf1[i].x / (double) lt;
            xm[1] += seqf1[i].y / (double) lt;
            xm[2] += seqf1[i].z / (double) lt;
            ym[0] += seqf2[i].x / (double) lt;
            ym[1] += seqf2[i].y / (double) lt;
            ym[2] += seqf2[i].z / (double) lt;
       }

       double xn[3] ,yn[3];
       for (int i1 = 0; i1 < lt; i1++) {
            xn[0] = seqf1[i1].x - xm[0];
            xn[1] = seqf1[i1].y - xm[1];
            xn[2] = seqf1[i1].z - xm[2];
            yn[0] = seqf2[i1].x - ym[0];
            yn[1] = seqf2[i1].y - ym[1];
            yn[2] = seqf2[i1].z - ym[2];
            for (int i = 0; i < 3; i++) {
                 for (int j = 0; j < 3; j++) {
                      c[i][j] = xn[i] * yn[j] + c[i][j];
                 }
            }
       }
       superpose_util(c, r, xm, ym);
}

void getSuperposeMatrix(const int& lt, const std::vector<double>& seqf1, const std::vector<double>& seqf2, double *r, double *xm, double *ym)
{
       double c[3][3];
       for (int i = 0; i < 3; i++) {
            xm[i] = 0.0;
            ym[i] = 0.0;
            for (int j = 0; j < 3; j++) c[i][j] = 0.0;
       }

       for (int i = 0; i < lt; i++) {
            for (int l = 0; l < 3; l++) {
                 xm[l] += seqf1[3 * i + l] / (double) lt;
                 ym[l] += seqf2[3 * i + l] / (double) lt;
            }
       }

       double xn[3] ,yn[3];
       for (int i1 = 0; i1 < lt; i1++) {
            for (int j = 0; j < 3; j++) {
                 xn[j] = seqf1[3 * i1 + j] - xm[j];
                 yn[j] = seqf2[3 * i1 + j] - ym[j];
            }
            for (int i = 0; i < 3; i++) {
                 for (int j = 0; j < 3; j++) {
                      c[i][j] = xn[i] * yn[j] + c[i][j];
                 }
            }
       }
       superpose_util(c, r, xm, ym);
}

void getSuperposeMatrix(const int& lt, const double* seqf1, const double* seqf2, double *r, double *xm, double *ym)
{
       double c[3][3];
       for (int i = 0; i < 3; i++) {
            xm[i] = 0.0;
            ym[i] = 0.0;
            for (int j = 0; j < 3; j++) c[i][j] = 0.0;
       }

       for (int i = 0; i < lt; i++) {
            for (int l = 0; l < 3; l++) {
                 xm[l] += seqf1[3 * i + l] / (double) lt;
                 ym[l] += seqf2[3 * i + l] / (double) lt;
            }
       }

       double xn[3] ,yn[3];
       for (int i1 = 0; i1 < lt; i1++) {
            for (int j = 0; j < 3; j++) {
                 xn[j] = seqf1[3 * i1 + j] - xm[j];
                 yn[j] = seqf2[3 * i1 + j] - ym[j];
            }
            for (int i = 0; i < 3; i++) {
                 for (int j = 0; j < 3; j++) {
                      c[i][j] = xn[i] * yn[j] + c[i][j];
                 }
            }
       }
       superpose_util(c, r, xm, ym);
}

static void superpose_util(double c[3][3], double *r, double *xm, double *ym)
{
       double a[4][4], s[4][4];
       a[0][0] = -c[0][0] + c[1][1] - c[2][2];
       a[0][1] = -c[0][1] - c[1][0];
       a[0][2] = -c[1][2] - c[2][1];
       a[0][3] =  c[0][2] - c[2][0];
       a[1][1] =  c[0][0] - c[1][1] - c[2][2];
       a[1][2] =  c[0][2] + c[2][0];
       a[1][3] =  c[1][2] - c[2][1];
       a[2][2] = -c[0][0] - c[1][1] + c[2][2];
       a[2][3] =  c[0][1] - c[1][0];
       a[3][3] =  c[0][0] + c[1][1] + c[2][2];

       for (int i = 1; i < 4; i++) {
            for (int j = 0; j < i; j++) {
                 a[i][j] = a[j][i];
            }
       }
       for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                  s[i][j] = ((j == i) ? 1.0 : 0.0);
            }
       }
        

       double s0 = 0.0;
       for (int i = 1; i < 4; i++) {
            for (int j = 0; j < i; j++) {
                 s0 += 2.0 * a[i][j] * a[i][j];
            }
       }

       double s1 = sqrt(s0);
       double s2 = XEPS / 4.0 * s1;
       double s3 = s1;

       int l = 0;
       do {
            s3 = s3 / 4.0;
            do {
                 l = 0;
                 for (int iq = 1; iq < 4; iq++) {
                      for (int ip = 0; ip < iq; ip++) {
                           if (fabs(a[ip][iq]) < s3) continue;
                           l = 1;
                           double v1 = a[ip][ip];
                           double v2 = a[ip][iq];
                           double v3 = a[iq][iq];
                           double u = 0.5 * (v1 - v3);
                           double g = 1.0;
                           if (fabs(u) >= XEPS) g = -v2 / sqrt(v2 * v2 + u * u) * u / fabs(u);
                           double st = g / sqrt(2.0*(1.0+sqrt(1.0-g*g)));
                           double ct = sqrt(1.0 - st * st);
                           for (int i = 0; i < 4; i++) {
                                g = a[i][ip] * ct - a[i][iq] * st;
                                a[i][iq] = a[i][ip] * st + a[i][iq] * ct;
                                a[i][ip] = g;
                                g = s[i][ip] * ct - s[i][iq] * st;
                                s[i][iq] = s[i][ip] * st + s[i][iq] * ct;
                                s[i][ip] = g;
                           }
                           for (int i = 0; i < 4; i++) {
                                a[ip][i] = a[i][ip];
                                a[iq][i] = a[i][iq];
                           }
                           g = 2.0 * v2 * st * ct;
                           a[ip][ip] = v1 * ct * ct + v3 * st * st - g;
                           a[iq][iq] = v1 * st * st + v3 * ct * ct + g;
                           a[ip][iq] = (v1 - v3) * st * ct + v2 * (ct * ct - st * st);
                           a[iq][ip] = a[ip][iq];
                      }
                 }
            } while (l == 1);
       } while (s3 > s2);

       l = 0;
       double evmax = a[0][0];
       for (int i = 1; i < 4; i++) {
            if (a[i][i] > evmax) {
                 l = i;
                 evmax = a[i][i];
            }
       }

       r[0] = -s[0][l] * s[0][l] + s[1][l] * s[1][l] - s[2][l] * s[2][l] + s[3][l] * s[3][l];
       r[1] = 2.0 * (s[2][l] * s[3][l] - s[0][l] * s[1][l]);
       r[2] = 2.0 * (s[1][l] * s[2][l] + s[0][l] * s[3][l]);
       r[3] = -2.0 * (s[0][l] * s[1][l] + s[2][l] * s[3][l]);
       r[4] = s[0][l] * s[0][l] - s[1][l] * s[1][l] - s[2][l] * s[2][l] + s[3][l] * s[3][l];
       r[5] = 2.0 * (s[1][l] * s[3][l] - s[0][l] * s[2][l]);
       r[6] = 2.0 * (s[1][l] * s[2][l] - s[0][l] * s[3][l]);
       r[7] = -2.0 * (s[0][l] * s[2][l] + s[1][l] * s[3][l]);
       r[8] = -s[0][l] * s[0][l] - s[1][l] * s[1][l] + s[2][l] * s[2][l] + s[3][l] * s[3][l];

       for (int i = 0; i < 3; i++) {
            s0 = 0.0;
            for (int j = 0; j < 3; j++) {
                 s0 = r[j * 3 + i] * xm[j] + s0;
            }
            r[9 + i] = ym[i] - s0;
       }
}
