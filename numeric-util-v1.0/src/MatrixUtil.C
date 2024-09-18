/*
FILE:     MatrixUtil.C
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
#include <string.h>
#include <math.h>

#include "MatrixUtil.h"

void matrix_init(double m[4][4])
{
       m[0][0] = 1.0; m[0][1] = 0.0; m[0][2] = 0.0; m[0][3] = 0.0;
       m[1][0] = 0.0; m[1][1] = 1.0; m[1][2] = 0.0; m[1][3] = 0.0;
       m[2][0] = 0.0; m[2][1] = 0.0; m[2][2] = 1.0; m[2][3] = 0.0;
       m[3][0] = 0.0; m[3][1] = 0.0; m[3][2] = 0.0; m[3][3] = 1.0;
}

void matrix_concatenate(double dest[4][4], const double src[4][4])
{
       double t[4][4];
     
       memcpy (t, dest, 16 * sizeof (double));
     
       for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                 dest[i][j] = t[i][0] * src[0][j] + t[i][1] * src[1][j] +
           	              t[i][2] * src[2][j] + t[i][3] * src[3][j];
            }
       }
}

void matrix_x_rotation(double m[4][4], const double& radians)
{
       matrix_init(m);
     
       double d = cos(radians);
       m[1][1] = d;
       m[2][2] = d;
       d = sin(radians);
       m[1][2] = d;
       m[2][1] = -d;
}

void matrix_y_rotation(double m[4][4], const double& radians)
{
       matrix_init(m);
     
       double d = cos(radians);
       m[0][0] = d;
       m[2][2] = d;
       d = sin(radians);
       m[2][0] = d;
       m[0][2] = -d;
}

void matrix_z_rotation(double m[4][4], const double& radians)
{
       matrix_init(m);
     
       double d = cos(radians);
       m[0][0] = d;
       m[1][1] = d;
       d = sin(radians);
       m[0][1] = d;
       m[1][0] = -d;
}

void matrix_translation(double m[4][4], const double& x, const double& y, const double& z)
{
       matrix_init(m);
     
       m[3][0] = x;
       m[3][1] = y;
       m[3][2] = z;
}

void matrix_scale(double m[4][4], const double& sx, const double& sy, const double& sz)
{
       matrix_init(m);
     
       m[0][0] = sx;
       m[1][1] = sy;
       m[2][2] = sz;
}

void matrix_orthographic_projection(double m[4][4], const double& left, const double& right, const double& bottom,
                                    const double& top, const double& near, const double& far)
{
       matrix_init(m);
     
       m[0][0] = 2.0 / (right - left);
       m[3][0] = - (right + left) / (right - left);
       m[1][1] = 2.0 / (top - bottom);
       m[3][1] = - (top + bottom) / (top - bottom);
       m[2][2] = -2.0 / (far - near);
       m[3][2] = - (far + near) / (far - near);
}

void matrix_transform(COORD &v, const double m[4][4])
{
       COORD old = v;
       matrix_transform(v, old, m);
}

void matrix_transform(COORD& new_coord, const COORD& old_coord, const double m[4][4])
{
       double w = old_coord.x * m[0][3] + old_coord.y * m[1][3] + old_coord .z * m[2][3] + m[3][3];
       new_coord.x = (old_coord.x * m[0][0] + old_coord.y * m[1][0] + old_coord.z * m[2][0] + m[3][0]) / w;
       new_coord.y = (old_coord.x * m[0][1] + old_coord.y * m[1][1] + old_coord.z * m[2][1] + m[3][1]) / w;
       new_coord.z = (old_coord.x * m[0][2] + old_coord.y * m[1][2] + old_coord.z * m[2][2] + m[3][2]) / w;
} 

void matrix_rotate(COORD &v, const double m[4][4])
{
       COORD old = v;
       matrix_rotate(v, old, m);
}

void matrix_rotate(COORD &new_coord, const COORD& old_coord, const double m[4][4])
{
       new_coord.x = old_coord.x * m[0][0] + old_coord.y * m[0][1] + old_coord.z * m[0][2];
       new_coord.y = old_coord.x * m[1][0] + old_coord.y * m[1][1] + old_coord.z * m[1][2];
       new_coord.z = old_coord.x * m[2][0] + old_coord.y * m[2][1] + old_coord.z * m[2][2];
}

void matrix_rotate_translate(COORD &v, const double m[4][4])
{
       COORD old = v;
       matrix_rotate_translate(v, old, m);
}

void matrix_rotate_translate(COORD &new_coord, const COORD& old_coord, const double m[4][4])
{
       new_coord.x = old_coord.x * m[0][0] + old_coord.y * m[0][1] + old_coord.z * m[0][2] + m[0][3];
       new_coord.y = old_coord.x * m[1][0] + old_coord.y * m[1][1] + old_coord.z * m[1][2] + m[1][3];
       new_coord.z = old_coord.x * m[2][0] + old_coord.y * m[2][1] + old_coord.z * m[2][2] + m[2][3];
}

void multi_matrix(COORD o[3], const COORD a[3], const COORD b[3])
{
       for (int i = 0; i < 3; i++) {
            o[i].x = a[0].x * b[i].x + a[1].x * b[i].y + a[2].x * b[i].z;
            o[i].y = a[0].y * b[i].x + a[1].y * b[i].y + a[2].y * b[i].z;
            o[i].z = a[0].z * b[i].x + a[1].z * b[i].y + a[2].z * b[i].z;
       }
}

void multi_matrix(double result[3][3], const double matrix1[3][3], const double matrix2[3][3])
{
       for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                 result[i][j] = 0.0;
                 for (int k = 0; k < 3; k++) {
                      result[i][j] += (matrix1[i][k] * matrix2[k][j]);
                 }
            }
       }
}

void multi_matrix_4(double result[4][4], const double matrix1[4][4], const double matrix2[4][4])
{
       for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                 result[i][j] = 0.0;
                 for (int k = 0; k < 4; k++) {
                      result[i][j] += (matrix1[i][k] * matrix2[k][j]);
                 }
            }
       }
}

void multi_matrix_vector(COORD &new_coord, const COORD matrix[3], const COORD& old_coord)
{
       new_coord.x = old_coord.x * matrix[0].x + old_coord.y * matrix[0].y + old_coord.z * matrix[0].z;
       new_coord.y = old_coord.x * matrix[1].x + old_coord.y * matrix[1].y + old_coord.z * matrix[1].z;
       new_coord.z = old_coord.x * matrix[2].x + old_coord.y * matrix[2].y + old_coord.z * matrix[2].z;
}

void multi_matrix_vector(COORD &new_coord, const double matrix[3][3], const COORD& old_coord)
{
       new_coord.x = old_coord.x * matrix[0][0] + old_coord.y * matrix[0][1] + old_coord.z * matrix[0][2];
       new_coord.y = old_coord.x * matrix[1][0] + old_coord.y * matrix[1][1] + old_coord.z * matrix[1][2];
       new_coord.z = old_coord.x * matrix[2][0] + old_coord.y * matrix[2][1] + old_coord.z * matrix[2][2];
}

void multi_matrix_vector(double result[3], const double matrix[3][3], const double vector[3])
{
       for (int i = 0; i < 3; i++) {
            result[i] = 0.0;
            for (int j = 0; j < 3; j++) {
                 result[i] += (matrix[i][j] * vector[j]);
            }
       }
}

void transpose_matrix(COORD out[3], const COORD in[3])
{
       out[0].x = in[0].x;
       out[0].y = in[1].x;
       out[0].z = in[2].x;
       out[1].x = in[0].y;
       out[1].y = in[1].y;
       out[1].z = in[2].y;
       out[2].x = in[0].z;
       out[2].y = in[1].z;
       out[2].z = in[2].z;
}

double determinant_of_matrix_3(const double matrix[3][3])
{
       return (matrix[0][0] * matrix[1][1] * matrix[2][2] - matrix[0][0] * matrix[1][2] * matrix[2][1] +
               matrix[0][1] * matrix[2][0] * matrix[1][2] - matrix[0][1] * matrix[1][0] * matrix[2][2] +
               matrix[0][2] * matrix[1][0] * matrix[2][1] - matrix[0][2] * matrix[1][1] * matrix[2][0]);
}
