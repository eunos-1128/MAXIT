/*
FILE:     AngleUtil.C
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

#include "AngleUtil.h"

static double ndb_dihed_angle(double x1, double y1, double x2, double y2);

double cal_angle(const COORD& coord1, const COORD& coord2, const COORD& coord3)
{
       double dx = coord1.x - coord2.x;
       double dy = coord1.y - coord2.y;
       double dz = coord1.z - coord2.z;
       double dsqij = dx * dx + dy * dy + dz * dz;

       dx = coord2.x - coord3.x;
       dy = coord2.y - coord3.y;
       dz = coord2.z - coord3.z;
       double dsqjk = dx * dx + dy * dy + dz * dz;

       dx = coord1.x - coord3.x;
       dy = coord1.y - coord3.y;
       dz = coord1.z - coord3.z;
       double dsqik = dx * dx + dy * dy + dz * dz;

       double dxy = sqrt(dsqij * dsqjk);
       double t = (dsqij + dsqjk - dsqik) / (2.0 * dxy);
       if (t > 1.0) t = 1.0;
       if (t <-1.0) t =-1.0;

       dxy = (acos(t) * 180.0) / acos(-1.0);

       return (dxy);
}

double cal_torsion(const COORD& coord1, const COORD& coord2, const COORD& coord3, const COORD& coord4)
{
       double dxik1 = coord1.x - coord3.x;
       double dxjk1 = coord2.x - coord3.x;
       double dxlk1 = coord4.x - coord3.x;

       double dyik1 = coord1.y - coord3.y;
       double dyjk1 = coord2.y - coord3.y;
       double dylk1 = coord4.y - coord3.y;

       double dzik1 = coord1.z - coord3.z;
       double dzjk1 = coord2.z - coord3.z;
       double dzlk1 = coord4.z - coord3.z;

       double dist = sqrt(dxjk1 * dxjk1 + dyjk1 * dyjk1 + dzjk1 * dzjk1);
       double cosa = dzjk1 / dist;
       if (cosa > 1.0) cosa = 1.0;
       if (cosa < -1.0) cosa = -1.0;

       double ddd = 1.0 - cosa * cosa;

       double dxik2 = dxik1;
       double dxlk2 = dxlk1;
       double dyik2 = dyik1;
       double dylk2 = dylk1;
       double costh = cosa;
       double sinth = 0.0;

       if (ddd <= 0.0 ) {
            dxik2 = dxik1;
            dxlk2 = dxlk1;
            dyik2 = dyik1;
            dylk2 = dylk1;
            costh = cosa;
            sinth = 0.0;
       } else {
            double yxd = dist * sqrt(ddd);
            if (yxd > XEPS) {
                 double cosph = dyjk1 / yxd;
                 double sinph = dxjk1 / yxd;
                 dxik2 = dxik1 * cosph - dyik1 * sinph;
                 dxlk2 = dxlk1 * cosph - dylk1 * sinph;

                 dyik2 = dxik1 * sinph + dyik1 * cosph;
                 double dyjk2 = dxjk1 * sinph + dyjk1 * cosph;
                 dylk2 = dxlk1 * sinph + dylk1 * cosph;
                 costh = cosa;
                 sinth = dyjk2 / dist;
            } else {
                 dxik2 = dxik1;
                 dxlk2 = dxlk1;
                 dyik2 = dyik1;
                 dylk2 = dylk1;
                 costh = cosa;
                 sinth = 0.0;
            }
       }
       double dyik3 = dyik2 * costh - dzik1 * sinth;
       double dylk3 = dylk2 * costh - dzlk1 * sinth;

       double angle = ndb_dihed_angle(dxlk2, dylk3, dxik2, dyik3);
       angle = (angle * 180.0) / PI;
       return (angle);
}

static double ndb_dihed_angle(double x1, double y1, double x2, double y2)
{
       if (fabs(x1) < XEPS && fabs(y1) < XEPS) return(0.0);
       if (fabs(x2) < XEPS && fabs(y2) < XEPS) return(0.0);

       double anorm = 1.0 / sqrt(x1 * x1 + y1 * y1);
       double bnorm = 1.0 / sqrt(x2 * x2 + y2 * y2);
       x1 = x1 * anorm;
       y1 = y1 * anorm;
       x2 = x2 * bnorm;
       y2 = y2 * bnorm;
       double sinth = x1 * y2 - y1 * x2;
       double costh = x1 * x2 + y1 * y2;
       if (costh > 1.0) costh = 1.0;
       if (costh < -1.0) costh = -1.0;
       double angle = acos(costh);
       if (sinth > 0 ) angle = TWOPI - angle;
       angle = -angle;
       return (angle);
}
