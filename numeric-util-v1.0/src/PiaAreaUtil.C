/*
FILE:     PiaAreaUtil.C
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

#include "PiaAreaUtil.h"

#define XBIG   1.0e+18
#define GAMUT  5.0e+8
#define GAMID  2.5e+8

typedef struct {
       long x;
       long y;
} point;

typedef struct {
       point ip;
       point rx;
       point ry;
       long inside;
} vertex;

static void _get_min_and_max(const std::vector<COORD>& a, double& minx, double& miny, double& maxx, double& maxy);
static std::vector<vertex> fit(const double& minx, const double& miny, const double& mid, const double& sclx, const double& scly,
                               const std::vector<COORD>& x, const long& fudge);
static bool ovl(const point& p, const point& q);
static long area(const point& a, const point& p, const point& q);
static double cross(vertex& a, const vertex& b, vertex& c, const vertex& d, const long& a1, const long& a2, const long& a3, const long& a4);
static double contrib(const point& f, const COORD& t, const long& w);
static double contrib(const COORD& f, const point& t, const long& w);
static double contrib(const point& f, const point& t, const long& w);
static double inness(const std::vector<vertex>& P, const std::vector<vertex>& Q);

double pia_area(const std::vector<COORD>& a, const std::vector<COORD>& b)
{
       if (a.size() < 3  || b.size() < 3) return 0.0;

       double minx = XBIG, miny = XBIG, maxx = -XBIG, maxy = -XBIG;
       _get_min_and_max(a, minx, miny, maxx, maxy);
       _get_min_and_max(b, minx, miny, maxx, maxy);

       double sclx = GAMUT / (maxx - minx);
       double scly = GAMUT / (maxy - miny);
       double ascale = sclx * scly;

       std::vector<vertex> ipa = fit(minx, miny, GAMID, sclx, scly, a, 0);
       std::vector<vertex> ipb = fit(minx, miny, GAMID, sclx, scly, b, 2);

       double out_s = 0.0;
       for (unsigned int j = 0; j < a.size(); ++j) {
            for (unsigned int k = 0; k < b.size(); ++k) {
                 if (ovl(ipa[j].rx, ipb[k].rx) && ovl(ipa[j].ry, ipb[k].ry)) {
                      long a1 = -area(ipa[j].ip, ipb[k].ip, ipb[k + 1].ip);
                      long a2 = area(ipa[j + 1].ip, ipb[k].ip, ipb[k + 1].ip);
                      if ((a1 < 0) == (a2 < 0)) {
                           long a3 = area(ipb[k].ip, ipa[j].ip, ipa[j + 1].ip);
                           long a4 = -area(ipb[k + 1].ip, ipa[j].ip, ipa[j + 1].ip);
                           if ((a3 < 0) == (a4 < 0)) {
                                if (a1 < 0)
                                     out_s += cross(ipa[j], ipa[j + 1], ipb[k], ipb[k + 1], a1, a2, a3, a4);
                                else out_s += cross(ipb[k], ipb[k + 1], ipa[j], ipa[j + 1], a3, a4, a1, a2);
                           }
                      }
                 }
            }
       }

       out_s += inness(ipa, ipb);
       out_s += inness(ipb, ipa);

       return fabs(out_s) / ascale;
}

static void _get_min_and_max(const std::vector<COORD>& a, double& minx, double& miny, double& maxx, double& maxy)
{
       for (std::vector<COORD>::const_iterator pos = a.begin(); pos != a.end(); ++pos) {
            if (minx > pos->x) minx = pos->x;
            if (miny > pos->y) miny = pos->y;
            if (maxx < pos->x) maxx = pos->x;
            if (maxy < pos->y) maxy = pos->y;
       }
}

static std::vector<vertex> fit(const double& minx, const double& miny, const double& mid, const double& sclx, const double& scly,
                               const std::vector<COORD>& x, const long& fudge)
{
       std::vector<vertex> ix(x.size() + 1);

       for (int c = (int) x.size() - 1; c >= 0; --c) {
            long t = (long) ((x[c].x - minx) * sclx - mid);
            ix[c].ip.x = (t & ~7) | fudge | (c & 1);
            t = (long) ((x[c].y - miny) * scly - mid);
            ix[c].ip.y = (t & ~7) | fudge;
       }

       ix[0].ip.y += x.size() & 1;
       ix[x.size()] = ix[0];

       point t1, t2;
       for (int c = (int) x.size() - 1; c >= 0; --c) {
            t1.x = ix[c].ip.x;
            t1.y = ix[c + 1].ip.x;
            t2.x = ix[c + 1].ip.x;
            t2.y = ix[c].ip.x;
            ix[c].rx = (ix[c].ip.x < ix[c + 1].ip.x) ? t1 : t2;

            t1.x = ix[c].ip.y;
            t1.y = ix[c + 1].ip.y;
            t2.x = ix[c + 1].ip.y;
            t2.y = ix[c].ip.y;
            ix[c].ry = (ix[c].ip.y < ix[c + 1].ip.y) ? t1 : t2;
            ix[c].inside = 0;
       }

       return ix;
}

static bool ovl(const point& p, const point& q)
{
       return p.x < q.y && q.x < p.y;
}

static long area(const point& a, const point& p, const point& q)
{
       return (p.x * q.y - p.y * q.x + a.x * (p.y - q.y) + a.y * (q.x - p.x));
}

static double cross(vertex& a, const vertex& b, vertex& c, const vertex& d, const long& a1, const long& a2, const long& a3, const long& a4)
{
       COORD dp;

       double r1 = (double) a1 / ((double) a1 + (double) a2);
       double r2 = (double) a3 / ((double) a3 + (double) a4);

       dp.x = (double) a.ip.x + r1 * ((double) b.ip.x - (double) a.ip.x);
       dp.y = (double) a.ip.y + r1 * ((double) b.ip.y - (double) a.ip.y);
       double out_s = contrib(dp, b.ip, 1);

       dp.x = (double) c.ip.x + r2 * ((double) d.ip.x - (double) c.ip.x);
       dp.y = (double) c.ip.y + r2 * ((double) d.ip.y - (double) c.ip.y);
       out_s += contrib(d.ip, dp, 1);

       ++a.inside;
       --c.inside;

       return out_s;
}

static double contrib(const point& f, const COORD& t, const long& w)
{
       return ((double) w * (t.x - (double) f.x) * (t.y + (double) f.y) * 0.5);
}

static double contrib(const COORD& f, const point& t, const long& w)
{
       return ((double) w * ((double) t.x - f.x) * ((double)  t.y + f.y) * 0.5);
}

static double contrib(const point& f, const point& t, const long& w)
{
       return ((double) w * ((double) t.x - (double) f.x) * ((double)  t.y + (double) f.y) * 0.5);
}

static double inness(const std::vector<vertex>& P, const std::vector<vertex>& Q)
{
       long s = 0;
       point p = P[0].ip;
       double out_s = 0;

       for (int c = (int) Q.size() - 2; c >= 0; --c) {
            if (Q[c].rx.x < p.x && p.x < Q[c].rx.y) {
                 bool sign = 0 < area(p, Q[c].ip, Q[c + 1].ip);
                 if (sign == (Q[c].ip.x < Q[c + 1].ip.x)) {
                      s += (sign ? -1 : 1);
                 }
            }
       }

       for (unsigned int j = 0; j < P.size() - 1; ++j) {
            if (s) out_s += contrib(P[j].ip, P[j + 1].ip, s);
            s += P[j].inside;
       }
       return out_s;
}
