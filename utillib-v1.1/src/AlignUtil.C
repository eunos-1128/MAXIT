/*
FILE:     AlignUtil.C
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

#include "AlignUtil.h"

#ifndef MAX
#define MAX(X, Y)   (((X) > (Y))? (X): (Y))
#endif
#ifndef MIN
#define MIN(X, Y)   (((X) < (Y))? (X): (Y))
#endif

void AlignUtil::initAssignment(std::vector<std::vector<int> >& ss, const int& size)
{
       ss.clear();
       ss.reserve(size);

       std::vector<int> data;
       data.clear();
       data.push_back(0);
       for (int i = 0; i < size; ++i) {
            data[0] = i;
            ss.push_back(data);
       }
}

void AlignUtil::Alignment(std::vector<std::vector<int> >& ss, const std::vector<std::vector<int> >& sa, const std::vector<std::vector<int> >& sb,
                          const double& gap_penalty, double (*sfunc)(const int&, const int&, void*), void* data)
{
       std::vector<SAVE> dr;
       dr.clear();
       int n = 60 * sb.size();
       dr.reserve(n);
       int origin = dynamic(sa, sb, gap_penalty, dr, sfunc, data);
       mu_trcback(origin, ss, sa, sb, dr);
}

void AlignUtil::Alignment(std::vector<std::vector<int> >& ss, const std::vector<std::vector<int> >& sa, const std::vector<std::vector<int> >& sb,
                          const double& vv, const double& uu, double (*sfunc)(const int&, const int&, void*), void* data,
                          const std::vector<bool>& linkage, const double& factor)
{
       std::vector<SAVE> dr;
       dr.clear();
       int n = 60 * sb.size();
       dr.reserve(n);
       int origin = dynamic(sa, sb, vv, uu, dr, sfunc, data, linkage, factor);
       mu_trcback(origin, ss, sa, sb, dr);
}

int AlignUtil::dynamic(const std::vector<std::vector<int> >& sa, const std::vector<std::vector<int> >& sb, const double gap_penalty, std::vector<SAVE>& dr,
                       double (*sfunc)(const int&, const int&, void*), void* data)
{
       int m = sa.size();
       int n = sb.size();
       int up, lw;
       stripe(m, n, up, lw);

       RECD *dd = new RECD[n + 1];
       RECD *p  = new RECD[n + 1];
       RECD *q  = new RECD[n + 1];

       int k = adr(0, 0, -1, dr);
       for (int i = 0; i <= n; i++) {
            p[i].val = 0; p[i].ptr = k;
            q[i].val = 0; q[i].ptr = k;
            dd[i].val = 0; dd[i].ptr = k;
       }
       int j = 0;
       double x;
       for (int i = 1; i <= m; i++) {
            int n1 = i + lw + 1;
            int n2 = i + up;
            j = MAX(n1, 1);
            int n9 = MIN(n2, n);
            RECD *d = dd +j - 1;
            RECD *f = q + j - 1;
            RECD *g = p + j - 1;
            if (j == 1) {
               d->val = 0; d->ptr = k;
               f->val = 0; f->ptr = k;
            }
            RECD gg = *g;
            RECD ff = *f;
            RECD dt = *d;
            for (d++, g++, f++; j <= n9; d++, g++, f++, j++) {
                 RECD nd = ff;
                 if (gg.val > nd.val) nd = gg;
                 nd.val -= gap_penalty;
                 gg = *g;
                 ff = *f;
                 RECD tt = *d;
                 if ((d-1)->val >= (f-1)->val) {
                      f->val = (d-1)->val; f->ptr = (d-1)->ptr;
                 } else {
                      f->val = (f-1)->val; f->ptr = (f-1)->ptr;
                 }
                 if ((x = d->val) >= g->val) {
                      g->val = x; g->ptr = d->ptr;
                 }
                 d->val = dt.val;
                 if (nd.val > d->val) {
                      d->val = nd.val;
                      d->ptr = adr(i - 1, j - 1, nd.ptr, dr);
                 } else d->ptr = adr(i - 1, j - 1, dt.ptr, dr);
                 dt = tt;
                 d->val += (*sfunc)(sa[i-1][0], sb[j-1][0], data);
            }
       }
       x = dd[0].val;
       for (int i = 1; i <= n; i++)
            if (dd[i].val > x) {
                 x = dd[i].val;
                 j = dd[i].ptr;
            }
       for (int i = 1; i <= n; i++)
            if (p[i].val > x) {
                 x = p[i].val;
                 j = p[i].ptr;
            }
       for (int i = 1; i <= n; i++)
            if (q[i].val > x) {
                 x = q[i].val;
                 j = q[i].ptr;
            }

       delete [] dd;
       delete [] p;
       delete [] q;

       return(adr(m, n, j, dr));
}

int AlignUtil::dynamic(const std::vector<std::vector<int> >& sa, const std::vector<std::vector<int> >& sb, const double& vv, const double& uu,
                       std::vector<SAVE>& dr, double (*sfunc)(const int&, const int&, void*), void *data, const std::vector<bool>& linkage,
                       const double& factor)
{
       int m = sa.size();
       int n = sb.size();

       RECD *dd = new RECD[n + 1];
       RECD *gg = new RECD[n + 1];
       char *dir = new char[m + n + 1];

       int origin = adr(0, 0, -1, dr);
       for (int i = 0; i < (m + n + 1); i++) dir[i] = 0;
       for (int i = 0; i <= n; i++) {
            dd[i].val = 0; dd[i].ptr = origin;
            gg[i].val = 0; gg[i].ptr = origin;
       }

       // double /* big_vv = vv;
       // if (!linkage.empty()) */ big_vv = 10 * vv;
       double big_vv = factor * vv;
       int k = 0;
       double x, max_val = 0, maximum = -1.0e38;
       RECD  nd, f, dt, dtt;
       for (int i = 0; i < m; i++) {
            int j  = 0;
            RECD *d = dd + j;
            RECD *g = gg + j;
            d->val = f.val = 0;
            d->ptr = f.ptr = origin;
            int r = j - i + m;
            dt = *d;
            for (d++, g++; j < n; d++, g++, r++, j++) {
                 if ((x = (d-1)->val - big_vv - uu) >= (f.val -= uu)) {
                      f.val = x; f.ptr = (d-1)->ptr;
                 }
                 double vv1 = vv;
                 if (!linkage.empty() && linkage[j]) vv1 = big_vv;
                 if ((x = d->val - vv1 - uu) >= (g->val -= uu)) {
                      g->val = x; g->ptr = d->ptr;
                 }
                 nd.val = maximum;
                 nd.ptr = -1;
                 if (f.val > nd.val) nd = f;
                 if (g->val > nd.val) nd = *g;
                 dt.val += (*sfunc)(sa[i][0], sb[j][0], data);
                 dtt = *d;
                 *d = dt;
                 if (nd.val > d->val) {
                      *d = nd;
                      dir[r] = 0;
                 } else if (!dir[r]) {
                      d->ptr = adr(i, j, d->ptr, dr);
                      dir[r] = 1;
                 }
                 if (d->val > max_val) {
                      max_val = d->val;
                      k = d->ptr;
                 }
                 dt = dtt;
            }
       }

       delete [] dd;
       delete [] gg;
       delete [] dir;

       return (adr(m, n, k, dr));
}

void AlignUtil::mu_trcback(const int& origin, std::vector<std::vector<int> >& ss, const std::vector<std::vector<int> >& a,
                           const std::vector<std::vector<int> >& b, const std::vector<SAVE>& dr)
{
       ss.clear();

       std::list<std::pair<int, int> > top;
       top.clear();

       int i = origin;
       while (i >= 0) {
            top.push_front(std::make_pair(dr[i].m, dr[i].n));
            i = dr[i].pp;
       }
       if (top.empty()) return;

       int n1 = 0;
       int n2 = 0;
       int j, m, n, r;
       std::list<std::pair<int, int> >::iterator pos = top.begin();
       std::list<std::pair<int, int> >::iterator pos1 = top.begin();
       pos1++;
       while (pos1 != top.end()) {
	    if ((i = pos->first) - (j = pos->second) > 
                (m = pos1->first) - (n = pos1->second) && 
                (r = m - i + j) != j) {
                 n1 += (n - r); 
                 top.insert(pos1, std::make_pair(m, r));
            } else if ((i = pos->first) - (j = pos->second) < 
               (m = pos1->first) - (n = pos1->second) && 
               (r = n + i - j) != i) {
                 n2 += (m - r); 
                 top.insert(pos1, std::make_pair(r, n));
            }
            pos = pos1;
            pos1++;
       }

       m = a.size() + n1;
       n = b.size() + n2;
       int n9 = MAX(m, n);
       ss.reserve(n9);
       std::vector<int> data;
       data.clear();
       data.push_back(-1);
       data.push_back(-1);
       for (i = 0; i < n9; ++i) ss.push_back(data);

       track(ss, top, a, b);
}

void AlignUtil::track(std::vector<std::vector<int> >& aa, std::list<std::pair<int, int> >& top,
                      const std::vector<std::vector<int> >& sa, const std::vector<std::vector<int> >& sb)
{
       int x, y, m, n, r;

       r = 0;
       std::list<std::pair<int, int> >::iterator pos = top.begin();
       std::list<std::pair<int, int> >::iterator pos1 = top.begin();
       pos1++;
       while (pos1 != top.end()) {
	    if ((x = pos->first) - (y = pos->second) == 
                (m = pos1->first) - (n = pos1->second)) {
		 for (int i = 0; i < m - x; ++i) {
                      aa[i + r][0] = sa[x + i][0];
                      aa[i + r][1] = sb[y + i][0];
		 }
		 r += (m - x);
	    } else if(x - y < m - n) {
		 for (int i = 0; i < m - x; ++i) {
                      aa[i + r][0] = sa[x + i][0];
		 }
		 r += (m - x);
	    } else {
		 for (int i = 0; i < n - y; ++i) {
                      aa[i + r][1] = sb[y + i][0];
		 }
		 r += (n - y);
	    }
            pos = pos1;
            pos1++;
       }
}

int AlignUtil::adr(const int& m, const int& n, const int& pp, std::vector<SAVE>& dr)
{
       int i_record = dr.size();

       SAVE s;
       s.m = m;
       s.n = n;
       s.pp = pp;
       dr.push_back(s);  

       return (i_record);
}

void AlignUtil::stripe(const int& a, const int& b, int& up, int& lw)
{
       int  p;

       up = b - a; lw = 0;
       if (up < lw) {p = up; up = lw; lw = p;}
       up += 25; lw -= 25;
       if ((p = b - 1) < up) up = p + 1;
       if ((p = 1 - a) > lw) lw = p - 1;
}
