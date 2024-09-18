/*
FILE:     Promotif_Assignment.C
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
#include <stdlib.h>
#include "Promotif.h"

#define BKDSZL    5
#define BLGSZS    2
#define NSTCHS   26
#define BRIDG1    4
#define BRIDG2    5

static const char stchap[27] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
static const char stchpa[27] = "abcdefghijklmnopqrstuvwxyz";

void Promotif::_assignSecondaryStruct()
{
       int turnpl[3] = {6, 0, 7};
       int turnsz[3] = {3, 4, 5};
       const char turnch[4] = "345";

       for (unsigned int i = 0; i < _residues.size(); ++i) {
            for (int j = 0; j < 2; ++j) {
                 for (int k = 0; k < 3; ++k) {
                      if (_residues[i].hbond[j]) {
                           int l = _residues[i].hbond[j] - 1;
                           if ((l - (int) i) == turnsz[k] &&
                               _residues[l].chn == _residues[i].chn) {
                                if (_residues[i].ssbond[turnpl[k]] == '<')
                                     _residues[i].ssbond[turnpl[k]] = 'X';
                                else _residues[i].ssbond[turnpl[k]] = '>';
                                _residues[l].ssbond[turnpl[k]] = '<';
                                for (l = (int) i + 1; l < (int) i + turnsz[k]; ++l) {
                                     if (_residues[l].ssbond[turnpl[k]] == ' ')
                                          _residues[l].ssbond[turnpl[k]] = turnch[k];
                                }
                                for (l = (int) i; l < (int) i + turnsz[k]; ++l)
                                     _residues[l].ssbond[8] = 'T';
                           }
                      }
                 }
            }
       }

       int bsrng[3][2] = { { 2, -1 }, { 1, 0 }, { 2, -1 } };
       int bsofst[3][4] = { { -1, 0, 0, 1 }, { 0, 0, 0, 0 }, { -1, -1, -2, 1 } };
       int pbrdg[3] = { 1, -1, -1 };
       int exrule[2][3] = { { -1, 1, -1 }, { 1, 1, -1 } };

       std::vector<int> brdg;
       brdg.reserve(_residues.size() * 4);
       brdg.clear();
       for (unsigned int i = 0; i < _residues.size() * 4; ++i) brdg.push_back(0);

       for (int i = 0; i < 3; ++i) {
            for (int j = bsrng[i][0] - 1; j < (int) _residues.size() + bsrng[i][1]; ++j) {
                 for (int k = 0; k < 2; ++k) {
                      int bndptr = _residues[j + bsofst[i][0]].hbond[k] - 1;
                      if (!i || bndptr > j) {
                           int brdgpt = bndptr + bsofst[i][1];
                           if ((bndptr+1) > -bsofst[i][2] && abs(j - brdgpt) > 2) {
                                for (int l = 0; l < 2; ++l) {
                                     if (_residues[bndptr + bsofst[i][2]].hbond[l] ==
                                        (j + bsofst[i][3] + 1)) {
                                          int brpl = 1;
                                          if (brdg[4 * j + brpl - 1]) brpl = 2;
                                          brdg[4 * j + brpl - 1] = pbrdg[i] * (brdgpt + 1);
                                          brdg[4 * j + brpl + 1] = exrule[0][i];
                                          int brpla = 1;
                                          if (brdg[4 * brdgpt + brpla - 1]) brpla = 2;
                                          brdg[4 * brdgpt + brpla - 1] = pbrdg[i] * (j + 1);
                                          brdg[4 * brdgpt + brpla + 1] = exrule[1][i];
                                     }
                                }
                           }
                      }
                 }
            }
       }

       int nxstrn = 0;
       std::vector<int> shtcod, strcd;
       shtcod.reserve(_residues.size());
       strcd.reserve(_residues.size() * 2);
       shtcod.clear();
       strcd.clear();
       for (unsigned int i = 0; i < _residues.size(); ++i) {
            shtcod.push_back(0);
            strcd.push_back(0);
            strcd.push_back(0);
            _residues[i].bridpt[0] = 0;
            _residues[i].bridpt[1] = 0;
       }

       for (unsigned int i = 0; i < _residues.size(); ++i) {
            for (int k = 0; k < 2; ++k) {
                 if (brdg[i * 4 + k] && abs(brdg[i * 4 + k]) > ((int) i + 1)) {
                      int brptr = abs(brdg[i * 4 + k]);
                      int brdir = brdg[i * 4 + k] / brptr;
                      int brptpl = 1;
                      if (abs(brdg[(brptr - 1) * 4]) == ((int) i + 1)) brptpl = 0;
                      int brpla = _residues.size();
                      if (((int) i + BKDSZL + 1) < brpla) brpla = i + BKDSZL + 1;
                      bool found = false;
                      for (int j = (int) i + 1; j < brpla; ++j) {
                           for (int l = 0; l < 2; ++l) {
                                if (brdg[j * 4 + l] * brdir > 0) {
                                     int newptr = abs(brdg[j * 4 + l]);
                                     if (abs(newptr - brptr) && ((newptr - brptr) /
                                         abs(newptr - brptr)) == brdir &&
                                        (((j - (int) i) > BLGSZS &&
                                         abs(newptr - brptr) <= BLGSZS) ||
                                         ((j - (int) i) <= BLGSZS &&
                                         abs(newptr - brptr) <= BKDSZL))) {
                                          int nwptpl = 1;
                                          if (abs(brdg[(newptr - 1) * 4]) == (j + 1))
                                               nwptpl = 0;
                                          for (int m = i; m <= j; ++m) {
                                               if (_residues[m].ssbond[2] == ' ')
                                                    _residues[m].ssbond[2] = 'E';
                                          }
                                          if (brdg[i * 4 + k + 2] < 0)
                                               _residues[i].ssbond[2] = 'e';
                                          if (brdg[j * 4 + l + 2] < 0)
                                               _residues[j].ssbond[2] = 'e';
                                          int stsht = brptr - 1;
                                          int finsht   = newptr - 1;
                                          if (stsht > finsht) {
                                               stsht = newptr - 1;
                                               finsht = brptr - 1;
                                          }
                                          for (int m = stsht; m <= finsht; ++m) {
                                               if (_residues[m].ssbond[2] == ' ')
                                                    _residues[m].ssbond[2] = 'E';
                                          }
                                          if (brdg[(brptr - 1) * 4 + brptpl + 2] < 0)
                                               _residues[brptr - 1].ssbond[2] = 'e';
                                          if (brdg[(newptr - 1) * 4 + nwptpl + 2] < 0)
                                               _residues[newptr - 1].ssbond[2] = 'e';
                                          if (!strcd[2 * i + k]) {
                                               nxstrn++;
                                               strcd[2 * i + k] = nxstrn * brdir;
                                               strcd[2 * (brptr - 1) + brptpl] =
                                                                  nxstrn * brdir;
                                          }
                                          strcd[2 * j + l] = strcd[2 * i + k];
                                          strcd[2*(newptr-1)+nwptpl] =
                                                             strcd[2*(brptr-1)+brptpl];
                                          found = true;
                                          break;
                                     }
                                }
                           }
                           if (found) break;
                      }
                      char wkch = 'b';
                      if (brdir < 0) wkch = 'B';
                      _residues[i].ssbond[3] = wkch;
                      _residues[brptr - 1].ssbond[3] = wkch;
                 }
            }
       }

       int stsht = -1;
       int finsht = -1;
       for (unsigned int i = 0; i < _residues.size(); ++i) {
            if (_residues[i].ssbond[2] == ' ') {
                 if (stsht >= 0 && finsht >= 0) {
                      finsht++;
                      _stbrg(strcd, brdg, stsht, finsht);
                 }
                 stsht = -1;
                 finsht = -1;
            } else {
                 if (stsht < 0) stsht = i;
                 finsht = i;
            }
       }
       if (stsht >= 0 && finsht >= 0) {
            finsht++;
            _stbrg(strcd, brdg, stsht, finsht);
       }

       int shtsum = 0;
       for (unsigned int i = 0; i < _residues.size(); ++i) {
            _fdshus(i, shtcod, stsht, finsht);
            if (stsht < 0) break;

            shtsum++;
            for (int j = stsht; j < finsht; ++j) {
                 shtcod[j] = shtsum;
                 for (int l = 0; l < 2; ++l) {
                      if (_residues[j].bridpt[l]) {
                           if (!shtcod[_residues[j].bridpt[l] - 1]) {
                                _setsht(_residues[j].bridpt[l] - 1, shtsum, shtcod);
                           }
                      }
                 }
            }

            bool found = true;
            while (found) {
                 found = _fdshps(finsht, shtsum, shtcod);
            }
            i = finsht;
       }

       for (unsigned int i = 0; i < _residues.size(); ++i) {
            if (shtcod[i])
                 _residues[i].ssbond1 = shtcod[i];
            else _residues[i].ssbond1 = 0;
       }

       shtsum = 0;
       for (unsigned int i = 0; i < _residues.size(); ++i) {
            for (int j = 0; j < 2; ++j) {
                 if (brdg[i * 4 + j] && !strcd[i * 2 + j] &&
                     abs(brdg[i * 4 + j]) > (int) (i + 1)) {
                      shtsum++;
                      int brptr = abs(brdg[i * 4 + j]);
                      int chcode = shtsum % NSTCHS;
                      if (!chcode) chcode = NSTCHS;
                      char strdch = stchpa[chcode - 1];
                      if (brdg[i * 4 + j] < 0)
                           strdch = stchap[chcode - 1];
                      int stchpl = BRIDG1;
                      if (_residues[i].ssbond[BRIDG1] != ' ')
                           stchpl = BRIDG2;
                      _residues[i].ssbond[stchpl] = strdch;
                      _residues[i].bridpt[stchpl - BRIDG1] = brptr;
                      stchpl = BRIDG1;
                      if (_residues[brptr-1].ssbond[BRIDG1] != ' ')
                           stchpl = BRIDG2;
                      _residues[brptr-1].ssbond[stchpl] = strdch;
                      _residues[brptr-1].bridpt[stchpl - BRIDG1] = i + 1;
                 }
            }
       }
}

void Promotif::_stbrg(const std::vector<int>& strcd, const std::vector<int>& brdg,
                      const int& stsht, const int& finsht)
{
       int j = 0, l = 0, lststd = 0, bstval = 0;

       while (1) {
            _whstrd(strcd, stsht, finsht, j, l, lststd, bstval);
            if (bstval) {
                 lststd = abs(strcd[2 * j + l]);
                 int chcode = lststd % NSTCHS;
                 if (!chcode) chcode = NSTCHS;
                 char strdch = stchpa[chcode - 1];
                 if (strcd[2 * j + l] < 0)
                      strdch = stchap[chcode - 1];
                 int stchpl = BRIDG1;
                 if (_residues[j].ssbond[BRIDG1] != ' ') stchpl = BRIDG2;
                 if (stchpl != BRIDG2) {
                      int thstpl = j;
                      do {
                           int nxstpl = _strdpl(thstpl, finsht, strcd[2 * j + l], strcd);
                           if (nxstpl < 0) break;
                           if (_residues[nxstpl].ssbond[BRIDG1] != ' ')
                                stchpl = BRIDG2;
                           thstpl = nxstpl;
                      } while (stchpl != BRIDG2);
                 }
                 _residues[j].ssbond[stchpl] = strdch;
                 _residues[j].bridpt[stchpl - BRIDG1] = abs(brdg[j * 4 + l]);
                 int thstpl = j;
                 while (1) {
                      int nxstpl = _strdpl(thstpl, finsht, strcd[2 * j + l], strcd);
                      if (nxstpl < 0) break;
                      for (int m = thstpl + 1; m < nxstpl; ++m)
                           _residues[m].ssbond[stchpl] = '*';
                      _residues[nxstpl].ssbond[stchpl] = strdch;
                      int brptpl = 0;
                      if (strcd[2 * nxstpl] != strcd[2 * j + l]) brptpl = 1;
                      _residues[nxstpl].bridpt[stchpl - BRIDG1] =
                                       abs(brdg[nxstpl * 4 + brptpl]);
                      thstpl = nxstpl;
                 }
            } else break;
       }
}

void Promotif::_whstrd(const std::vector<int>& strcd, const int& stsht, const int& finsht,
                       int &stres, int &stpl, const int& lststd, int &bstval)
{
       bstval = 0;
       for (int i = stsht; i < finsht; ++i) {
            for (int j = 0; j < 2; ++j) {
                 if (abs(strcd[2 * i + j]) > lststd) {
                      if (!bstval) {
                           bstval = abs(strcd[2 * i + j]);
                           stres = i;
                           stpl = j;
                      } else {
                           if (abs(strcd[2 * i + j]) < bstval) {
                                bstval = abs(strcd[2 * i + j]);
                                stres = i;
                                stpl = j;
                           }
                      }
                 }
            }
       }
}

int Promotif::_strdpl(const int& stpl, const int& finpl, const int& thcode,
                      const std::vector<int>& strand)
{
       for (int i = stpl + 1; i < finpl; ++i) {
            if (strand[2 * i] == thcode) return i;
            if (strand[2*i+1] == thcode) return i;
       }
       return -1;
}

void Promotif::_fdshus(const int& stpt, const std::vector<int>& shtcod, int &stsht,
                       int &finsht)
{
       stsht = -1;
       for (int i = stpt; i < (int) _residues.size(); ++i) {
            if (_residues[i].ssbond[2] != ' ' && !shtcod[i]) {
                 stsht = i;
                 break;
            }
       }
       if (stsht < 0) return;
       finsht = -1;
       for (int i = stsht; i < (int) _residues.size(); ++i) {
            if (_residues[i].ssbond[2] == ' ') break;
            finsht = i;
       }
       finsht++;
}

void Promotif::_setsht(const int& refpl, const int& shtsum, std::vector<int>& shtcod)
{
       int finsht = -1;
       for (int i = refpl; i < (int) _residues.size(); ++i) {
            if (_residues[i].ssbond[2] == ' ') break;
            finsht = i;
       }
       int stsht = 0;
       for (int i = refpl; i >= 0; i--) {
            if (_residues[i].ssbond[2] == ' ') break;
            stsht = i;
       }
       for (int i = stsht; i <= finsht; ++i) shtcod[i] = -shtsum;
}

bool Promotif::_fdshps(const int& stpt, const int& shtsum, std::vector<int>& shtcod)
{
       bool found = false;
       int stsht = -1;
       for (int i = stpt; i < (int) _residues.size(); ++i) {
            if (shtcod[i] == (-shtsum)) {
                 found = true;
                 stsht = i;
                 break;
            }
       }
       if (!found) return false;

       int finsht = -1;
       for (int i = stsht; i < (int) _residues.size(); ++i) {
            if (shtcod[i] != (-shtsum)) break;
            finsht = i;
       }
       for (int i = stsht; i <= finsht; ++i) {
            shtcod[i] = shtsum;
            for (int j = 0; j < 2; ++j) {
                 if (_residues[i].bridpt[j]) {
                      if (!shtcod[_residues[i].bridpt[j] - 1]) {
                           _setsht(_residues[i].bridpt[j] - 1, shtsum, shtcod);
                      }
                 }
            }
       }
       return true;
}
