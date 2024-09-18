/*
FILE:     Promotif_Helix.C
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

#include "Promotif.h"
#include "utillib.h"

void Promotif::_assignHELIXandSHEET()
{
       _hlxset(0, 4, 'H', 'h');

       for (unsigned int i = 0; i < _residues.size(); ++i) {
            if (_residues[i].ssbond[2] != ' ') {
                 if (_residues[i].ssbond[9] == ' ')
                      _residues[i].ssbond[9] = 'E';
                 if (_residues[i].ssbond[10] == ' ' ||
                     _residues[i].ssbond[10] == 'h' ||
                     _residues[i].ssbond[10] == 'e')
                      _residues[i].ssbond[10] = 'E';
                 if (_residues[i].ssbond[2] == 'e') {
                      if ((i && _residues[i-1].ssbond[10] == ' ') ||
                          _residues[i-1].ssbond[10] == 'B')
                           _residues[i-1].ssbond[10] = 'e';
                      if (((i+1) < _residues.size() && i &&
                          _residues[i+1].ssbond[10] == ' ') ||
                          _residues[i-1].ssbond[10] == 'B')
                           _residues[i+1].ssbond[10] = 'e';
                 }
                 if (_residues[i].ssbond[3] != ' ') {
                      if (_residues[i].ssbond[9] == ' ')
                           _residues[i].ssbond[9] = 'B';
                      if (_residues[i].ssbond[10] == ' ')
                           _residues[i].ssbond[10] = 'B';
                 }
            }
            if (i && _residues[i-1].chnbreak) _residues[i].ssbond[10] = ' ';
       }

       _hlxset(6, 3, 'G', 'g');
       _hlxset(7, 5, 'I', 'i');
       _hlxres('G', 'g', 'T', 't', 3);
       _hlxres('I', 'i', 'T', 't', 5);
}

void Promotif::_hlxset(const int& hlxpos, const int& hlxsz, const char& hlxch,
                       const char& altch)
{
       for (int i = 1; i < (int) _residues.size() - hlxsz; ++i) {
            if (_residues[i].ssbond[hlxpos] == '>' ||
                _residues[i].ssbond[hlxpos] == 'X') {
                 if (_residues[i-1].ssbond[hlxpos] == '>' ||
                     _residues[i-1].ssbond[hlxpos] == 'X') {
                      for (int j = 0; j < hlxsz; ++j) {
                           if (_residues[i+j].ssbond[9] == ' ')
                                _residues[i+j].ssbond[9] = hlxch;
                           if (_residues[i+j].ssbond[10] == ' ' ||
                               _residues[i+j].ssbond[10] == altch)
                                _residues[i+j].ssbond[10] = hlxch;
                      }
                      if (_residues[i-1].ssbond[10] == ' ')
                           _residues[i-1].ssbond[10] = altch;
                      if (_residues[i+hlxsz].ssbond[10] == ' ')
                           _residues[i+hlxsz].ssbond[10] = altch;
                 }
            }
       }
}

void Promotif::_hlxres(const char& hlxch, const char& altch, const char& turnch,
                       const char& alturn, const int& hlxsz)
{
       int index[2];
       index[0] = 9; index[1] = 10;
       for (int i = 0; i < 2; ++i) {
            int stpos = -1;
            for (int j = 0; j < (int) _residues.size(); ++j) {
                 if (_residues[j].ssbond[index[i]] == hlxch ||
                     _residues[j].ssbond[index[i]] == altch) {
                      if (stpos < 0) stpos = j;
                 } else {
                      if (stpos >= 0) {
                           if ((j - stpos) < hlxsz) {
                                for (int k = stpos; k < j; ++k) {
                                     if (_residues[k].ssbond[index[i]] == hlxch)
                                          _residues[k].ssbond[index[i]] = turnch;
                                     else _residues[k].ssbond[index[i]] = alturn;
                                }
                           }
                           stpos = -1;
                      }
                 }
            }
            if (stpos >= 0) {
                 if (((int) _residues.size() - stpos) < hlxsz) {
                      for (int k = stpos; k < (int) _residues.size(); ++k) {
                           if (_residues[k].ssbond[index[i]] == hlxch)
                                _residues[k].ssbond[index[i]] = turnch;
                           else _residues[k].ssbond[index[i]] = alturn;
                      }
                 }
            }
       }
}

void Promotif::_findHelices()
{
       int n = _residues.size() / 4 + 1;
       std::vector<int> start, stop;
       start.reserve(n);
       stop.reserve(n);
       for (int i = 0; i < n; ++i) {
            start.push_back(0);
            stop.push_back(0);
       }
       int start_count = 0;
       int stop_count = 0;
       for (unsigned int i = 0; i < _residues.size(); ++i) {
            if (i == 0 || (i && _residues[i-1].chnbreak)) {
                 if (_residues[i].ssbond[10] == 'H' ||
                     _residues[i].ssbond[10] == 'h' ||
                     _residues[i].ssbond[10] == 'G' ||
                     _residues[i].ssbond[10] == 'g')
                      start[start_count++] = i;
            } else {
                 if ((_residues[i].ssbond[10] == 'H' &&
                      _residues[i-1].ssbond[10] != 'h' &&
                      _residues[i-1].ssbond[10] != 'H') ||
                     (_residues[i].ssbond[10] == 'h' &&
                      _residues[i-1].ssbond[10] != 'H') ||
                     (_residues[i].ssbond[10] == 'G' &&
                      _residues[i-1].ssbond[10] != 'g' &&
                      _residues[i-1].ssbond[10] != 'G') ||
                     (_residues[i].ssbond[10] == 'g' &&
                      _residues[i-1].ssbond[10] != 'G')) {
                      if (start_count) {
                           if (start[start_count-1] != (int) i)
                                start[start_count++] = i;
                      } else start[start_count++] = i;
                 } else if (i < (_residues.size() - 1) &&
                     _residues[i].ssbond[10] == 'h' &&
                     _residues[i+1].ssbond[10] == 'H') {
                      if (start_count) {
                           if (start[start_count-1] != (int) i)
                                start[start_count++] = i;
                      } else start[start_count++] = i;
                 }
            }
       }

       for (unsigned int i = 0; i < _residues.size(); ++i) {
            if (i == (_residues.size() - 1) || _residues[i].chnbreak) {
                 if (_residues[i].ssbond[10] == 'H' ||
                     _residues[i].ssbond[10] == 'h' ||
                     _residues[i].ssbond[10] == 'G' ||
                     _residues[i].ssbond[10] == 'g') {
                      if (stop_count) {
                           if (stop[stop_count-1] != (int) i)
                                stop[stop_count++] = i;
                      } else stop[stop_count++] = i;
                 }
            } else {
                 if ((_residues[i].ssbond[10] == 'H' &&
                      _residues[i+1].ssbond[10] != 'h' &&
                      _residues[i+1].ssbond[10] != 'H') ||
                     (_residues[i].ssbond[10] == 'h' &&
                      _residues[i+1].ssbond[10] != 'H') ||
                     (_residues[i].ssbond[10] == 'G' &&
                      _residues[i+1].ssbond[10] != 'g' &&
                      _residues[i+1].ssbond[10] != 'G') ||
                     (_residues[i].ssbond[10] == 'g' &&
                      _residues[i+1].ssbond[10] != 'G')) {
                      if (stop_count) {
                           if (stop[stop_count-1] != (int) i)
                                stop[stop_count++] = i;
                      } else stop[stop_count++] = i;
                 } else if (i && _residues[i-1].ssbond[10] == 'H' &&
                     _residues[i].ssbond[10] == 'h') {
                      if (stop_count) {
                           if (stop[stop_count-1] != (int) i)
                                stop[stop_count++] = i;
                      } else stop[stop_count++] = i;
                 }
            }
       }
       if (start_count != stop_count) {
            std::string error = "Start_count=" + String::IntToString(start_count)
                              + " Stop_count=" + String::IntToString(stop_count)
                              + ".\n Counting does not match. Assign HELIX failed.\n";
            _logIo->messageError(error);
            return;
       }

       for (int i = 1; i < start_count; ++i) {
            if (start[i] == stop[i] && (stop[i-1]+1) == start[i] &&
                !_residues[stop[i-1]].chnbreak) {
                 if (((_residues[start[i]].ssbond[10] == 'G' ||
                       _residues[start[i]].ssbond[10] == 'g') &&
                      (_residues[start[i-1]].ssbond[10] == 'G' ||
                       _residues[start[i-1]].ssbond[10] == 'g')) ||
                     ((_residues[start[i]].ssbond[10] == 'H' ||
                       _residues[start[i]].ssbond[10] == 'h') &&
                      (_residues[start[i-1]].ssbond[10] == 'H' ||
                       _residues[start[i-1]].ssbond[10] == 'h'))) {
                      stop[i-1] = stop[i];
                      for (int j = i; j < start_count - 1; ++j) {
                           start[j] = start[j + 1];
                           stop[j] = stop[j + 1];
                      }
                      start_count--;
                 }
            }
       }

       for (int i = 0; i < start_count - 1; ++i) {
            if (start[i] == stop[i] && !_residues[start[i]].chnbreak &&
                start[i+1] == (stop[i]+1)) {
                 if (((_residues[start[i]].ssbond[10] == 'G' ||
                       _residues[start[i]].ssbond[10] == 'g') &&
                      (_residues[start[i+1]].ssbond[10] == 'G' ||
                       _residues[start[i+1]].ssbond[10] == 'g')) ||
                     ((_residues[start[i]].ssbond[10] == 'H' ||
                       _residues[start[i]].ssbond[10] == 'h') &&
                      (_residues[start[i+1]].ssbond[10] == 'H' ||
                      _residues[start[i+1]].ssbond[10] == 'h'))) {
                      stop[i] = stop[i+1];
                      for (int j = i + 1; j < start_count - 1; ++j) {
                           start[j] = start[j + 1];
                           stop[j] = stop[j + 1];
                      }
                      start_count--;
                 }
            }
       }

       std::string character_1, character_2;
       if (start_count > 6084) {
            character_1 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ123456789";
            character_2 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ123456789";
       } else {
            character_1 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
            character_2 = "123456789";
       }
       std::set<std::string> id_set;
       id_set.clear();

       _HELIX helix;
       helix.comment.clear();
       for (int i = 0; i < start_count; ++i) {
            helix.mol_index = _mol->index();
            helix.ID = get_next_id(character_1, character_1, character_2, id_set);
            helix.initRes = _residues[start[i]].res->index();
            helix.endRes = _residues[stop[i]].res->index();
            helix.helixClass = 1;
            if (_residues[stop[i]].ssbond[10] == 'G' ||
                _residues[stop[i]].ssbond[10] == 'g') helix.helixClass = 5;
            _helices.push_back(helix);
       }
}
