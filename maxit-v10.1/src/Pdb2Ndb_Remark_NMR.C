/*
FILE:     Pdb2Ndb_Remark_NMR.C
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

#include "Maxit.h"
#include "utillib.h"

static void get_name_and_version(std::string& name, std::string& version, const std::string& program);
static void clean_author_name(std::string& authors, const std::string& name);
static void get_program_list(std::vector<std::string>& data, const std::string& string_buff);

void Maxit::_pdb_to_ndb_update_SFTWAR()
{
       std::string program, version, authors;
       std::vector<std::string> data, data1;

       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find("REFMET");
       if (ppos != _pdb_records.end()) {
            for (std::list<std::vector<std::string> >::const_iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 if ((*lpos)[3].empty() && (*lpos)[4].empty()) continue;

                 get_wordarray_delimit_by_string(data, (*lpos)[3], " AND ");
                 get_wordarray_delimit_by_string(data1, (*lpos)[4], " AND ");
                 if (data.size() == 1) {
                      get_wordarray_delimit_by_string(data, (*lpos)[3], ",");
                      get_wordarray_delimit_by_string(data1, (*lpos)[4], ",");
                      if (data.size() == 1) {
                           data1.clear();
                           data1.push_back((*lpos)[4]);
                      }
                 }
                 if (data1.size() < data.size()) {
                      for (unsigned int i = 0; i < data.size() - data1.size(); ++i) {
                           data1.push_back("");
                      }
                 }
                 for (unsigned int i = 0; i < data.size(); ++i) {
                      String::StripAndCompressWs(data[i]);
                      String::StripAndCompressWs(data1[i]);
                      if (data[i].empty() && data1[i].empty()) continue;
                      authors = data1[i];
                      get_name_and_version(program, version, data[i]);
                      clean_author_name(authors, program);
                      _addNewRecord("SFTWAR");
                      _updateRecordBack("SFTWAR", 2, "refinement");
                      _updateRecordBack("SFTWAR", 3, program);
                      _updateRecordBack("SFTWAR", 4, version);
                      _updateRecordBack("SFTWAR", 5, authors);
                 }
            }
            _pdb_records.erase("REFMET");
       }
 
       std::string value;
       _getRecordFront("NMRSFT", 2, value);
       if (!value.empty()) {
            get_program_list(data, value);
            for (std::vector<std::string>::const_iterator
                 pos = data.begin(); pos != data.end(); ++pos) {
                 get_name_and_version(program, version, *pos);
                 _addNewRecord("SFTWAR");
                 _updateRecordBack("SFTWAR", 2, "structure solution");
                 _updateRecordBack("SFTWAR", 3, program);
                 _updateRecordBack("SFTWAR", 4, version);
            }
            _pdb_records.erase("NMRSFT");
       }
}

void Maxit::_pdb_to_ndb_update_NMRSDT_and_NMREXP()
{
       std::vector<std::string> data, data1;

       int num_solution = 0;
       std::string value;
       _getRecordFront("NMRSDT", 2, value);
       if (!value.empty()) {
            get_wordlist_from_string_separated_by_delimits(data, value, ";", "MHZ", false);
            num_solution = data.size();
            _pdb_records.erase("NMRSDT");
            for (unsigned int i = 0; i < data.size(); ++i) {
                 _addNewRecord("NMRSDT");
                 _updateRecordBack("NMRSDT", 1, String::IntToString(i + 1));
                 _updateRecordBack("NMRSDT", 2, data[i]);
            }
       }

       _getRecordFront("NMRECD", 2, value);
       if (!value.empty()) {
            get_wordlist_from_string_separated_by_delimits(data, value, ";", "MHZ", false);
            int exp_id = 0;
            for (unsigned int i = 0; i < data.size(); ++i) {
                 get_wordlist_from_string_separated_by_delimits(data1, data[i], ",", "MHZ", false);
                 for (unsigned int j = 0; j < data1.size(); ++j) {
                      exp_id++;
                      _addNewRecord("NMREXP");
                      _updateRecordBack("NMREXP", 1, String::IntToString(exp_id));
                      _updateRecordBack("NMREXP", 2, String::IntToString(i + 1));
                      _updateRecordBack("NMREXP", 3, data1[j]);
                      int k = 1;
                      if ((int) data.size() == num_solution) k = i + 1;
                      _updateRecordBack("NMREXP", 4, String::IntToString(k));
                 }
            }
            _pdb_records.erase("NMRECD");
       }
}

void Maxit::_pdb_to_ndb_update_NMRSPM()
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find("NMRSPM");
       if (ppos == _pdb_records.end()) return;

       const std::vector<std::string>& Field = ppos->second.front();

       std::vector<std::string> data, data1, data2;

       get_wordlist_from_string_separated_by_delimits(data,  Field[2], ";", "MHZ", false);
       get_wordlist_from_string_separated_by_delimits(data1,  Field[3], ";", "MHZ", false);
       get_wordlist_from_string_separated_by_delimits(data2,  Field[4], ";", "MHZ", false);

       int max_number = (int) data.size();
       if ((int) data1.size() > max_number) max_number = (int) data1.size();
       if ((int) data2.size() > max_number) max_number = (int) data2.size();

       _pdb_records.erase("NMRSPM");
            
       for (int i = 0; i < max_number; ++i) {
            _addNewRecord("NMRSPM");
            _updateRecordBack("NMRSPM", 1, String::IntToString(i + 1));
            if (i < (int) data.size())  _updateRecordBack("NMRSPM", 2, data[i]);
            if (i < (int) data1.size()) _updateRecordBack("NMRSPM", 3, data1[i]);
            if (i < (int) data2.size()) _updateRecordBack("NMRSPM", 4, data2[i]);
       }
}

void Maxit::_pdb_to_ndb_update_NMRSMP()
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find("NMRSMP");
       if (ppos == _pdb_records.end()) return;

       const std::vector<std::string>& Field = ppos->second.front();

       std::vector<std::string> data, data1, data2, data3;

       get_wordlist_from_string_separated_by_delimits(data,  Field[3], ";", "MHZ", false);
       get_wordlist_from_string_separated_by_delimits(data1,  Field[4], ";", "MHZ", false);
       get_wordlist_from_string_separated_by_delimits(data2,  Field[5], ";", "MHZ", false);
       get_wordlist_from_string_separated_by_delimits(data3,  Field[6], ";", "MHZ", false);

       int max_number = (int) data.size();
       if ((int) data1.size() > max_number) max_number = (int) data1.size();
       if ((int) data2.size() > max_number) max_number = (int) data2.size();
       if ((int) data3.size() > max_number) max_number = (int) data3.size();

       _pdb_records.erase("NMRSMP");
       
       for (int i = 0; i < max_number; ++i) {
            _addNewRecord("NMRSMP");
            _updateRecordBack("NMRSMP", 1, String::IntToString(i + 1));
            _updateRecordBack("NMRSMP", 2, String::IntToString(i + 1));
            if (i < (int) data.size())  _updateRecordBack("NMRSMP", 3, data[i]);
            if (i < (int) data1.size()) _updateRecordBack("NMRSMP", 4, data1[i]);
            if (i < (int) data2.size()) _updateRecordBack("NMRSMP", 5, data2[i]);
            if (i < (int) data3.size()) _updateRecordBack("NMRSMP", 6, data3[i]);
       }
}

void Maxit::_pdb_to_ndb_update_EMSFTW()
{
       std::vector<std::string> data;
       std::string value;
       _getRecordFront("CMRFMT", 1, value);
       if (value.empty()) return;

       get_wordarray(data, value, ", ");
       for (std::vector<std::string>::const_iterator pos = data.begin(); pos != data.end(); ++pos) {
            _addNewRecord("EMSFTW");
            _updateRecordBack("EMSFTW", 2, *pos);
       }
}

void Maxit::_pdb_to_ndb_update_EMPIXL()
{
       std::vector<std::string> data;
       std::string value;
       _getRecordFront("CMRNST", 2, value);
       if (value.empty()) return;

       get_wordarray(data, value, ", ");
       if (data.size() != 1 && data.size() != 3) return;

       _addNewRecord("EMPIXL");
       if (data.size() == 3) {
            _updateRecordBack("EMPIXL", 1, "primary");
            _updateRecordBack("EMPIXL", 2, data[0]);
            _updateRecordBack("EMPIXL", 3, data[1]);
            _updateRecordBack("EMPIXL", 4, data[2]);
       } else {
            _updateRecordBack("EMPIXL", 1, "primary");
            _updateRecordBack("EMPIXL", 2, data[0]);
            _updateRecordBack("EMPIXL", 3, data[0]);
            _updateRecordBack("EMPIXL", 4, data[0]);
       }
}

static void get_name_and_version(std::string& name, std::string& version, const std::string& program)
{
       name = program;
       version.clear();

       std::vector<std::string> data;
       get_wordarray(data, program, " ");
       if (data.size() > 1) {
            bool found = false;
            std::string cs = data[data.size()-1];
            for (unsigned int i = 0; i < cs.size(); ++i) {
                 if (isdigit(cs[i])) {
                      found = true;
                      break;
                 }
            }
            if (found) {
                 version = data[data.size()-1];
                 std::string::size_type p = name.find(version);
                 name.erase(p);
                 String::StripAndCompressWs(name);
            }
       }
}

static void clean_author_name(std::string& authors, const std::string& name)
{
       std::string buff = "(" + name + ")";
       std::string::size_type p = authors.find(buff);
       if (p != std::string::npos) {
            authors.erase(p);
            String::StripAndCompressWs(authors);
            return;
       }
       p = authors.find(name);
       if (p != std::string::npos) {
            authors.erase(p);
            String::StripAndCompressWs(authors);
       }
}

static void get_program_list(std::vector<std::string>& wordlist, const std::string& string_buff)
{
       wordlist.clear();
       std::string cs = string_buff;
       String::StripAndCompressWs(cs);
       if (cs.empty()) return;

       std::vector<std::string> data, data1, data2;
       get_wordarray(data, cs, ",");
       for (unsigned int i = 0; i < data.size(); ++i) {
            get_wordarray_delimit_by_string(data1, data[i], " AND ");
            for (unsigned int j = 0; j < data1.size(); ++j) {
                 get_wordarray_delimit_by_string(data2, data1[j], " and ");
                 for (unsigned int k = 0; k < data2.size(); ++k) {
                      String::StripAndCompressWs(data2[k]);
                      if (!data2[k].empty()) wordlist.push_back(data2[k]);
                 }
            }
       }
}
