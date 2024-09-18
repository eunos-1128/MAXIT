/*
FILE:     NdbRefn.C
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
#include <sys/types.h>
#include <sys/stat.h>
#include <ctype.h>

#include "CifFile.h"
#include "NdbRefn.h"
#include "NdbRefn_global.h"
#include "utillib.h"

void NdbRefn::initialize()
{
       _all_refs.clear();
       _abbrev_mapping.clear();
       _issn_mapping.clear();
}

bool NdbRefn::Read(LogUtil& logIo, const std::string& path)
{
       std::string binaryfile = path + "/data/binary/";
       binaryfile += refbinfile;
       struct stat statbuf;
       if (stat(binaryfile.c_str(), &statbuf) != 0) {
            std::string error = "NdbRefn::Read: Can't find "
                              + binaryfile + "\n";
            logIo.message(error.c_str());
            return false;
       }

       CifFile* fobjR = new CifFile(READ_MODE, binaryfile);
       if (!fobjR) {
            std::string error = "NdbRefn::Read: Can't open "
                              + binaryfile + "\n";
            logIo.message(error.c_str());
            return false;
       }

       Block &block = fobjR->GetBlock(fobjR->GetFirstBlockName());

       ISTable *t = getTablePtr(block, "ndb_refn");
       if (!t) {
            std::string error = "NdbRefn::Read: Can't find ndb_refn category in "
                              + binaryfile + "\n";
            logIo.message(error.c_str());
            return false;
       }

       read_ndb_refn_table(t);
       create_indices();

       delete fobjR;

       return true;
}

bool NdbRefn::Read(FILE* log, const std::string& path)
{
       std::string binaryfile = path + "/data/binary/";
       binaryfile += refbinfile;
       struct stat statbuf;
       if (stat(binaryfile.c_str(), &statbuf) != 0) {
            fprintf(log, "NdbRefn::Read: Can't find %s\n", binaryfile.c_str());
            return false;
       }

       CifFile* fobjR = new CifFile(READ_MODE, binaryfile);
       if (!fobjR) {
            fprintf(log, "NdbRefn::Read: Can't open %s\n", binaryfile.c_str());
            return false;
       }

       Block &block = fobjR->GetBlock(fobjR->GetFirstBlockName());

       ISTable *t = getTablePtr(block, "ndb_refn");
       if (!t) {
            fprintf(log, "NdbRefn::Read: Can't find ndb_refn category in %s\n",
                         binaryfile.c_str());
            return false;
       }

       read_ndb_refn_table(t);
       create_indices();

       delete fobjR;

       return true;
}

std::map<std::string, std::string> NdbRefn::GetJournalInfo(const std::string& abbrev, const bool& is_epublished,
                                                           const std::string& publish_year)
{
       std::map<std::string, std::string> result_map;
       result_map.clear();

       if (abbrev.empty()) return result_map;

       std::string abbrev_key = get_abbrev_key(abbrev);
       if (_abbrev_mapping.find(abbrev_key) == _abbrev_mapping.end()) return result_map;

       std::pair<std::multimap<std::string, unsigned int>::const_iterator,
                 std::multimap<std::string, unsigned int>::const_iterator>
           range = _abbrev_mapping.equal_range(abbrev_key);
     
       std::vector<std::map<std::string, std::string> > results, essn_results, issn_results;
       results.clear();
       essn_results.clear();
       issn_results.clear();
       for (std::multimap<std::string, unsigned int>::const_iterator
            pos = range.first; pos != range.second; ++pos) {
            results.push_back(_all_refs[pos->second]);
            if (_all_refs[pos->second]["issn"] == "ESSN")
                 essn_results.push_back(_all_refs[pos->second]);
            else issn_results.push_back(_all_refs[pos->second]);
       }
       if (results.empty()) return result_map;

       if (is_epublished && essn_results.size() == 1) return essn_results[0];

       if (!is_epublished && issn_results.size() == 1) return issn_results[0];

       if (results.size() == 1) return results[0];

       if (is_epublished) {
            if (!essn_results.empty()) return __getJournalInfo(essn_results, publish_year);
       } else {
            if (!issn_results.empty()) return __getJournalInfo(issn_results, publish_year);
       }

       return __getJournalInfo(results, publish_year);
}

std::map<std::string, std::string> NdbRefn::GetJournalInfo(const std::string& abbrev, const std::string& issn,
                                                           const std::string& publish_year)
{
       std::map<std::string, std::string> result_map;
       result_map.clear();

       if (issn.empty()) return GetJournalInfo(abbrev, true, publish_year);

       if (_issn_mapping.find(issn) == _issn_mapping.end()) return GetJournalInfo(abbrev, true, publish_year);

       std::pair<std::multimap<std::string, unsigned int>::const_iterator,
                 std::multimap<std::string, unsigned int>::const_iterator>
           range = _issn_mapping.equal_range(issn);
     
       std::string abbrev_key = get_abbrev_key(abbrev);

       std::vector<std::map<std::string, std::string> > results, key_results;
       results.clear();
       key_results.clear();
       for (std::multimap<std::string, unsigned int>::const_iterator
            pos = range.first; pos != range.second; ++pos) {
            results.push_back(_all_refs[pos->second]);
            std::string key = get_abbrev_key(_all_refs[pos->second]["publication"]);
            if (key == abbrev_key) key_results.push_back(_all_refs[pos->second]);
       }
       if (results.empty()) return result_map;

       if (key_results.size() == 1) return key_results[0];

       if (results.size() == 1) return results[0];

       if (!key_results.empty()) return __getJournalInfo(key_results, publish_year);

       return __getJournalInfo(results, publish_year);
}

void NdbRefn::read_ndb_refn_table(ISTable* t)
{
       std::map<std::string, std::string> data_mapping;
       std::string cs;

       const std::vector<std::string>& columnnames = t->GetColumnNames();
       unsigned int rowNo = t->GetNumRows();
       for (unsigned int i = 0; i < rowNo; ++i) {
            data_mapping.clear();
            for (std::vector<std::string>::const_iterator
                 pos = columnnames.begin(); pos != columnnames.end(); ++pos) {
                 get_value_clean(cs, t, i, *pos);
                 if (cs.empty()) continue;
                 data_mapping.insert(std::make_pair(*pos, cs));
            }
            if (data_mapping.empty()) continue;
            if (data_mapping.find("publication") == data_mapping.end()) continue;
            _all_refs.push_back(data_mapping);
       }
}

void NdbRefn::create_indices()
{
       for (unsigned int i = 0; i < _all_refs.size(); ++i) {
            std::string abbrev_key = get_abbrev_key(_all_refs[i]["publication"]);
            _abbrev_mapping.insert(std::make_pair(abbrev_key, i));
            set_string_to_mixed(_all_refs[i]["publication"]);
            std::map<std::string, std::string>::const_iterator
                mpos = _all_refs[i].find("issn_code");
            if (mpos == _all_refs[i].end()) continue;
            _issn_mapping.insert(std::make_pair(mpos->second, i));
       }
}

std::string NdbRefn::get_abbrev_key(const std::string& abbrev)
{
       std::string abbrev_key;
       String::UpperCase(abbrev, abbrev_key);
       RemoveNonAlnum(abbrev_key);
       return abbrev_key;
}

void NdbRefn::set_string_to_mixed(std::string& str)
{
       String::LowerCase(str);
       str[0] = toupper(str[0]);
       for (unsigned int i = 1; i < str.size(); ++i) {
            if ((str[i - 1] < 'a' || str[i - 1] > 'z') &&
                (str[i - 1] < 'A' || str[i - 1] > 'Z')) {
                 str[i] = toupper(str[i]);
            }
       }
       if (str == "Proc.Natl.Acad.Sci.Usa") str = "Proc.Natl.Acad.Sci.USA";
}

std::map<std::string, std::string> NdbRefn::__getJournalInfo(const std::vector<std::map<std::string, std::string> >& results,
                                                             const std::string& publish_year)
{
       if (publish_year.empty()) return results[0];
 
       for (std::vector<std::map<std::string, std::string> >::const_iterator
            pos = results.begin(); pos != results.end(); ++pos) {
            std::map<std::string, std::string>::const_iterator
                mpos = pos->find("end");
            if (mpos != pos->end()) {
                 if (publish_year <= mpos->second) return *pos;
            }
            mpos = pos->find("start");
            if (mpos != pos->end()) {
                 if (publish_year >= mpos->second) return *pos;
            }
       }
       return results[0];
}
