/*
FILE:     RmsRecord.C
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
#include <stdexcept>

#include "RmsRecord.h"
#include "utillib.h"

#define rmsbinfile            "ndb_cif_rms.odb"

RmsRecord::RmsRecord()
{
       _programs.clear();
       _abbrevTypes.clear();
       _rmsRecords.clear();
}

RmsRecord::~RmsRecord()
{
       _programs.clear();
       _abbrevTypes.clear();
       _rmsRecords.clear();
}

const std::string RmsRecord::findAbbreviation(const std::string& name) const
{
       std::string cs;
       std::vector<std::string> data;

       String::LowerCase(name, cs);
       get_wordarray(data, cs, ",");

       for (std::vector<std::string>::const_iterator
            pos = data.begin(); pos != data.end(); ++pos) {
            std::map<std::string, std::string>::const_iterator
                mpos = _programs.find(*pos);
            if (mpos != _programs.end()) return mpos->second;

            for (mpos = _programs.begin(); mpos != _programs.end(); ++mpos) {
                 unsigned int size = pos->size();
                 if (mpos->first.size() < size) size = mpos->first.size();
                 if (pos->substr(0, size) == mpos->first.substr(0, size)) {
                      if (mpos->first == "refmac" && pos->substr(0, 8) == "refmac 5")
                           return "r";
                      else return mpos->second;
                 }
            }
       }

       return "o";
}

const std::vector<std::string> RmsRecord::findTypes(const std::string& abbrev) const
{
       std::vector<std::string> empty_vector;
       empty_vector.clear();

       std::map<std::string, std::vector<std::string> >::const_iterator
           mpos = _abbrevTypes.find(abbrev);
       if (mpos != _abbrevTypes.end()) return mpos->second;

       return empty_vector;
}

const std::map<std::string, std::pair<std::string, int> >
       RmsRecord::findItemTokenMapping(const std::string& type) const
{
       std::map<std::string, std::pair<std::string, int> > empty_map;
       empty_map.clear();

       std::map<std::string, std::map<std::string, std::pair<std::string, int> > >
          ::const_iterator mpos = _rmsRecords.find(type);
       if (mpos != _rmsRecords.end()) return mpos->second;

       return empty_map;
}

bool RmsRecord::Read(LogUtil& logIo, const std::string& path) 
{
       std::string binaryfile = path + "/data/binary/";
       binaryfile += rmsbinfile;

       struct stat statbuf;
       if (stat(binaryfile.c_str(), &statbuf) != 0) {
            std::string error = "RmsRecord::Read: Can't find " + binaryfile + "\n";
            logIo.message(error.c_str());
            return false;
       }

       CifFile* fobjR = new CifFile(READ_MODE, binaryfile);
       if (!fobjR) {
            std::string error = "RmsRecord::Read: Can't open " + binaryfile + "\n";
            logIo.message(error.c_str());
            return false;
       }

       Block &block = fobjR->GetBlock(fobjR->GetFirstBlockName());

       bool found_error = false;
       if (!block.IsTablePresent("ndb_cif_program")) {
            std::string error = "RmsRecord::Read: Can't find ndb_cif_program category in "
                              + binaryfile + "\n";
            logIo.message(error.c_str());
            found_error = true;
       }
       if (!block.IsTablePresent("ndb_cif_rms")) {
            std::string error = "RmsRecord::Read: Can't find ndb_cif_rms category in "
                              + binaryfile + "\n";
            logIo.message(error.c_str());
            found_error = true;
       }
       if (found_error) {
            delete fobjR; return false;
       }

       std::string cs, cs1;

       ISTable *t = getTablePtr(block, "ndb_cif_program");

       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value(cs, t, i, "full_name");
            get_value(cs1, t, i, "abbreviation");
            _programs.insert(std::make_pair(cs, cs1));
       }

       t = getTablePtr(block, "ndb_cif_rms");

       
       std::map<std::string, std::pair<std::string, std::string> > read_item_mapping;
       read_item_mapping.clear();
       read_item_mapping.insert(std::make_pair("dev_ideal",
                    std::make_pair("dev_ideal_token", "dev_ideal_field")));
       read_item_mapping.insert(std::make_pair("dev_ideal_target",
                    std::make_pair("dev_ideal_target_token", "dev_ideal_target_field")));
       read_item_mapping.insert(std::make_pair("weight",
                    std::make_pair("ndb_weight_token", "ndb_weight_field")));
       read_item_mapping.insert(std::make_pair("number",
                    std::make_pair("ndb_count_token", "ndb_count_field")));
       read_item_mapping.insert(std::make_pair("pdbx_restraint_function",
                    std::make_pair("ndb_function_token", "ndb_function_field")));

       std::vector<std::string> tmp_vector;

       std::map<std::string, std::pair<std::string, int> > write_item_mapping;

       std::string type;
       rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(type, t, i, "type");

            write_item_mapping.clear();
            for (std::map<std::string, std::pair<std::string, std::string> >::const_iterator
                 impos = read_item_mapping.begin(); impos != read_item_mapping.end(); ++impos) {
                 get_value_clean(cs, t, i, impos->second.first);
                 if (cs.empty()) continue;
                 get_value_clean(cs1, t, i, impos->second.second);
                 if (cs1.empty()) continue;
                 int field = atoi(cs1.c_str()) - 1;
                 write_item_mapping.insert(std::make_pair(impos->first,
                                           std::make_pair(cs, field)));
            }
            _rmsRecords.insert(std::make_pair(type, write_item_mapping));

            std::string abbrev = type.substr(0, 1);
            std::map<std::string, std::vector<std::string> >::iterator
                mpos = _abbrevTypes.find(abbrev);
            if (mpos == _abbrevTypes.end()) {
                 tmp_vector.clear();
                 tmp_vector.push_back(type);
                 _abbrevTypes.insert(std::make_pair(abbrev, tmp_vector));
            } else mpos->second.push_back(type);
       }

       if (fobjR) delete fobjR;

       return true;
}
