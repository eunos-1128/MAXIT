/*
FILE:     StandardUtil.C
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
#include <sys/stat.h>
#include <math.h>

#include "StandardUtil.h"
#include "StandardUtil_global.h"
#include "utillib.h"

#define SIGMA   6

static void add_values(ISTable* bond, int&row, const std::string& name,
                       const std::vector<_MEAN_STD>& values,
                       const std::vector<std::string>& items, const std::string& val,
                       const std::string& esd, const int& precision);

void StandardUtil::initialize()
{
       _sigma_value = SIGMA;
       _standard_bonds.clear();
       _standard_angles.clear();
       _empty_standard.clear();
}

void StandardUtil::setSigma(const std::string& sigma)
{
       if (sigma.empty() || !String::IsNumber(sigma)) return;
       _sigma_value = atof(sigma.c_str());
}

bool StandardUtil::Read(LogUtil& logIo, const std::string& path)
{
       std::string ciffile = path + "/data/ascii/bond_angle_values.cif";
       struct stat statbuf;
       if (stat(ciffile.c_str(), &statbuf) != 0) {
            std::string error = "StandardUtil::Read: Can't find "
                              + ciffile + "\n";
            logIo.message(error.c_str());
            return false;
       }

       std::string /* odbfile, */ cs;
       // get_temp_filename(odbfile, "/tmp");
       CifFile *fobjR = read_cif_file("", ciffile, cs);
       if (!cs.empty()) {
            std::string error = "StandardUtil::Read: Read file " + ciffile
                    + "  error:\n" + cs + "\n";
            logIo.message(error.c_str());
            if (fobjR) delete fobjR;
            // remove(odbfile.c_str());
            return false;
       }

       Block &block = fobjR->GetBlock(fobjR->GetFirstBlockName());

       ISTable *Table = getTablePtr(block, "chem_comp_bond");
       if (!Table) {
            std::string error =
               "StandardUtil::Read: Can't find chem_comp_bond category in "
                              + ciffile + "\n";
            logIo.message(error.c_str());
            delete fobjR;
            // remove(odbfile.c_str());
            return false;
       }

       _MEAN_STD mean_std;
       std::vector<_MEAN_STD> m_vector;

       int rowNo = Table->GetNumRows();
       std::string id, val, esd, type;
       vector<string> data, items; 
       items.clear();
       items.push_back("atom_id_1");
       items.push_back("atom_id_2");
       items.push_back("atom_id_3");
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(id, Table, i, "comp_id");
            get_value_clean(val, Table, i, "value_dist");
            get_value_clean(esd, Table, i, "value_dist_esd");
            get_value_clean(type, Table, i, "type");
            if (type.empty()) type = ".";
            data.clear();
            for (int j = 0; j < 2; ++j) {
                 get_value_clean(cs, Table, i, items[j]);
                 data.push_back(cs);
            }

            std::map<std::string, std::vector<_MEAN_STD> >::iterator
                mpos = _standard_bonds.find(id);
            if (mpos == _standard_bonds.end()) {
                 mean_std.ats = data;
                 mean_std.mean_std.clear();
                 mean_std.mean_std.push_back(std::make_pair(std::make_pair(atof(val.c_str()),
                                             atof(esd.c_str())), type));
                 m_vector.clear();
                 m_vector.push_back(mean_std);
                 _standard_bonds.insert(std::make_pair(id, m_vector));
            } else {
                 unsigned int j = mpos->second.size() - 1;
                 unsigned int count = 0;
                 for (unsigned int k = 0; k < data.size(); ++k) {
                      if (data[k] == mpos->second[j].ats[k]) count++;
                 }
                 if (count == data.size()) {
                      mpos->second[j].mean_std.push_back(std::make_pair(std::make_pair(atof(val.c_str()),
                                                         atof(esd.c_str())), type));
                 } else {
                      mean_std.ats = data;
                      mean_std.mean_std.clear();
                      mean_std.mean_std.push_back(std::make_pair(std::make_pair(atof(val.c_str()),
                                                  atof(esd.c_str())), type));
                      mpos->second.push_back(mean_std);
                 }
            }
       }

       Table = getTablePtr(block, "chem_comp_angle");
       if (!Table) {
            std::string error =
               "StandardUtil::Read: Can't find chem_comp_angle category in "
                              + ciffile + "\n";
            logIo.message(error.c_str());
            delete fobjR;
            // remove(odbfile.c_str());
            initialize();
            return false;
       }

       rowNo = Table->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(id, Table, i, "comp_id");
            get_value_clean(val, Table, i, "value_angle");
            get_value_clean(esd, Table, i, "value_angle_esd");
            get_value_clean(type, Table, i, "type");
            if (type.empty()) type = ".";
            data.clear();
            for (int j = 0; j < 3; ++j) {
                 get_value_clean(cs, Table, i, items[j]);
                 data.push_back(cs);
            }

            std::map<std::string, std::vector<_MEAN_STD> >::iterator
                mpos = _standard_angles.find(id);
            if (mpos == _standard_angles.end()) {
                 mean_std.ats = data;
                 mean_std.mean_std.clear();
                 mean_std.mean_std.push_back(std::make_pair(std::make_pair(atof(val.c_str()),
                                             atof(esd.c_str())), type));
                 m_vector.clear();
                 m_vector.push_back(mean_std);
                 _standard_angles.insert(std::make_pair(id, m_vector));
            } else {
                 unsigned int j = mpos->second.size() - 1;
                 unsigned int count = 0;
                 for (unsigned int k = 0; k < data.size(); ++k) {
                      if (data[k] == mpos->second[j].ats[k]) count++;
                 }
                 if (count == data.size()) {
                      mpos->second[j].mean_std.push_back(std::make_pair(std::make_pair(atof(val.c_str()),
                                                         atof(esd.c_str())), type));
                 } else {
                      mean_std.ats = data;
                      mean_std.mean_std.clear();
                      mean_std.mean_std.push_back(std::make_pair(std::make_pair(atof(val.c_str()),
                                                                 atof(esd.c_str())), type));
                      mpos->second.push_back(mean_std);
                 }
            }
       }

       delete fobjR;
       // remove(odbfile.c_str());

       return true;
}

const std::vector<_MEAN_STD>& StandardUtil::GetBond(const std::string& ResName)
{
       std::map<std::string, std::vector<_MEAN_STD> >::const_iterator
           mpos = _standard_bonds.find(ResName);
       if (mpos != _standard_bonds.end())
            return mpos->second;
       else return _empty_standard;
}

const std::vector<_MEAN_STD>& StandardUtil::GetAngle(const std::string& ResName)
{
       std::map<std::string, std::vector<_MEAN_STD> >::const_iterator
           mpos = _standard_angles.find(ResName);
       if (mpos != _standard_angles.end())
            return mpos->second;
       else return _empty_standard;
}

bool StandardUtil::isOutlier(const double& val, const std::string& type, double& ept,
                             double& std, const std::vector<std::pair<std::pair<float,
                             float>, std::string> >& mean_std)
{
       _getValue(val, type, ept, std, mean_std);

       double range = _sigma_value * std;
       double high  = ept + range;
       double low   = ept - range;
       if (val < low || val > high) return true;
       return false;
}

void StandardUtil::_getValue(const double& val, const std::string& type, double& ept,
                             double& std, const std::vector<std::pair<std::pair<float,
                             float>, std::string> >& mean_std)
{
       ept = mean_std[0].first.first;
       std = mean_std[0].first.second;
       if (type == mean_std[0].second || mean_std.size() == 1) return;

       for (unsigned int i = 1; i < mean_std.size(); ++i) {
            if (type == mean_std[i].second) {
                 ept = mean_std[i].first.first;
                 std = mean_std[i].first.second;
                 return;
            }
       }

       _getClosestValue(val, ept, std, mean_std);
}

void StandardUtil::_getClosestValue(const double& val, double& ept, double& std,
                                    const std::vector<std::pair<std::pair<float,
                                    float>, std::string> >& mean_std)
{
       double diff = fabs(val - ept);
       for (unsigned int i = 1; i < mean_std.size(); ++i) {
            if (fabs(val - mean_std[i].first.first) < diff) {
                 ept = mean_std[i].first.first;
                 std = mean_std[i].first.second;
                 diff = fabs(val - ept);
            }
       }
}

void StandardUtil::Check(const std::string& filename)
{
       ISTable *bond = new ISTable("chem_comp_bond");
       bond->AddColumn("comp_id");
       bond->AddColumn("atom_id_1");
       bond->AddColumn("atom_id_2");
       bond->AddColumn("value_dist");
       bond->AddColumn("value_dist_esd");
       bond->AddColumn("type");

       ISTable *angle = new ISTable("chem_comp_angle");
       angle->AddColumn("comp_id");
       angle->AddColumn("atom_id_1");
       angle->AddColumn("atom_id_2");
       angle->AddColumn("atom_id_3");
       angle->AddColumn("value_angle");
       angle->AddColumn("value_angle_esd");
       angle->AddColumn("type");

       int bond_row = 0;
       int angle_row = 0;

       std::vector<std::string> items;
       items.clear();
       items.push_back("atom_id_1");
       items.push_back("atom_id_2");
       for (std::map<std::string, std::vector<_MEAN_STD> >::const_iterator
            mpos = _standard_bonds.begin(); mpos != _standard_bonds.end(); ++mpos) {
            add_values(bond, bond_row, mpos->first, mpos->second, items, "value_dist",
                       "value_dist_esd", 3);
       }

       items.clear();
       items.push_back("atom_id_1");
       items.push_back("atom_id_2");
       items.push_back("atom_id_3");
       for (std::map<std::string, std::vector<_MEAN_STD> >::const_iterator
            mpos = _standard_angles.begin(); mpos != _standard_angles.end(); ++mpos) {
            add_values(angle, angle_row, mpos->first, mpos->second, items, "value_angle",
                       "value_angle_esd", 1);
       }

       CifFile *fobj = create_fobj("", "STANDARD");
       Block& block = fobj->GetBlock("STANDARD");
       block.WriteTable(bond);
       block.WriteTable(angle);

       fobj->SetQuoting(CifFile::eDOUBLE);
       fobj->Write(filename);
       delete fobj;
}

static void add_values(ISTable* t, int&row, const std::string& name,
                       const std::vector<_MEAN_STD>& values,
                       const std::vector<std::string>& items, const std::string& val,
                       const std::string& esd, const int& precision)
{
       for (unsigned int i = 0; i < values.size(); ++i) {
            for (unsigned int j = 0; j < values[i].mean_std.size(); ++j) {
                 t->AddRow();
                 t->UpdateCell(row, "comp_id", name);
                 for (unsigned int k = 0; k < items.size(); ++k) {
                      t->UpdateCell(row, items[k], values[i].ats[k]);
                 }
                 t->UpdateCell(row, val, FloatToString(values[i].mean_std[j].first.first,
                                  0, precision, false, false));
                 t->UpdateCell(row, esd, FloatToString(values[i].mean_std[j].first.second,
                                  0, precision, false, false));
                 t->UpdateCell(row, "type", values[i].mean_std[j].second);
                 row++;
            }
       }
}
