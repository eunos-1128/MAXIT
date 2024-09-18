/*
FILE:     CifUpdate_Special_Util.C
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

#include "AnnotationObj.h"
#include "utillib.h"

void AnnotationObj::_cif_update_diffrn_source(Block& block)
{
       ISTable *t = _getTablePtr(block, "diffrn_source");
       if (!t) {
            deleteTable(block, "diffrn_radiation_wavelength");
            return;
       }

       std::set<std::string> unique_data_set;
       unique_data_set.clear();

       std::vector<std::string> data, data_set;
       std::string cs, cs1, cs2;

       data_set.clear();

       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean_upper(cs, t, i, "source");
            if (cs == "SYNCHROTRON") t->UpdateCell(i, "pdbx_synchrotron_y_n", "Y");

            get_value_clean_upper(cs, t, i, "type");
            if (!cs.empty()) {
                 std::string::size_type p = cs.find(" BEAMLINE ");
                 if (p != std::string::npos) {
                      get_value_clean(cs, t, i, "type");
                      t->UpdateCell(i, "pdbx_synchrotron_site", cs.substr(0, p));
                      t->UpdateCell(i, "pdbx_synchrotron_beamline", cs.substr(p + 10));
                 }
            }

            get_value_clean(cs1, t, i, "pdbx_wavelength");
            get_value_clean(cs2, t, i, "pdbx_wavelength_list");
            if (cs1.empty() && cs2.empty()) continue;

            cs = cs1 + " " + cs2;
            get_wordarray(data, cs, ",;- ");
            for (std::vector<std::string>::const_iterator
                 pos = data.begin(); pos != data.end(); ++pos) {
                 if (unique_data_set.find(*pos) != unique_data_set.end()) continue;
                 unique_data_set.insert(*pos);
                 data_set.push_back(*pos);
            }
       }
       block.WriteTable(t);

       if (data_set.empty()) {
            deleteTable(block, "diffrn_radiation_wavelength");
            return;
       }

       ISTable* t1 = _newTablePtr("diffrn_radiation_wavelength");
       for (std::vector<std::string>::const_iterator pos = data_set.begin(); pos != data_set.end(); ++pos) {
            bool not_a_number = false;
            cs = *pos;
            for (unsigned int j = 0; j < cs.size(); ++j) {
                 if (cs[j] != '.' && !isdigit(cs[j])) {
                      not_a_number = true;
                      break;
                 }
            }
            if (not_a_number) continue;

            int irow = t1->GetNumRows();
            t1->AddRow();
            t1->UpdateCell(irow, "id", String::IntToString(irow + 1));
            t1->UpdateCell(irow, "wavelength", *pos);
            t1->UpdateCell(irow, "wt", "1.0");
       }
       block.WriteTable(t1);
}

void AnnotationObj::_cif_update_diffrn_radiation(Block& block)
{
       if (_experiment_type & EXPERIMENT_TYPE_XRAY)     _experimental_method_count++;
       if (_experiment_type & EXPERIMENT_TYPE_NEUTRON)  _experimental_method_count++;
       if (_experiment_type & EXPERIMENT_TYPE_FIBER)    _experimental_method_count++;
       if (_experiment_type & EXPERIMENT_TYPE_ELECTRON) _experimental_method_count++;
       if (_experiment_type & EXPERIMENT_TYPE_EC) _experimental_method_count++;
       // if (_experiment_type & EXPERIMENT_TYPE_EM)  _experimental_method_count++;
       // if (_experiment_type & EXPERIMENT_TYPE_SOLN_SCT) _experimental_method_count++;

       ISTable *t = _getTablePtr(block, "diffrn_radiation");
       if (!t) {
            if (_experimental_method_count) _default_diffrn_id = "1";
            return;
       }

       std::string diffrn_id, cs, cs1;

       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(diffrn_id, t, i, "diffrn_id");
            if (diffrn_id.empty()) {
                 diffrn_id = String::IntToString(i + 1);
                 t->UpdateCell(i, "diffrn_id", diffrn_id);
            }
            if (!_default_diffrn_id.empty()) _default_diffrn_id += ",";
            _default_diffrn_id += diffrn_id;

            get_value_clean(cs, t, i, "wavelength_id");
            if (cs.empty()) t->UpdateCell(i, "wavelength_id", "1"); 

            get_value_clean_upper(cs, t, i, "pdbx_monochromatic_or_laue_m_l");
            if (cs == "N")
                 t->UpdateCell(i, "pdbx_monochromatic_or_laue_m_l", "M");
            else if (cs == "Y")
                 t->UpdateCell(i, "pdbx_monochromatic_or_laue_m_l", "L");

            get_value_clean(cs, t, i, "pdbx_scattering_type");
            if (!cs.empty()) continue;
            get_value_clean_upper(cs, t, i, "pdbx_diffrn_protocol");
            cs1.clear();
            if (cs == "SINGLE WAVELENGTH" || cs == "MAD")
                 cs1 = "x-ray";
            else if (cs == "LAUE")
                 cs1 = "neutron";
            if (cs1.empty()) {
                 if (_experiment_type == EXPERIMENT_TYPE_XRAY)
                      cs1 = "x-ray";
                 else if (_experiment_type == EXPERIMENT_TYPE_NEUTRON)
                      cs1 = "neutron";
                 else if ((_experiment_type == EXPERIMENT_TYPE_ELECTRON) || (_experiment_type == EXPERIMENT_TYPE_EC) ||
                          (_experiment_type == EXPERIMENT_TYPE_EM) || (_experiment_type == EXPERIMENT_TYPE_ET))
                      cs1 = "electron";
            }
            if (!cs1.empty()) {
                 t->UpdateCell(i, "pdbx_scattering_type", cs1);
                 if (cs1 == "x-ray") _scattering_type_mapping.insert(std::make_pair("X-RAY DIFFRACTION", diffrn_id));
                 else if (cs1 == "neutron") _scattering_type_mapping.insert(std::make_pair("NEUTRON DIFFRACTION", diffrn_id));
                 else if (cs1 == "electron") {
                      if (_experiment_type == EXPERIMENT_TYPE_ELECTRON)
                           _scattering_type_mapping.insert(std::make_pair("ELECTRON DIFFRACTION", diffrn_id));
                      else if (_experiment_type & EXPERIMENT_TYPE_EC)
                           _scattering_type_mapping.insert(std::make_pair("ELECTRON CRYSTALLOGRAPHY", diffrn_id));
                      else if (_experiment_type == EXPERIMENT_TYPE_EM)
                           _scattering_type_mapping.insert(std::make_pair("ELECTRON MICROSCOPY", diffrn_id));
                      else if (_experiment_type == EXPERIMENT_TYPE_ET)
                           _scattering_type_mapping.insert(std::make_pair("ELECTRON TOMOGRAPHY", diffrn_id));
                 }
            }
       }
       block.WriteTable(t);
}

void AnnotationObj::_cif_update_refine_ls_restr_ncs(Block &block)
{
       ISTable *t = _getTablePtr(block, "refine_ls_restr_ncs");
       if (!t) return;

       std::set<std::string> keyitems;
       keyitems.clear();
       keyitems.insert("dom_id");
       keyitems.insert("pdbx_ens_id");
       keyitems.insert("pdbx_refine_id");
       keyitems.insert("pdbx_ordinal");

       const std::vector<std::string>& ColumnNames = t->GetColumnNames();
       unsigned int rowNo = t->GetNumRows();
       std::vector<unsigned int> rows;
       rows.clear();
       for (unsigned int i = 0; i < rowNo; ++i) {
            if (is_empty_row(t, ColumnNames, i, keyitems)) {
                 rows.push_back(i);
            }
       }
       if (!rows.empty()) t->DeleteRows(rows);

       std::string cs, cs1;

       rowNo = t->GetNumRows();
       if (rowNo == 0) {
            deleteTable(block, "refine_ls_restr_ncs");
            return;
       }

       for (unsigned int i = 0; i < rowNo; ++i) {
            t->UpdateCell(i, "pdbx_ordinal", String::IntToString(i + 1));
            get_value_clean(cs, t, i, "pdbx_ens_id");
            get_value_clean(cs1, t, i, "dom_id");
            if (cs.empty() && !cs1.empty()) t->UpdateCell(i, "pdbx_ens_id", cs1);
       }

       if (rowNo == 1) {
            get_value_clean(cs, t, 0, "dom_id");
            if (cs.empty()) t->UpdateCell(0, "dom_id", "1");
            get_value_clean(cs, t, 0, "pdbx_ens_id");
            if (cs.empty()) t->UpdateCell(0, "pdbx_ens_id", "1");
       }

       block.WriteTable(t);
}

void AnnotationObj::_cif_update_refine(Block &block)
{
       ISTable *t = _getTablePtr(block, "refine");
       if (!t) return;

       std::string cs;

       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value(cs, t, i, "pdbx_refine_id");
            if (cs.empty() || cs == "1") t->UpdateCell(i, "pdbx_refine_id", _default_exp_type);
            get_value(cs, t, i, "pdbx_ls_sigma_F");
            std::string::size_type p = cs.find("COMPLETENESS FOR RANGE (%) :");
            if (p != std::string::npos) {
                 std::string cs1 = cs.substr(0, p);
                 String::StripAndCompressWs(cs1);    
                 t->UpdateCell(i, "pdbx_ls_sigma_F", cs1);
                 cs1 = cs.substr(p);
                 String::StripAndCompressWs(cs1);
                 if (!cs1.empty() && cs1 != "NULL" && cs1 != "N/A") {
                      get_value(cs, t, i, "ls_percent_reflns_obs");
                      if (cs.empty()) t->UpdateCell(i, "ls_percent_reflns_obs", cs1);
                 }
            }

            get_value(cs, t, i, "B_iso_mean");
            p = cs.find("OVERALL ANISOTROPIC B VALUE");
            if (p != std::string::npos) {
                 cs.erase(p);
                 String::StripAndCompressWs(cs);
                 t->UpdateCell(i, "B_iso_mean", cs);
            }

            get_value(cs, t, i, "aniso_B[2][3]");
            p = cs.find("ESTIM");
            if (p != std::string::npos) {
                 cs.erase(p);
                 String::StripAndCompressWs(cs);
                 t->UpdateCell(i, "aniso_B[2][3]", cs);
            }

            _remove_infinite(t, "ls_d_res_low", i);
            _remove_comma(t, "ls_number_reflns_obs", i);
            _remove_comma(t, "ls_number_reflns_all", i);
       }
       block.WriteTable(t);

       if (!(_experiment_type & EXPERIMENT_TYPE_XRAY)) return;

       ISTable *t1 = _getTablePtr(block, "refine_hist");
       if (!t1) {
            t1 = _newTablePtr("refine_hist");
            t1->AddRow();
       }

       get_value_clean(cs, t, 0, "ls_d_res_high");
       if (!cs.empty()) t1->UpdateCell(0, "d_res_high", cs);

       get_value_clean(cs, t, 0, "ls_d_res_low");
       if (!cs.empty()) t1->UpdateCell(0, "d_res_low", cs);

       get_value_clean(cs, t1, 0, "cycle_id");
       if (cs.empty()) t1->UpdateCell(0, "cycle_id", "LAST");

       get_value_clean(cs, t1, 0, "pdbx_refine_id");
       if (cs.empty()) t1->UpdateCell(0, "pdbx_refine_id", "X-RAY DIFFRACTION");

       int total_num = 0;

       get_value_clean(cs, t1, 0, "pdbx_number_atoms_protein");
       if (!cs.empty()) total_num += atoi(cs.c_str());
       if (cs.empty() && !_molecules.empty()) {
            int num = _molecules[0]->count_atom_no("ATOMP");
            total_num += num;
            t1->UpdateCell(0, "pdbx_number_atoms_protein", String::IntToString(num));
       }

       get_value_clean(cs, t1, 0, "pdbx_number_atoms_nucleic_acid");
       if (!cs.empty()) total_num += atoi(cs.c_str());
       if (cs.empty() && !_molecules.empty()) {
            int num = _molecules[0]->count_atom_no("ATOMN");
            total_num += num;
            t1->UpdateCell(0, "pdbx_number_atoms_nucleic_acid", String::IntToString(num));
       }

       get_value_clean(cs, t1, 0, "pdbx_number_atoms_ligand");
       if (!cs.empty()) total_num += atoi(cs.c_str());
       if (cs.empty() && !_molecules.empty()) {
            std::set<std::string> token_set;
            token_set.clear();
            token_set.insert("HETAD"); token_set.insert("HETAC"); token_set.insert("HETAI");
            token_set.insert("HETIC"); token_set.insert("HETAIN"); token_set.insert("ATOMS");
            int num = _molecules[0]->count_atom_no(token_set);
            total_num += num;
            t1->UpdateCell(0, "pdbx_number_atoms_ligand", String::IntToString(num));
       }

       get_value_clean(cs, t1, 0, "number_atoms_solvent");
       if (!cs.empty()) total_num += atoi(cs.c_str());
       if (cs.empty() && !_molecules.empty()) {
            int num = _molecules[0]->count_atom_no("HETAS");
            total_num += num;
            t1->UpdateCell(0, "number_atoms_solvent", String::IntToString(num));
       }

       t1->UpdateCell(0, "number_atoms_total", String::IntToString(total_num));

       block.WriteTable(t1);
}

void AnnotationObj::_remove_infinite(ISTable *t, const std::string &item, const int& row)
{
       std::string cs;
       get_value_upper(cs, t, row, item);
       if (cs.find("INFINI") != std::string::npos) t->UpdateCell(row, item, "");
}

void AnnotationObj::_remove_comma(ISTable *t, const std::string &item, const int& row)
{
       std::string cs;
       get_value(cs, t, row, item);
       std::string::size_type p = cs.find(",");
       if (p == std::string::npos) return;

       while (p != std::string::npos) {
            cs.erase(p, 1);
            p = cs.find(",");
       }
       t->UpdateCell(row, item, cs);
}
