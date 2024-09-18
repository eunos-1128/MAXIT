/*
FILE:     ConnectDic.C
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

/*!
** \file ConnectDic.C
**
** \brief Implementation file for ConnectDic class.
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "ConnectDic.h"
#include "ConnectDic_global.h"
#include "Exceptions.h"
#include "utillib.h"

ConnectDic::ConnectDic()
{
       clear();
       init();
}

ConnectDic::~ConnectDic()
{
       if (_fObj) delete _fObj;
       clear();
}

void ConnectDic::clear()
{
       _empty_string.clear();
       _empty_set.clear();
       _empty_vector.clear();
       _drugs.clear();
       _ConnMap.clear();
       _VariantConnMap.clear();
       _fObj = NULL;
       _rcsbrootpath.clear();
       _cvscomponentpath.clear();
       _additional_comp_path.clear();
       _peptide_type.clear();
       _nucleic_type.clear();
       _dna_type.clear();
       _rna_type.clear();
       _atom_type.clear();
       _tailTerminalAtoms.clear();
       _headTerminalAtoms.clear();
       _userProvidedCCDFiles.clear();
}

void ConnectDic::init()
{
       for (int i = 0; i < NUM_PEPTIDE_TYPE; ++i)
            _peptide_type.insert(_peptide_tokens[i]);

       for (int i = 0; i < NUM_DNA_TYPE; ++i) {
            _nucleic_type.insert(_dna_tokens[i]);
            _dna_type.insert(_dna_tokens[i]);
       }

       for (int i = 0; i < NUM_RNA_TYPE; ++i) {
            _nucleic_type.insert(_rna_tokens[i]);
            _rna_type.insert(_rna_tokens[i]);
       }

       for (int i = 0; i < NUM_ATOM_TYPE; ++i)
            _atom_type.insert(_atom_tokens[i]);
}

void ConnectDic::add_drug(Block &block, const bool& with_preferred_name)
{
       ConnectFormat drug;
       drug.clear();
       if (!drug.Read(block, true, with_preferred_name)) return;

       _ConnMap.insert(std::make_pair(drug.drugname(), _drugs.size()));
       _VariantConnMap.insert(std::make_pair(drug.getMetaData("three_letter_code"), _drugs.size()));
       _drugs.push_back(drug);
}

void ConnectDic::add_drug(const std::string& drugID, const std::string& compfile, const bool& with_preferred_name)
{
       struct stat statbuf;
       if (stat(compfile.c_str(), &statbuf) != 0) return;

       std::string cs;
       CifFile *fobj = read_cif_file("", compfile, cs);
       if (cs.empty() && fobj) {
            if (fobj->IsBlockPresent(drugID)) {
                 Block &block = fobj->GetBlock(drugID);
                 add_drug(block, with_preferred_name);
            }
       }
       if (fobj) delete fobj;
}

std::set<std::string> ConnectDic::_get_terminal_atoms(const std::string& drugID, const std::string& searchType, const std::string& excludeType)
{
       std::set<std::string> terminalAtoms;
       terminalAtoms.clear();

       const std::vector<std::string>& sAtoms = find_terminal_atoms(drugID, searchType);
       if (!sAtoms.empty()) {
            for (std::vector<std::string>::const_iterator vpos = sAtoms.begin(); vpos != sAtoms.end(); ++vpos) terminalAtoms.insert(*vpos);
            return terminalAtoms;
       }

       if (terminalAtoms.empty()) {
            terminalAtoms = get_connected_atoms(drugID);
            if (!terminalAtoms.empty()) {
                 const std::vector<std::string>& eAtoms = find_terminal_atoms(drugID, excludeType);
                 if (!eAtoms.empty()) {
                      for (std::vector<std::string>::const_iterator vpos = eAtoms.begin(); vpos != eAtoms.end(); ++vpos) terminalAtoms.erase(*vpos);
                 }
            }
       }
       return terminalAtoms;
}

void ConnectDic::setRCSBROOT(const std::string& path)
{
       if (path.empty()) return;

       struct stat statbin;
       _rcsbrootpath = path;
       std::string additional_comp_path = path;
       if (_rcsbrootpath[_rcsbrootpath.size() - 1] != '/') {
            _rcsbrootpath += "/data/binary/";
            additional_comp_path += "/data/ccd";
       } else {
            _rcsbrootpath += "data/binary/";
            additional_comp_path += "data/ccd";
       }

       if (stat(additional_comp_path.c_str(), &statbin) == 0) setAdditionalCompPath(additional_comp_path);

       if (stat(_rcsbrootpath.c_str(), &statbin) != 0) {
            _rcsbrootpath.clear();
            return;
       }
       std::string filename = _rcsbrootpath + componentbinaryfile;
       if (stat(filename.c_str(), &statbin) != 0) {
            _rcsbrootpath.clear();
       }
}

void ConnectDic::setCOMP_PATH(const std::string& path)
{
       if (path.empty()) return;

       struct stat statbin;
       _cvscomponentpath = path;
       if (_cvscomponentpath[_cvscomponentpath.size() - 1] != '/')
            _cvscomponentpath += "/";
       if (stat(_cvscomponentpath.c_str(), &statbin) != 0) {
            _cvscomponentpath.clear();
            return;
       }
       for (int i = 0; i < NUM_HASH; ++i) {
            std::string hashpath = _cvscomponentpath + _hash_path[i];
            if (stat(hashpath.c_str(), &statbin) != 0) {
                 _cvscomponentpath.clear();
                 return;
            }
       }
}

void ConnectDic::setAdditionalCompPath(const std::string& path)
{
       if (path.empty()) return;

       struct stat statbin;
       _additional_comp_path = path;
       if (_additional_comp_path[_additional_comp_path.size() - 1] != '/')
            _additional_comp_path += "/";
       if (stat(_additional_comp_path.c_str(), &statbin) != 0) {
            _additional_comp_path.clear();
            return;
       }
       for (int i = 0; i < NUM_HASH; ++i) {
            std::string hashpath = _additional_comp_path + _hash_path[i];
            if (stat(hashpath.c_str(), &statbin) != 0) {
                 _additional_comp_path.clear();
                 return;
            }
       }
}

const std::string& ConnectDic::getCOMP_PATH() const
{
       return _cvscomponentpath;
}

const std::string& ConnectDic::getAdditionalCompPath() const
{
       return _additional_comp_path;
}

void ConnectDic::readUserProvidedCCDictionary(const std::string& cc_directory_path)
{
       struct stat statbin;
       if (stat(cc_directory_path.c_str(), &statbin) != 0) return;

       std::vector<std::string> file_list, data;
       get_file_list_from_directory(file_list, cc_directory_path);
       if (file_list.empty()) return;

       for (std::vector<std::string>::const_iterator fpos = file_list.begin(); fpos != file_list.end(); ++fpos) {
            std::string filename = cc_directory_path + "/" + *fpos;
            if (stat(filename.c_str(), &statbin) != 0) continue;

            get_wordarray(data, *fpos, ".");
            if ((data.size() != 2) || (data[1] != "cif")) continue;

            _userProvidedCCDFiles.insert(std::make_pair(data[0], filename));
       }
}

void ConnectDic::readUserProvidedCCDictionaryAll(const std::string& cc_directory_path)
{
       struct stat statbin;
       if (stat(cc_directory_path.c_str(), &statbin) != 0) return;

       std::vector<std::string> file_list, blockNames;
       get_file_list_from_directory(file_list, cc_directory_path);
       if (file_list.empty()) return;

       for (std::vector<std::string>::const_iterator fpos = file_list.begin(); fpos != file_list.end(); ++fpos) {
            std::string filename = cc_directory_path + "/" + *fpos;
            if (stat(filename.c_str(), &statbin) != 0) continue;

            CifFile *fobj = get_fobj(stdout, "", filename);
            if (!fobj) continue;

            fobj->GetBlockNames(blockNames);
            for (std::vector<std::string>::const_iterator bpos = blockNames.begin(); bpos != blockNames.end(); ++bpos) {
                 Block& block = fobj->GetBlock(*bpos);
                 add_drug(block);
            }
            delete fobj;
       }
}

bool ConnectDic::OpenFile()
{
       if (_rcsbrootpath.empty() && _cvscomponentpath.empty()) return false;

       if (!_rcsbrootpath.empty()) {
            std::string filename = _rcsbrootpath + variantbinaryfile;
            struct stat statbin;
            if (stat(filename.c_str(), &statbin) == 0) {
                 CifFile *fobj = new CifFile(READ_MODE, filename);
                 if (fobj) {
                      std::vector<std::string> blockNames;
                      fobj->GetBlockNames(blockNames);
                      for (unsigned int i = 0; i < blockNames.size(); ++i) {
                           if (blockNames[i] == "INDEX" ||
                               blockNames[i] == "OBSOLETE") continue;
               
                           Block &block = fobj->GetBlock(blockNames[i]);
                           add_drug(block);
                      }
                      delete fobj;
                 }
            }

            filename = _rcsbrootpath + componentbinaryfile;
            if (stat(filename.c_str(), &statbin) == 0) {
                 _fObj = new CifFile(READ_MODE, filename);
                 if (!_fObj) return false;
            }
       }

       return true;
}

const ConnectFormat& ConnectDic::find_drug(const std::string& drugID, const std::string& compfile)
{
       if (drugID.empty() && compfile.empty()) throw NotFoundException();

       std::map<std::string, int>::const_iterator pos = _ConnMap.find(drugID);
       if (pos != _ConnMap.end()) return _drugs[pos->second];

       if (!compfile.empty()) {
            add_drug(drugID, compfile);
            pos = _ConnMap.find(drugID);
            if (pos != _ConnMap.end())
                 return _drugs[pos->second];
            else throw NotFoundException();
       }

       std::map<std::string, std::string>::const_iterator mpos = _userProvidedCCDFiles.find(drugID);
       if (mpos != _userProvidedCCDFiles.end()) {
            add_drug(drugID, mpos->second, true);
            pos = _ConnMap.find(drugID);
            if (pos != _ConnMap.end()) return _drugs[pos->second];
       }

       std::string hashID = drugID.substr(0, 1);
       if (drugID.size() > 3) hashID = drugID.substr(drugID.size() - 2);

       if (!_additional_comp_path.empty()) {
            std::string filename = _additional_comp_path + hashID + "/" + drugID + "/" + drugID + ".cif";
            add_drug(drugID, filename, true);
            pos = _ConnMap.find(drugID);
            if (pos != _ConnMap.end()) return _drugs[pos->second];
       }

       if (!_cvscomponentpath.empty()) {
            std::string filename = _cvscomponentpath + hashID + "/" + drugID + "/" + drugID + ".cif";
            add_drug(drugID, filename);
            pos = _ConnMap.find(drugID);
            if (pos != _ConnMap.end()) return _drugs[pos->second];
       }

       if (_fObj && _fObj->IsBlockPresent(drugID)) {
            Block &block = _fObj->GetBlock(drugID);
            add_drug(block);
            pos = _ConnMap.find(drugID);
            if (pos != _ConnMap.end()) return _drugs[pos->second];
       }

       throw NotFoundException();
}

const ConnectFormat& ConnectDic::find_author_defined_drug(const std::string& drugID)
{
       if (drugID.empty()) throw NotFoundException();

       std::map<std::string, int>::const_iterator pos = _ConnMap.find("Author_defined_" + drugID);
       if (pos != _ConnMap.end()) return _drugs[pos->second];

       throw NotFoundException();
}

void ConnectDic::add_author_defined_drug(const std::string& drugID, const ConnectFormat& drug)
{
       _ConnMap.insert(std::make_pair("Author_defined_" + drugID, _drugs.size()));
       _drugs.push_back(drug);
}

const std::string ConnectDic::find_residue_type(const std::string& residueName)
{
       std::string drugID, type;
       String::UpperCase(residueName, drugID);
       if (drugID == "UNL" || drugID == "UVL" || drugID == "UVR") return "HETAIN";

       try {
            const ConnectFormat& drug = find_drug(drugID);
            String::UpperCase(drug.getMetaData("type"), type);
            std::string AtomType = drug.getMetaData("pdbx_type");
            if (_peptide_type.find(type) != _peptide_type.end())
                 return "ATOMP";
            else if (_nucleic_type.find(type) != _nucleic_type.end())
                 return "ATOMN";
            else if (_atom_type.find(AtomType) != _atom_type.end())
                 return AtomType;
            return "HETAIN";
       } catch (const std::exception& exc) {
            return "unknown";
       }
}

const int ConnectDic::find_monomer_type(const std::string& residueName)
{
       std::string drugID, type;
       String::UpperCase(residueName, drugID);

       try {
            const ConnectFormat& drug = find_drug(drugID);
            String::UpperCase(drug.getMetaData("type"), type);
            if (_peptide_type.find(type) != _peptide_type.end())
                 return MONOMER_TYPE_A_A;
            else if (_dna_type.find(type) != _dna_type.end())
                 return MONOMER_TYPE_DNA;
            else if (_rna_type.find(type) != _rna_type.end())
                 return MONOMER_TYPE_RNA;

            return MONOMER_TYPE_UNK;

       } catch (const std::exception& exc) {
            return MONOMER_TYPE_UNK;
       }

}

const std::string ConnectDic::getAtomType(const std::string& drugID, const std::string& atomName)
{
       try { 
            const ConnectFormat& drug = find_drug(drugID);
            const AtomFormat& atom = drug.find_atom(atomName);
            return atom.atomtype();
       } catch (const std::exception& exc) {
            return "";
       }
}

const double ConnectDic::get_molecule_weight(const std::string& drugID)
{
       try { 
            const ConnectFormat& drug = find_drug(drugID);
/*
            if ((drugID == "ARG") || (drugID == "HIS") || (drugID == "LYS"))
                 return (atof(drug.getMetaData("formula_weight").c_str()) - 1.00794 * atof(drug.getMetaData("pdbx_formal_charge").c_str()));
            else */ return (atof(drug.getMetaData("formula_weight").c_str()));
       } catch (const std::exception& exc) {
            return 0.0;
       }
}

bool ConnectDic::is_L_aa_residue(const std::string& drugID)
{
       try { 
            const ConnectFormat& drug = find_drug(drugID);
            if (String::IsEqual(drug.getMetaData("type"), "L-peptide linking", Char::eCASE_INSENSITIVE)) return true;
            return false;
       } catch (const std::exception& exc) {
            return false;
       }
}

bool ConnectDic::is_D_aa_residue(const std::string& drugID)
{
       try { 
            const ConnectFormat& drug = find_drug(drugID);
            if (String::IsEqual(drug.getMetaData("type"), "D-peptide linking", Char::eCASE_INSENSITIVE)) return true;
            return false;
       } catch (const std::exception& exc) {
            return false;
       }
}

const std::string ConnectDic::find_cc_metadata(const std::string& residueName, const std::string& cifitem)
{
       std::string drugID;
       String::UpperCase(residueName, drugID);

       try {
            const ConnectFormat& drug = find_drug(drugID);
            return drug.getMetaData(cifitem);
       } catch (const std::exception& exc) {}
       return "";
}

const std::string& ConnectDic::find_terminal_atom(const std::string& residueName, const std::string& type)
{
       std::string drugID;
       String::UpperCase(residueName, drugID);

       try {
            const ConnectFormat& drug = find_drug(drugID);
            return drug.getFirstTerminalConnectedAtomByType(type);
       } catch (const std::exception& exc) {}
       return _empty_string;
}

const std::vector<std::string>& ConnectDic::find_terminal_atoms(const std::string& residueName, const std::string& type)
{
       std::string drugID;
       String::UpperCase(residueName, drugID);

       try {
            const ConnectFormat& drug = find_drug(drugID);
            return drug.getTerminalConnectedAtomByType(type);
       } catch (const std::exception& exc) {}
       return _empty_vector;
}

std::string ConnectDic::getTailTerminalAtoms(const std::string& type, const std::string& residueName, std::set<std::string>& terminalAtoms)
{
       std::string standardTerminalAtomName = "C";
       std::string standardTerminalAtomType = "C";
       std::string excludeType = "N";
       if (type == "ATOMN") {
            standardTerminalAtomName = "O3'";
            standardTerminalAtomType = "O";
            excludeType = "P";
       }

       std::map<std::string, std::set<std::string> >::const_iterator mpos = _tailTerminalAtoms.find(residueName);
       if (mpos != _tailTerminalAtoms.end()) terminalAtoms = mpos->second;
       else {
            terminalAtoms = _get_terminal_atoms(residueName, standardTerminalAtomType, excludeType);
            _tailTerminalAtoms.insert(std::make_pair(residueName, terminalAtoms));
       }

       return standardTerminalAtomName;
}

std::string ConnectDic::getHeadTerminalAtoms(const std::string& type, const std::string& residueName, std::set<std::string>& terminalAtoms)
{
       std::string standardTerminalAtomName = "N";
       std::string standardTerminalAtomType = "N";
       std::string excludeType = "C";
       if (type == "ATOMN") {
            standardTerminalAtomName = "P";
            standardTerminalAtomType = "P";
            excludeType = "O";
       }

       std::map<std::string, std::set<std::string> >::const_iterator mpos = _headTerminalAtoms.find(residueName);
       if (mpos != _headTerminalAtoms.end()) terminalAtoms = mpos->second;
       else {
            terminalAtoms = _get_terminal_atoms(residueName, standardTerminalAtomType, excludeType);
            _headTerminalAtoms.insert(std::make_pair(residueName, terminalAtoms));
       }

       return standardTerminalAtomName;
}

const bool ConnectDic::is_terminal_atom(const std::string& residueName, const std::string& atomName)
{
       std::string drugID;
       String::UpperCase(residueName, drugID);

       try {
            const ConnectFormat& drug = find_drug(drugID);
            return drug.isTerminalAtom(atomName);
       } catch (const std::exception& exc) {}
       return false;
}

void ConnectDic::find_drugs(const std::string &drugID, std::vector<ConnectFormat>& drug_list)
{
       drug_list.clear();

       try { 
            const ConnectFormat& drug = find_drug(drugID);
            drug_list.push_back(drug);

            std::pair<std::multimap<std::string, int>::iterator, std::multimap<std::string, int>::iterator> range = _VariantConnMap.equal_range(drugID);
            std::multimap<std::string, int>::iterator pos;
            for (pos = range.first; pos != range.second; ++pos) {
                 if (_drugs[pos->second].drugname() == drug.drugname()) continue;
                 drug_list.push_back(_drugs[pos->second]);
            }
       } catch (const std::exception& exc) {}
}

void ConnectDic::find_variants(const std::string &drugID, std::vector<ConnectFormat>& drug_list)
{
       drug_list.clear();

       try { 
            const ConnectFormat& drug = find_drug(drugID);

            std::pair<std::multimap<std::string, int>::iterator, std::multimap<std::string, int>::iterator> range = _VariantConnMap.equal_range(drugID);
            std::multimap<std::string, int>::iterator pos;
            for (pos = range.first; pos != range.second; ++pos) {
                 if (_drugs[pos->second].drugname() == drug.drugname()) continue;
                 drug_list.push_back(_drugs[pos->second]);
            }
       } catch (const std::exception& exc) {}
}

const std::set<std::string>& ConnectDic::get_connected_atoms(const std::string& drugID)
{
       try { 
            const ConnectFormat& drug = find_drug(drugID);
            return drug.getTerminalConnectedAtoms();
       } catch (const std::exception& exc) {}
       return _empty_set;
}

bool ConnectDic::BinaryConvertor(const std::string& textfile, const std::string& binaryfile)
{
       std::string cs;
       _fObj = read_cif_file(binaryfile, textfile, cs);
       if (!cs.empty()) {
            fprintf(stdout, "Read %s failed\n\n%s\n", textfile.c_str(), cs.c_str());
            if (_fObj) delete _fObj;
            return false;
       }

       ISTable *Table = new ISTable("index");
       Table->AddColumn("id");
       Table->AddColumn("formula");
       Table->AddColumn("subcomponent_list");
       Table->AddColumn("three_letter_code");
       ISTable *Table1 = new ISTable("obsolete");
       Table1->AddColumn("id");
       Table1->AddColumn("supercede");

       int row = 0;
       int row1 = 0;
       std::vector<std::string> blockNames;
       _fObj->GetBlockNames(blockNames);
       for (unsigned int i = 0; i < blockNames.size(); ++i) {
            Block &block = _fObj->GetBlock(blockNames[i]);
            ISTable *t = getTablePtr(block, "chem_comp");
            if (!t) continue;

            Table->AddRow();
            get_value_clean(cs, t, 0, "id");
            Table->UpdateCell(row, "id", cs);
            get_value_clean(cs, t, 0, "formula");
            Table->UpdateCell(row, "formula", cs);
            get_value_clean(cs, t, 0, "pdbx_subcomponent_list");
            Table->UpdateCell(row, "subcomponent_list", cs);
            get_value_clean(cs, t, 0, "three_letter_code");
            Table->UpdateCell(row, "three_letter_code", cs);
            row++;

            get_value_clean_upper(cs, t, 0, "pdbx_release_status");
            if (cs == "OBS") {
                 Table1->AddRow();
                 get_value_clean(cs, t, 0, "id");
                 Table1->UpdateCell(row1, "id", cs);
                 get_value_clean(cs, t, 0, "pdbx_replaced_by");
                 Table1->UpdateCell(row1, "supercede", cs);
                 row1++;
            }
       }

       _fObj->AddBlock("INDEX");
       Block &block = _fObj->GetBlock("INDEX");
       block.WriteTable(Table);

       _fObj->AddBlock("OBSOLETE");
       Block &block1 = _fObj->GetBlock("OBSOLETE");
       block1.WriteTable(Table1);

       _fObj->Close();

       return true;
}



bool ConnectDic::hasValidPhosphorylGroup(const std::string& drugID)
{
       try { 
            const ConnectFormat& drug = find_drug(drugID);
            if (drug.hasValidPhosphorylGroup()) return true;
            return false;
       } catch (const std::exception& exc) {
            return false;
       }
}
