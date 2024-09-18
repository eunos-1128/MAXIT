/*
FILE:     Ndb2Pdb.C
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

#include "Maxit.h"
#include "NdbRefn.h"
#include "NdbToken.h"
#include "PdbWrite.h"
#include "TypeDef.h"
#include "utillib.h"

extern const char *_pdb_source_tokens[NUM_PDB_SOURCE][2];

void Maxit::ndb_to_pdb(const bool& biological_flag)
{
       _ndb_to_pdb_proc_header();
       // if (_title_suppression) _pdb_records.erase("TITLE");
       _ndb_to_pdb_update_CAVEAT();
       if (!biological_flag) _ndb_to_pdb_processing_compnd_source();
       if (_molecules.size() > 1) {
            _updateRecordFront("NUMMDL", 2, String::IntToString(_molecules.size()), true);
       }
       _ndb_to_pdb_processing_Ca_and_P_atom_only();
       _ndb_to_pdb_processing_authors();
       _ndb_to_pdb_update_SITE(); // has to be done before generating all remarks

       std::map<std::string, std::string> id_abbrev_mapping;
       _ndb_to_pdb_update_REF(id_abbrev_mapping);
       _ndb_to_pdb_update_REFN(id_abbrev_mapping);
       _ndb_to_pdb_update_PMID();
       _ndb_to_pdb_update_CRYST1(); 
       if (!_pdb_header_only_flag) _ndb_to_pdb_update_ORIGX();
       _updateRecords(NDB_FILE_FORMAT_PDB);
       if (!biological_flag) _ndb_to_pdb_get_remarks(); 

       _ndb_to_pdb_update_EXPDTA();
       if (!biological_flag && !_pdb_header_only_flag) _ndb_to_pdb_update_DBREF();
       _ndb_to_pdb_update_SEQRES();
       if (!biological_flag && !_pdb_header_only_flag) _ndb_to_pdb_update_MODRES();
       _ndb_to_pdb_get_HET_info();
       if (!_pdb_header_only_flag) {
            _ndb_to_pdb_update_HELIX();
            _ndb_to_pdb_update_SHEET();
            _ndb_to_pdb_update_TURN();
            _ndb_to_pdb_update_SSBOND();
            _ndb_to_pdb_update_LINK();
            _ndb_to_pdb_update_SLTBRG();
       }
       _ndb_to_pdb_update_CISPEP();
       _ndb_to_pdb_update_OBSLTE();
}

void Maxit::write_pdb_file(const bool& segid_flag, const bool& checking_linkage_flag)
{
       std::list<std::pair<RCSB::Atom*, RCSB::Atom*> > atom_connects;
       _ndb_to_pdb_get_atom_connects(checking_linkage_flag, atom_connects);

       if (_pdb_header_only_flag) {
            std::set<std::string> removed_tokens;
            removed_tokens.clear();
            removed_tokens.insert("DBREF");
            removed_tokens.insert("DBREF1");
            removed_tokens.insert("DBREF2");
            removed_tokens.insert("SEQADV");
            removed_tokens.insert("MODRES");
            removed_tokens.insert("SCALE");

            for (std::set<std::string>::const_iterator spos = removed_tokens.begin(); spos != removed_tokens.end(); ++spos) {
                 std::map<std::string, std::list<std::vector<std::string> > >::const_iterator ppos = _pdb_records.find(*spos);
                 if (ppos == _pdb_records.end()) continue;
                 _pdb_records.erase(*spos);
            }
       }

       PdbWrite writer;
       writer.setLog(_logIo);
       writer.setCCDic(_ccDic);
       writer.setRecord(&_pdb_records);
       if (!_pdb_header_only_flag) writer.setMolecule(&_molecules);
       writer.setSEGIDFlag(segid_flag);
       if (!_original_entry_ids.empty()) writer.setOriginalEntryIds(_original_entry_ids);
       if (!atom_connects.empty() && !_pdb_header_only_flag) writer.setConnect(atom_connects);
       writer.setFileName(_output_filename);
       writer.setOutputFormat(NDB_FILE_FORMAT_PDB);
       writer.WriteRecord();
}

void Maxit::_ndb_to_pdb_proc_header()
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator ppos = _pdb_records.find("HEADER");
       if (ppos == _pdb_records.end()) return;

       std::string value;
       _getRecordFront("PDBFIL", 1, value);
       if (value.empty()) _getRecordFront("NDBFIL", 1, value);
       if (value.empty()) _getRecordFront("BMRBID", 1, value);
       /* if (!value.empty()) */ _updateRecordFront("HEADER", 3, value);
}

void Maxit::_ndb_to_pdb_update_CAVEAT()
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find("CAVEAT");
       if (ppos == _pdb_records.end()) return;

       std::string pdbid;
       _getRecordFront("PDBFIL", 1, pdbid);
       if (pdbid.empty()) pdbid = ppos->second.front()[2];

       std::string value = "";
       for (std::list<std::vector<std::string> >::const_iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            if ((*lpos)[3].empty()) continue;
            if (!value.empty()) value += " ";
            value += (*lpos)[3];
       }

       _pdb_records.erase("CAVEAT");
       _updateRecordFront("CAVEAT", 2, pdbid);
       _updateRecordFront("CAVEAT", 3, value);
}

void Maxit::_ndb_to_pdb_processing_compnd_source()
{
       _pdb_records.erase("COMPND");
       _pdb_records.erase("SOURCE");

       int SerialNo = 0; 
       for (std::map<int, Entity>::const_iterator
            epos = _entities.begin(); epos != _entities.end(); ++epos) {
            if (epos->second.chain_type() != "ATOMN" &&
                epos->second.chain_type() != "ATOMP") continue;
            SerialNo++;
            std::string type = epos->second.getValue("src_method");
            _ndb_to_pdb_processing_compnd(SerialNo, type, epos->second);
            _ndb_to_pdb_processing_source(SerialNo, type, epos->second);
       }

       if (!SerialNo) {
            for (std::map<int, Entity>::const_iterator
                 epos = _entities.begin(); epos != _entities.end(); ++epos) {
                 SerialNo++;
                 std::string type = epos->second.getValue("src_method");
                 _ndb_to_pdb_processing_compnd(SerialNo, type, epos->second);
                 _ndb_to_pdb_processing_source(SerialNo, type, epos->second);
                 break;
            }
       }

       // remove last ';'
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("COMPND");
       if (ppos != _pdb_records.end()) {
            std::vector<std::string>& Field = ppos->second.back();
            if (Field[2][Field[2].size() - 1] == ';') Field[2].erase(Field[2].size() - 1);
       }
       ppos = _pdb_records.find("SOURCE");
       if (ppos != _pdb_records.end()) {
            std::vector<std::string>& Field = ppos->second.back();
            if (Field[2][Field[2].size() - 1] == ';') Field[2].erase(Field[2].size() - 1);
       }
}

void Maxit::_ndb_to_pdb_processing_compnd(const int& Mol_ID, const std::string& type,
                                          const Entity& entity)
{
       std::vector<std::pair<std::string, std::string> > key_item_pair;
       key_item_pair.clear();
       key_item_pair.push_back(std::make_pair("MOLECULE: ", "pdbx_description"));
       key_item_pair.push_back(std::make_pair("CHAIN: ", ""));
       key_item_pair.push_back(std::make_pair("FRAGMENT: ", "pdbx_fragment"));
       key_item_pair.push_back(std::make_pair("SYNONYM: ", "com_name"));
       key_item_pair.push_back(std::make_pair("EC: ", "pdbx_ec"));
       key_item_pair.push_back(std::make_pair("ENGINEERED: ", ""));
       key_item_pair.push_back(std::make_pair("MUTATION: ","pdbx_mutation"));
       // key_item_pair.push_back(std::make_pair("BIOLOGICAL_UNIT: ", ""));
       key_item_pair.push_back(std::make_pair("OTHER_DETAILS: ", "details"));

       std::vector<std::string> compnd;
       compnd.clear();
       compnd.push_back("MOL_ID: " + String::IntToString(Mol_ID));

       for (std::vector<std::pair<std::string, std::string> >::const_iterator
            pos = key_item_pair.begin(); pos != key_item_pair.end(); ++pos) {
            if (pos->first == "CHAIN: ") {
                 std::string chainID = "";
                 const std::vector<std::string>& PDB_chainIDs = entity.PDB_chainID();
                 for (std::vector<std::string>::const_iterator
                      vpos = PDB_chainIDs.begin(); vpos != PDB_chainIDs.end(); ++vpos) {
                      if (!chainID.empty()) chainID += ", ";
                      chainID += *vpos;
                 }
                 if (chainID.empty()) chainID = "NULL";
                 compnd.push_back(pos->first + chainID);
            } else if (pos->first == "ENGINEERED: ") {
                 if (type == "man" || type == "syn") compnd.push_back(pos->first + "YES");
            } else {
                 std::string value = entity.getValue(pos->second);
                 if (value.empty() && pos->first != "MOLECULE: ") continue;
                 if (pos->first == "MUTATION: ") value = "YES";
                 String::UpperCase(value);
                 String::StripAndCompressWs(value);
                 compnd.push_back(pos->first + value);
            }
       }

       _ndb_to_pdb_add_compnd_or_source("COMPND", compnd);
}

void Maxit::_ndb_to_pdb_processing_source(const int& Mol_ID, const std::string& type, const Entity& entity)
{
       std::vector<std::string> source, data;
       source.clear();
       source.push_back("MOL_ID: " + String::IntToString(Mol_ID));
       if (type == "syn") source.push_back("SYNTHETIC: YES");

       const std::vector<std::pair<std::string, std::map<std::string, std::string> > >& s_info = entity.source();
       std::vector<std::vector<std::string> > values;
       values.clear();
       for (std::vector<std::pair<std::string, std::map<std::string, std::string> > >::const_iterator
            pos = s_info.begin(); pos != s_info.end(); ++pos) {
            if (pos->first != type) continue;
            data.clear();
            data.push_back("");
            for (int i = 1; i < NUM_PDB_SOURCE; ++i) {
                 std::map<std::string, std::string>::const_iterator
                     mpos = pos->second.find(_pdb_source_tokens[i][1]);
                 if (mpos == pos->second.end() || mpos->second == "?")
                      data.push_back("");
                 else data.push_back(mpos->second);
            }
            values.push_back(data);
       }

       std::set<std::string> tmp_set;
       std::string val;
       for (int i = 1; i < NUM_PDB_SOURCE; ++i) {
            tmp_set.clear();
            val.clear();
            for (unsigned int j = 0; j < values.size(); ++j) {
                 if (values[j][i].empty()) continue;
                 if (tmp_set.find(values[j][i]) != tmp_set.end()) continue;
                 tmp_set.insert(values[j][i]);
                 if (!val.empty()) val += ", ";
                 val += values[j][i];
            }
            String::StripAndCompressWs(val);
            if (val.empty()) continue;

            std::string buffer = _pdb_source_tokens[i][0];
            source.push_back(buffer + " " + val);
       }

       _ndb_to_pdb_add_compnd_or_source("SOURCE", source);
}

void Maxit::_ndb_to_pdb_add_compnd_or_source(const std::string& token, const std::vector<std::string>& data)
{
       for (unsigned int i = 0; i < data.size(); ++i) {
            std::string cs = data[i];
            String::StripTrailingWs(cs);
            if (cs[cs.size() - 1] != ';') cs += ";"; 
            _addNewRecord(token);
            _updateRecordBack(token, 2, cs);
            _updateRecordBack(token, 3, String::IntToString(i + 1));
       }
}

void Maxit::_ndb_to_pdb_update_EXPDTA()
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("EXPDTA");
       if (ppos == _pdb_records.end()) return;

       if (ppos->second.size() < 2) return;

       std::string cs = "";
       for (std::list<std::vector<std::string> >::iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            if (!cs.empty()) cs += "; ";
            cs += (*lpos)[2];
       }

       std::list<std::vector<std::string> >::iterator lpos = ppos->second.begin();
       (*lpos)[2] = cs;

       lpos++;
       while (lpos != ppos->second.end()) lpos = ppos->second.erase(lpos);
}

void Maxit::_ndb_to_pdb_processing_Ca_and_P_atom_only()
{
       if (_molecules.empty()) return;

       std::string Ca_list, P_list;
       Ca_list.clear(); P_list.clear();
       RCSB::Chain* chain = _molecules[0]->GetFirstChain();
       while (chain) {
            if (chain->chain_type() == "ATOMP" || chain->chain_type() == "ATOMN") {
                 if (chain->ca_or_p_atom_only() == CAAtom_ONLY) {
                      if (!Ca_list.empty()) Ca_list += ", ";
                      Ca_list += chain->PDB_ChainID();
                 } else if (chain->ca_or_p_atom_only() == PAtom_ONLY) {
                      if (!P_list.empty()) P_list += ", ";
                      P_list += chain->PDB_ChainID();
                 }
            }
            chain = _molecules[0]->GetNextChain();
       }
       if (Ca_list.empty() && P_list.empty()) return;

       std::string value;
       _getRecordFront("MDLTYP", 2, value);
       if (!Ca_list.empty() && value.find("CA ATOMS ONLY, CHAIN") == std::string::npos) {
            if (!value.empty()) value += "; ";
            value += "CA ATOMS ONLY, CHAIN " + Ca_list;
       }
       if (!P_list.empty() && value.find("P ATOMS ONLY, CHAIN") == std::string::npos) {
            if (!value.empty()) value += "; ";
            value += "P ATOMS ONLY, CHAIN " + P_list;
       }

       _updateRecordFront("MDLTYP", 2, value);
}

void Maxit::_ndb_to_pdb_processing_authors()
{
       _removeWhiteSpaceBetweenName("AUTHOR", 2); 
       _removeWhiteSpaceBetweenName("AUTH", 4); 

       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
            ppos = _pdb_records.find("AUTHOR");
       if (ppos == _pdb_records.end()) {
            std::string value;
            _getRecordFront("AUTH", 4, value);
            if (!value.empty()) {
                 _addNewRecord("AUTHOR");
                 _updateRecordBack("AUTHOR", 2, value);
            }   
       }
}

void Maxit::_ndb_to_pdb_update_REF(std::map<std::string, std::string>& id_abbrev_mapping)
{
       id_abbrev_mapping.clear();
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("REF");
       if (ppos == _pdb_records.end()) return;

       int vField = 5;
       for (std::list<std::vector<std::string> >::iterator
            lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            id_abbrev_mapping.insert(std::make_pair((*lppos)[1], (*lppos)[4]));
            if (!(*lppos)[vField + 1].empty() && (*lppos)[vField].empty())
                 (*lppos)[vField] = "V.";
            else if ((*lppos)[vField + 1].empty())
                 (*lppos)[vField].clear();
       }
}

void Maxit::_ndb_to_pdb_update_REFN(const std::map<std::string, std::string>& id_abbrev_mapping)
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("REFN");
       if (ppos == _pdb_records.end()) return;

       NdbRefn::initialize();
       NdbRefn::Read(*_logIo, _rcsbroot);

       int astmField = 3;
       int countryField = 5;
       int issnField = 6;
       for (std::list<std::vector<std::string> >::iterator
            lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            (*lppos)[astmField].clear();
            (*lppos)[astmField + 1].clear();
            (*lppos)[countryField].clear();

            std::map<std::string, std::string>::const_iterator
                mpos = id_abbrev_mapping.find((*lppos)[1]);
            if (mpos != id_abbrev_mapping.end()) {
                 std::map<std::string, std::string> jrnlInfo = NdbRefn::GetJournalInfo(mpos->second,
                                       (*lppos)[issnField + 1]);
                 if (!jrnlInfo.empty()) (*lppos)[issnField] = jrnlInfo["issn"];
            }

            if (!(*lppos)[issnField + 1].empty() && (*lppos)[issnField].empty())
                 (*lppos)[issnField] = "ISSN";
            else if ((*lppos)[issnField + 1].empty())
                 (*lppos)[issnField].clear();

            if ((*lppos)[issnField] == "ISBN") {
                 (*lppos)[issnField].clear();
                 (*lppos)[issnField + 1].clear();
            }
       }
}

void Maxit::_ndb_to_pdb_update_PMID()
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("PMID");
       if (ppos == _pdb_records.end()) return;

       for (std::list<std::vector<std::string> >::iterator
            lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            if ((*lppos)[4] == "-1") {
                 lppos = ppos->second.erase(lppos);
                 --lppos;
            }
       }
       if (ppos->second.empty()) _pdb_records.erase("PMID");
}

void Maxit::_ndb_to_pdb_update_DBREF()
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
            ppos = _pdb_records.find("DBREF");
       if (ppos == _pdb_records.end()) return;

       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("DBREF");

       for (std::list<std::vector<std::string> >::iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            if ((*lpos)[1].size() > 4) (*lpos)[1].clear();
            if ((int) (*lpos)[8].size() > ndbformat.FieldList[8].FieldWidth ||
                (int) (*lpos)[9].size() > ndbformat.FieldList[9].FieldWidth ||
                (int) (*lpos)[10].size() > ndbformat.FieldList[10].FieldWidth ||
                (int) (*lpos)[12].size() > ndbformat.FieldList[12].FieldWidth) {
                 std::vector<std::string> data = *lpos;
                 data[0] = "DBREF1";
                 data[8] = data[9];
                 (*lpos)[0] = "DBREF2";
                 (*lpos)[3] = (*lpos)[8];
                 (*lpos)[4] = (*lpos)[10];
                 (*lpos)[5] = (*lpos)[11];
                 (*lpos)[6] = (*lpos)[12];
                 (*lpos)[7] = (*lpos)[13];
                 lpos = ppos->second.insert(lpos, data);
                 lpos++;
            }
       }
}

void Maxit::_ndb_to_pdb_update_SEQRES()
{
       if (_molecules.empty()) return;

       _pdb_records.erase("SEQRES");

       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("SEQRES");

       std::list<std::string> seqs;
       RCSB::Chain* chain = _molecules[0]->GetFirstChain();
       while (chain) {
            if ((chain->chain_type() != "ATOMN") && (chain->chain_type() != "ATOMP")) {
                 chain = _molecules[0]->GetNextChain();
                 continue;
            }
            chain->get_seq(seqs);
            if (seqs.empty()) {
                 chain = _molecules[0]->GetNextChain();
                 continue;
            }

            int serialNo = 1;
            _addNewRecord("SEQRES");
            _updateRecordBack("SEQRES", 1, String::IntToString(serialNo));
            _updateRecordBack("SEQRES", 2, chain->PDB_ChainID());
            _updateRecordBack("SEQRES", 3, String::IntToString(seqs.size()));
            int seqField = 4;
            for (std::list<std::string>::const_iterator
                 lpos = seqs.begin(); lpos != seqs.end(); ++lpos) {
                 if (seqField == ndbformat.NumField) {
                      serialNo++;
                      _addNewRecord("SEQRES");
                      _updateRecordBack("SEQRES", 1, String::IntToString(serialNo));
                      _updateRecordBack("SEQRES", 2, chain->PDB_ChainID());
                      _updateRecordBack("SEQRES", 3, String::IntToString(seqs.size()));
                      seqField = 4;
                 }
                 _updateRecordBack("SEQRES", seqField, *lpos);
                 seqField++;
            }
            chain = _molecules[0]->GetNextChain();
       }
}

void Maxit::_ndb_to_pdb_update_CISPEP()
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
            ppos = _pdb_records.find("CISPEP");
       if (ppos == _pdb_records.end()) return;

       if (_molecules.size() > 1) return;

       for (std::list<std::vector<std::string> >::iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            (*lpos)[10] = "0";
       }
}

void Maxit::_ndb_to_pdb_update_CRYST1()
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
            ppos = _pdb_records.find("CRYST1");
       if (ppos != _pdb_records.end()) return;

       _addNewRecord("CRYST1");
       _updateRecordBack("CRYST1", 1, "1.000");
       _updateRecordBack("CRYST1", 2, "1.000");
       _updateRecordBack("CRYST1", 3, "1.000");
       _updateRecordBack("CRYST1", 4, "90.00");
       _updateRecordBack("CRYST1", 5, "90.00");
       _updateRecordBack("CRYST1", 6, "90.00");
       _updateRecordBack("CRYST1", 7, "P 1");
       _updateRecordBack("CRYST1", 8, "1");
}

void Maxit::_ndb_to_pdb_update_ORIGX()
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
            ppos = _pdb_records.find("ORIGX");
       if (ppos != _pdb_records.end()) return;

       _addNewRecord("ORIGX");
       _updateRecordBack("ORIGX", 1, "1.000000");
       _updateRecordBack("ORIGX", 2, "0.000000");
       _updateRecordBack("ORIGX", 3, "0.000000");
       _updateRecordBack("ORIGX", 4, "0.00000");

       _addNewRecord("ORIGX");
       _updateRecordBack("ORIGX", 1, "0.000000");
       _updateRecordBack("ORIGX", 2, "1.000000");
       _updateRecordBack("ORIGX", 3, "0.000000");
       _updateRecordBack("ORIGX", 4, "0.00000");

       _addNewRecord("ORIGX");
       _updateRecordBack("ORIGX", 1, "0.000000");
       _updateRecordBack("ORIGX", 2, "0.000000");
       _updateRecordBack("ORIGX", 3, "1.000000");
       _updateRecordBack("ORIGX", 4, "0.00000");
}

void Maxit::_ndb_to_pdb_update_OBSLTE()
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
            ppos = _pdb_records.find("OBSLTE");
       if (ppos == _pdb_records.end()) return;

       for (std::list<std::vector<std::string> >::iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            if (String::IsEqual((*lpos)[4], "NONE", Char::eCASE_INSENSITIVE)) (*lpos)[4].clear();
       }
}
