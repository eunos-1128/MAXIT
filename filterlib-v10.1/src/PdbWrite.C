/*
FILE:     PdbWrite.C
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
#include <math.h>
#include <sys/stat.h>

#include "Element.h"
#include "PdbWrite.h"
#include "SeqCodeUtil.h"
#include "TypeDef.h"
#include "utillib.h"

#define NUM_AUX_ATOM   3

static const char *_aux_atoms[NUM_AUX_ATOM] = { "SIGATM", "ANISOU", "SIGUIJ" };

#define NUM_RESIDUE_NAME_LABEL  12

static const char *_residue_name_labels[NUM_RESIDUE_NAME_LABEL] = {
       "Beg_NDB_res_name",
       "End_NDB_res_name",
       "H_bond_cur_str_NDB_res_name",
       "H_bond_prev_str_NDB_res_name",
       "Het_Group_ID",
       "NDB_Res_Name_CYS_1",
       "NDB_Res_Name_CYS_2",
       "Residue_Name",
       "Residue_Name_1",
       "Residue_Name_2",
       "Residue_Name_3",
       "Residue_Name_4"
};

#define NUM_CHAIN_ID_LABEL  12

static const char *_chain_id_labels[NUM_CHAIN_ID_LABEL] = {
       "Beg_NDB_Strand_ID",
       "End_NDB_Strand_ID",
       "H_bond_cur_str_NDB_Strand_ID",
       "H_bond_prev_str_NDB_Strand_ID",
       "Het_Strand_ID",
       "NDB_Strand_ID_CYS_1",
       "NDB_Strand_ID_CYS_2",
       "Strand_ID",
       "Strand_ID_1",
       "Strand_ID_2",
       "Strand_ID_3",
       "Strand_ID_4"
};

PdbWrite::PdbWrite()
{
       _init();
}

PdbWrite::~PdbWrite()
{
       _logIo = NULL;
       _ccDic = NULL;
       _pdb_records = NULL;
       _molecules = NULL;
       _segid_flag = false;
       _original_entry_ids.clear();
       _atom_connects.clear();
       _filename.clear();
       _atomField_mapping.clear();
       _field_length_mapping.clear();
       _affected_field_set.clear();
}

void PdbWrite::setLog(LogUtil *logPt)
{
       _logIo = logPt;
}

void PdbWrite::setCCDic(ConnectDic *ccdic)
{
       _ccDic = ccdic;
}

void PdbWrite::setRecord(std::map<std::string, std::list<std::vector<std::string> > >* records)
{
       _pdb_records = records;
}

void PdbWrite::setMolecule(std::vector<RCSB::Molecule*>* mols)
{
       _molecules = mols;
}

void PdbWrite::setSEGIDFlag(const bool& segid_flag)
{
       _segid_flag = segid_flag;
}

void PdbWrite::setOriginalEntryIds(const std::vector<std::string>& ids)
{
       _original_entry_ids = ids;
}

void PdbWrite::setConnect(const std::list<std::pair<RCSB::Atom*, RCSB::Atom*> >& connects)
{
       _atom_connects = connects;
}

void PdbWrite::setFileName(const std::string& filename)
{
       _filename = filename;
}

void PdbWrite::setOutputFormat(const int& format)
{
       _output_format = format;
}

void PdbWrite::setExtendedFlag()
{
       _extended_flag = true;
}

void PdbWrite::setFieldLabel(const std::string& fieldLabel)
{
       _affected_field_set.insert(fieldLabel);
}

void PdbWrite::WriteRecord()
{
       _write();
}

void PdbWrite::WriteRecord(FILE *fp)
{
       _write(fp);
}

void PdbWrite::WriteRecord(std::string& text)
{
       if (_output_format == NDB_FILE_FORMAT_PDB) NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_PDB);

       text.clear();
       std::vector<std::string> lines;
       const std::vector<std::string>& ndbTokens = NdbToken::getNdbTokens();
       for (std::vector<std::string>::const_iterator tpos = ndbTokens.begin(); tpos != ndbTokens.end(); ++tpos) {
            if (*tpos == "ATOMN") {
                 WriteCoordinates(text);
                 continue;
            } else if (*tpos == "JRNL") {
                 _writeJrnlRecords(text);
                 continue;
            }

            if (!_pdb_records) continue;

            std::map<std::string, std::list<std::vector<std::string> > >::iterator pos = _pdb_records->find(*tpos);
            if (pos == _pdb_records->end()) continue;

            const ndb_token_format& ndbformat = NdbToken::getTokenFormat(*tpos);
            writeGeneralRecords(lines, ndbformat, pos->second);
            for (std::vector<std::string>::const_iterator vpos = lines.begin(); vpos != lines.end(); ++vpos) {
                 text += *vpos + "\n";
            }
       }

       if (_output_format == NDB_FILE_FORMAT_PDB) {
            // print_master_card(fp, numRemark);
       }

       std::vector<std::string> FieldInfo;
       FieldInfo.clear();
       FieldInfo.push_back("END");
       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("END");
       writeGeneralRecord(lines, ndbformat, FieldInfo);
       for (std::vector<std::string>::const_iterator vpos = lines.begin(); vpos != lines.end(); ++vpos) {
            text += *vpos + "\n";
       }

       if (_output_format == NDB_FILE_FORMAT_PDB) NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_NDB);
}

void PdbWrite::WriteRecord(const std::string& filename, const int format, std::map<std::string, std::list<std::vector<std::string> > >* records,
                           std::vector<RCSB::Molecule*>* mols, const std::list<std::pair<RCSB::Atom*, RCSB::Atom*> >& connects)
{
       _filename = filename;
       _output_format = format;
       _pdb_records = records;
       _molecules = mols;
       _atom_connects = connects;
       _write();
}

void PdbWrite::WriteRecord(const std::string& filename, const std::vector<RCSB::Residue*>& residues)
{
       if (_output_format == NDB_FILE_FORMAT_PDB) NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_PDB);

       FILE *fp = fopen(filename.c_str(), "w");

       std::vector<std::string> FieldInfo;
       for (std::vector<RCSB::Residue*>::const_iterator rpos = residues.begin(); rpos != residues.end(); ++rpos) {
            RCSB::Atom* atom = (*rpos)->GetFirstAtom();
            while (atom) {
                 WriteAtom(fp, "HETAIN", atom, FieldInfo, false);
                 atom = (*rpos)->GetNextAtom();          
            }
       }

       FieldInfo.clear();
       FieldInfo.push_back("END");
       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("END");
       writeGeneralRecord(fp, ndbformat, FieldInfo);

       fclose (fp);

       if (_output_format == NDB_FILE_FORMAT_PDB) NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_NDB);
}

void PdbWrite::WriteRecord(const std::string& filename, const std::vector<RCSB::Atom*>& atoms, const std::list<std::pair<RCSB::Atom*, RCSB::Atom*> >& connects)
{
       if (_output_format == NDB_FILE_FORMAT_PDB) NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_PDB);

       FILE *fp = fopen(filename.c_str(), "w");

       std::vector<std::string> FieldInfo;
       for (std::vector<RCSB::Atom*>::const_iterator apos = atoms.begin(); apos != atoms.end(); ++apos) {
            WriteAtom(fp, "HETAIN", *apos, FieldInfo, false);
       }

       if (!connects.empty() && (_output_format == NDB_FILE_FORMAT_PDB)) {
            _atom_connects = connects;

            std::map<std::string, std::list<std::vector<std::string> > > records;
            records.clear();
            _pdb_records = &records;
            _getConnectRecord();

            std::map<std::string, std::list<std::vector<std::string> > >::iterator mpos = records.find("CONECT");
            if (mpos != records.end()) {
                 const ndb_token_format& ndbformat = NdbToken::getTokenFormat("CONECT");
                 writeGeneralRecords(fp, ndbformat, mpos->second);
            }
       }

       FieldInfo.clear();
       FieldInfo.push_back("END");
       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("END");
       writeGeneralRecord(fp, ndbformat, FieldInfo);

       fclose (fp);

       if (_output_format == NDB_FILE_FORMAT_PDB) NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_NDB);
}

void PdbWrite::WriteRecord(std::string& text, const std::list<std::pair<std::string, RCSB::Atom*> >& atom_list)
{
       if (_output_format == NDB_FILE_FORMAT_PDB) NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_PDB);

       text.clear();
       std::vector<std::string> FieldInfo, lines;
       for (std::list<std::pair<std::string, RCSB::Atom*> >::const_iterator apos = atom_list.begin(); apos != atom_list.end(); ++apos) {
            WriteAtom(lines, apos->first, apos->second, FieldInfo, false);
            for (std::vector<std::string>::const_iterator vpos = lines.begin(); vpos != lines.end(); ++vpos) {
                 text += *vpos + "\n";
            }
       }
/*
       FieldInfo.clear();
       FieldInfo.push_back("END");
       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("END");
       writeGeneralRecord(lines, ndbformat, FieldInfo);
       for (std::vector<std::string>::const_iterator vpos = lines.begin(); vpos != lines.end(); ++vpos) {
            text += *vpos + "\n";
       }
*/
       if (_output_format == NDB_FILE_FORMAT_PDB) NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_NDB);
}

void PdbWrite::WriteRecord(const std::string& filename, const std::list<std::pair<std::string, RCSB::Atom*> >& atom_list,
                           const std::list<std::pair<RCSB::Atom*, RCSB::Atom*> >& connects)
{
       if (_output_format == NDB_FILE_FORMAT_PDB) NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_PDB);

       FILE *fp = fopen(filename.c_str(), "w");

       std::vector<std::string> FieldInfo;
       for (std::list<std::pair<std::string, RCSB::Atom*> >::const_iterator apos = atom_list.begin(); apos != atom_list.end(); ++apos) {
            WriteAtom(fp, apos->first, apos->second, FieldInfo, false);
       }

       if (!connects.empty() && (_output_format == NDB_FILE_FORMAT_PDB)) {
            _atom_connects = connects;

            std::map<std::string, std::list<std::vector<std::string> > > records;
            records.clear();
            _pdb_records = &records;
            _getConnectRecord();

            std::map<std::string, std::list<std::vector<std::string> > >::iterator mpos = records.find("CONECT");
            if (mpos != records.end()) {
                 const ndb_token_format& ndbformat = NdbToken::getTokenFormat("CONECT");
                 writeGeneralRecords(fp, ndbformat, mpos->second);
            }
       }

       fclose (fp);

       if (_output_format == NDB_FILE_FORMAT_PDB) NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_NDB);
}

void PdbWrite::writeGeneralRecords(FILE *fp, std::map<std::string, std::list<std::vector<std::string> > >* records)
{
       _pdb_records = records;
       _write(fp, true);
}

void PdbWrite::writeGeneralRecords(std::vector<std::string>& lines, const ndb_token_format& ndbformat, std::list<std::vector<std::string> >& records,
                                   const bool& full_record_flag)
{
       lines.clear();

       if (_output_format == NDB_FILE_FORMAT_PDB && ndbformat.NdbOnly) return;
       if (records.empty()) return;

       std::string prevSeqNo = "-1";
       int PdbContNo = 1;
       for (std::list<std::vector<std::string> >::iterator pos = records.begin(); pos != records.end(); ++pos) {
            if ((*pos)[0] == "HETNAM" || (*pos)[0] == "HETSYN") PdbContNo = 1;
            if ((*pos)[0] == "DBREF1") {
                 const ndb_token_format& ndbformat1 = NdbToken::getTokenFormat("DBREF1");
                 _writeGeneralRecord(lines, ndbformat1, *pos, prevSeqNo, PdbContNo, full_record_flag);
            } else if ((*pos)[0] == "DBREF2") {
                 const ndb_token_format& ndbformat1 = NdbToken::getTokenFormat("DBREF2");
                 _writeGeneralRecord(lines, ndbformat1, *pos, prevSeqNo, PdbContNo, full_record_flag);
            } else _writeGeneralRecord(lines, ndbformat, *pos, prevSeqNo, PdbContNo, full_record_flag);
       }
}

void PdbWrite::writeGeneralRecords(FILE *fp, const ndb_token_format& ndbformat, std::list<std::vector<std::string> >& records,
                                   const bool& full_record_flag)
{
       if (_output_format == NDB_FILE_FORMAT_PDB && ndbformat.NdbOnly) return;
       if (records.empty()) return;

       std::vector<std::string> lines;
       writeGeneralRecords(lines, ndbformat, records, full_record_flag);
       for (std::vector<std::string>::const_iterator pos = lines.begin(); pos != lines.end(); ++pos) {
            fprintf(fp, "%s\n", pos->c_str());
       }
}

void PdbWrite::writeGeneralRecord(std::vector<std::string>& lines, const ndb_token_format& ndbformat, std::vector<std::string>& record,
                                   const bool& full_record_flag)
{
       lines.clear();

       if (_output_format == NDB_FILE_FORMAT_PDB && ndbformat.NdbOnly) return;
       if (record.empty()) return;

       std::string prevSeqNo = "-1";
       int PdbContNo = 1;
       if (record[0] == "DBREF1") {
            const ndb_token_format& ndbformat1 = NdbToken::getTokenFormat("DBREF1");
            _writeGeneralRecord(lines, ndbformat1, record, prevSeqNo, PdbContNo, full_record_flag);
       } else if (record[0] == "DBREF2") {
            const ndb_token_format& ndbformat1 = NdbToken::getTokenFormat("DBREF2");
            _writeGeneralRecord(lines, ndbformat1, record, prevSeqNo, PdbContNo, full_record_flag);
       } else _writeGeneralRecord(lines, ndbformat, record, prevSeqNo, PdbContNo, full_record_flag);
}

void PdbWrite::writeGeneralRecord(FILE *fp, const ndb_token_format& ndbformat, std::vector<std::string>& record, const bool& full_record_flag)
{
       if (_output_format == NDB_FILE_FORMAT_PDB && ndbformat.NdbOnly) return;
       if (record.empty()) return;

       std::vector<std::string> lines;
       writeGeneralRecord(lines, ndbformat, record, full_record_flag);
       for (std::vector<std::string>::const_iterator pos = lines.begin(); pos != lines.end(); ++pos) {
            fprintf(fp, "%s\n", pos->c_str());
       }
}

void PdbWrite::WriteCoordinates(std::string& text)
// writing to file only, all other data manipulation were stripped out.
{
       if (!_molecules) return;

       if (_output_format == NDB_FILE_FORMAT_PDB) NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_PDB);

       std::vector<std::string> FieldInfo, lines;
       for (unsigned int l = 0; l < _molecules->size(); ++l) { 
            if (_molecules->size() > 1) {
                 FieldInfo.clear();
                 FieldInfo.push_back("MODEL");
                 FieldInfo.push_back(String::IntToString((*_molecules)[l]->Mol_ID()));
                 const ndb_token_format& ndbformat = NdbToken::getTokenFormat("MODEL");
                 writeGeneralRecord(lines, ndbformat, FieldInfo);
                 for (std::vector<std::string>::const_iterator vpos = lines.begin(); vpos != lines.end(); ++vpos) {
                      text += *vpos + "\n";
                 }
            }

            WriteMolecule(text, (*_molecules)[l]);

            if (_molecules->size() > 1) {
                 FieldInfo.clear();
                 FieldInfo.push_back("ENDMDL");
                 const ndb_token_format& ndbformat = NdbToken::getTokenFormat("ENDMDL");
                 writeGeneralRecord(lines, ndbformat, FieldInfo);
                 for (std::vector<std::string>::const_iterator vpos = lines.begin(); vpos != lines.end(); ++vpos) {
                      text += *vpos + "\n";
                 }
            }
       } 

       if (_output_format == NDB_FILE_FORMAT_PDB) NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_NDB);
}

void PdbWrite::WriteCoordinates(FILE *fp)
// writing to file only, all other data manipulation were stripped out.
{
       if (!_molecules) return;

       if (_output_format == NDB_FILE_FORMAT_PDB) NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_PDB);

       std::vector<std::string> FieldInfo;
       for (unsigned int l = 0; l < _molecules->size(); ++l) { 
            if (_molecules->size() > 1) {
                 FieldInfo.clear();
                 FieldInfo.push_back("MODEL");
                 FieldInfo.push_back(String::IntToString((*_molecules)[l]->Mol_ID()));
                 const ndb_token_format& ndbformat = NdbToken::getTokenFormat("MODEL");
                 writeGeneralRecord(fp, ndbformat, FieldInfo);
            }

            _num_ter = 0;
            _num_atom = 0;
            WriteMolecule(fp, (*_molecules)[l]);

            if (l == 0 && _output_format == NDB_FILE_FORMAT_PDB) {
                 _getConnectRecord();
                 _getMasterRecord();
            }

            if (_molecules->size() > 1) {
                 FieldInfo.clear();
                 FieldInfo.push_back("ENDMDL");
                 const ndb_token_format& ndbformat = NdbToken::getTokenFormat("ENDMDL");
                 writeGeneralRecord(fp, ndbformat, FieldInfo);
            }
       } 

       if (_output_format == NDB_FILE_FORMAT_PDB) NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_NDB);
}

void PdbWrite::WriteMolecule(std::string& text, RCSB::Molecule* molecule)
{
       int atomSerialNo = 0;
       std::vector<std::string> FieldInfo, lines;
       std::vector<RCSB::Residue*> residue_list;
       RCSB::Chain* chain = molecule->GetFirstChain();
       while (chain) {
            chain->GetFirstResidueList(residue_list);
            while (!residue_list.empty()) {
                 for (std::vector<RCSB::Residue*>::const_iterator rpos = residue_list.begin(); rpos != residue_list.end(); ++rpos) {
                      RCSB::Atom* atom = (*rpos)->GetFirstAtom();
                      while (atom) {
                           atomSerialNo++;
                           atom->set_atnum(atomSerialNo);
                           WriteAtom(lines, chain->chain_type(), atom, FieldInfo);
                           for (std::vector<std::string>::const_iterator vpos = lines.begin(); vpos != lines.end(); ++vpos) {
                                text += *vpos + "\n";
                           }
                           atom = (*rpos)->GetNextAtom();
                      }
                 }
                 chain->GetNextResidueList(residue_list);
            }

            if (!chain->empty_chain() && (chain->chain_type() == "ATOMN" || chain->chain_type() == "ATOMP")) {
                 atomSerialNo++;
                 if (atomSerialNo > 99999)
                      FieldInfo[1] = "99999";
                 else FieldInfo[1] = String::IntToString(atomSerialNo);
                 _write_ter_card(text, FieldInfo);
            }

            chain = molecule->GetNextChain();
       }
}

void PdbWrite::WriteMolecule(FILE *fp, RCSB::Molecule* molecule)
{
       std::set<std::string> index_set;
       int atomSerialNo = 0;
       std::vector<std::string> FieldInfo;
       std::vector<RCSB::Residue*> residue_list;
       RCSB::Chain* chain = molecule->GetFirstChain();
       while (chain) {
            chain->GetFirstResidueList(residue_list);
            while (!residue_list.empty()) {
                 for (std::vector<RCSB::Residue*>::const_iterator rpos = residue_list.begin(); rpos != residue_list.end(); ++rpos) {
                      index_set.clear();
                      RCSB::Atom* atom = (*rpos)->GetFirstAtom();
                      while (atom) {
                           std::string index = atom->getAtomIndex();
                           if (index_set.find(index) == index_set.end()) {
                                if (!atom->is_hydrogen()) _num_atom++;
                                index_set.insert(index);
                           }
                           atomSerialNo++;
                           atom->set_atnum(atomSerialNo);
                           WriteAtom(fp, chain->chain_type(), atom, FieldInfo);
                           atom = (*rpos)->GetNextAtom();
                      }
                 }
                 chain->GetNextResidueList(residue_list);
            }

            if (!chain->empty_chain() && (chain->chain_type() == "ATOMN" || chain->chain_type() == "ATOMP")) {
                 _num_ter++;
                 atomSerialNo++;
                 if (atomSerialNo > 99999)
                      FieldInfo[1] = "99999";
                 else FieldInfo[1] = String::IntToString(atomSerialNo);
                 _write_ter_card(fp, FieldInfo);
            }

            chain = molecule->GetNextChain();
       }
}

void PdbWrite::WriteAtom(std::vector<std::string>& records, const std::string& chain_type, RCSB::Atom* atom, std::vector<std::string>& FieldInfo,
                         const bool include_auxiliary_atoms)
{
       records.clear();

       FieldInfo = atom->getValue();
       if (atoi(FieldInfo[1].c_str()) > 99999) FieldInfo[1] = "99999";
       if (FieldInfo[6].empty()) FieldInfo[6] = FieldInfo[15];

       // set proper atom record token
       if (!chain_type.empty()) {
            if (_output_format != NDB_FILE_FORMAT_NDB) {
                 if ((chain_type == "ATOMN" || chain_type == "ATOMP") && SeqCodeUtil::is_a_standard_residue(FieldInfo[4]))
                      FieldInfo[0] = "ATOM"; 
                 else FieldInfo[0] = "HETATM";
            } else FieldInfo[0] = chain_type;
       }

       // set proper Atom Name Field
       FieldInfo[2] = printAtomNameField(_ccDic, FieldInfo[13], FieldInfo[2], FieldInfo[4]);
       FieldInfo[18] = FieldInfo[2];

       // set PDB nomenclature
       if (_output_format == NDB_FILE_FORMAT_PDB) {
            FieldInfo[5] = FieldInfo[17];
            FieldInfo[6] = FieldInfo[15];
            reformatting_charge_to_suffix(FieldInfo[14]);
       }

       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("ATOMN");
       if (_affected_field_set.empty()) {
            writeGeneralRecord(records, ndbformat, FieldInfo);
       } else {
            ndb_token_format updated_ndbformat = ndbformat;
            updateNdbTokenFormat(_affected_field_set, updated_ndbformat);
            writeGeneralRecord(records, updated_ndbformat, FieldInfo);
       }

       if (!include_auxiliary_atoms) {
            if (_segid_flag) {
                 std::string segid = FormattedFieldValue(FieldInfo[5], 3, 4, 0, true, true);
                 for (std::vector<std::string>::iterator
                      vpos = records.begin(); vpos != records.end(); ++vpos) {
                      vpos->replace(20, 2, "  ");
                      vpos->replace(72, 4, segid);
                 }
            }
            return;
       }

       std::vector<std::string> Field, lines;
       
       for (int i = 0; i < NUM_AUX_ATOM; ++i) {
            atom->getAuxiliaryValue(_aux_atoms[i], Field);
            if (Field.empty()) continue;; 

            Field[1] = FieldInfo[1];
            Field[2] = FieldInfo[2];
            Field[18] = FieldInfo[18];
            if (_output_format == NDB_FILE_FORMAT_PDB) {
                 Field[5] = FieldInfo[17];
                 Field[6] = FieldInfo[15];
                 if (Field[0] == "ANISOU" || Field[0] == "SIGUIJ") {
                      Field[14] = FieldInfo[13];
                      Field[15] = FieldInfo[14];
                 }
            }

            const ndb_token_format& ndbformat1 = NdbToken::getTokenFormat(_aux_atoms[i]);
            writeGeneralRecord(lines, ndbformat1, Field);
            for (std::vector<std::string>::const_iterator vpos = lines.begin(); vpos != lines.end(); ++vpos) {
                 records.push_back(*vpos);
            }
       }

       if (_segid_flag) {
            std::string segid = FormattedFieldValue(FieldInfo[5], 3, 4, 0, true, true);
            for (std::vector<std::string>::iterator
                 vpos = records.begin(); vpos != records.end(); ++vpos) {
                 vpos->replace(20, 2, "  ");
                 vpos->replace(72, 4, segid);
            }
       }
}

void PdbWrite::WriteAtom(FILE* fp, const std::string& chain_type, RCSB::Atom* atom, std::vector<std::string>& FieldInfo, const bool include_auxiliary_atoms)
{
       std::vector<std::string> lines;
       WriteAtom(lines, chain_type, atom, FieldInfo, include_auxiliary_atoms);
       for (std::vector<std::string>::const_iterator vpos = lines.begin(); vpos != lines.end(); ++vpos) {
            fprintf(fp, "%s\n", vpos->c_str());
       }
}

void PdbWrite::ResetAtomCounting()
{
       _num_atom = 0;
       _num_ter = 0;
}

void PdbWrite::writeMasterRecord(FILE* fp)
{
       _getMasterRecord();

       std::map<std::string, std::list<std::vector<std::string> > >::iterator pos = _pdb_records->find("MASTER");
       if (pos != _pdb_records->end()) {
            const ndb_token_format& ndbformat = NdbToken::getTokenFormat("MASTER");
            writeGeneralRecords(fp, ndbformat, pos->second);
       }

       std::vector<std::string> FieldInfo;
       FieldInfo.clear();
       FieldInfo.push_back("END");
       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("END");
       writeGeneralRecord(fp, ndbformat, FieldInfo);
}

void PdbWrite::getAffectedFieldSet(const ndb_token_format& ndbformat, const std::list<std::vector<std::string> >& records, std::set<std::string>& fieldSet)
{
       fieldSet.clear();

       if (records.empty()) return;

       for (int i = 1; i < ndbformat.NumField; ++i) {
            std::map<std::string, std::pair<int, int> >::const_iterator mpos = _field_length_mapping.find(ndbformat.FieldList[i].FieldLabel);
            if (mpos == _field_length_mapping.end()) continue;

            for (std::list<std::vector<std::string> >::const_iterator lpos = records.begin(); lpos != records.end(); ++lpos) {
                 if ((int) (*lpos)[i].size() >= mpos->second.first) {
                      fieldSet.insert(ndbformat.FieldList[i].FieldLabel);
                      break;
                 }
            }
       }
}

void PdbWrite::updateNdbTokenFormat(const std::set<std::string>& fieldSet, ndb_token_format& ndbformat)
{
       if (fieldSet.empty()) return;

       int off_set = 0;
       for (int i = 1; i < ndbformat.NumField; ++i) {
            ndbformat.FieldList[i].FieldStcol += off_set;

            if (fieldSet.find(ndbformat.FieldList[i].FieldLabel) == fieldSet.end()) continue;
            std::map<std::string, std::pair<int, int> >::const_iterator mpos = _field_length_mapping.find(ndbformat.FieldList[i].FieldLabel);
            if (mpos == _field_length_mapping.end()) continue;
            off_set += mpos->second.second - ndbformat.FieldList[i].FieldWidth;
            ndbformat.FieldList[i].FieldWidth = mpos->second.second;
       }
}

// private

void PdbWrite::_init()
{
       _logIo = NULL;
       _ccDic = NULL;
       _pdb_records = NULL;
       _molecules = NULL;
       _segid_flag = false;
       _extended_flag = false;
       _original_entry_ids.clear();
       _atom_connects.clear();
       _filename.clear();
       _output_format = NDB_FILE_FORMAT_PDB;
       _num_remark = 0;
       _num_atom = 0;
       _num_ter = 0;

       _atomField_mapping.clear();

       // define corresponding atom/residue fields for LINK, SLTBRG, HYDBND & SHEET records
       std::map<int, int> tmp_map;
       tmp_map.clear();
       tmp_map.insert(std::make_pair(1, 3));
       tmp_map.insert(std::make_pair(7, 9));
       _atomField_mapping.insert(std::make_pair("LINK", tmp_map));
       _atomField_mapping.insert(std::make_pair("SLTBRG", tmp_map));

       tmp_map.clear();
       tmp_map.insert(std::make_pair(1, 3));
       tmp_map.insert(std::make_pair(7, 3));
       tmp_map.insert(std::make_pair(12, 14));
       _atomField_mapping.insert(std::make_pair("HYDBND", tmp_map));

       tmp_map.clear();
       tmp_map.insert(std::make_pair(13, 14));
       tmp_map.insert(std::make_pair(18, 19));
       _atomField_mapping.insert(std::make_pair("SHEET", tmp_map));

       // key: field name
       // pair.first: field length for checking
       // pair.second: field length for changing
       _field_length_mapping.clear();

       for (int i = 0; i < NUM_RESIDUE_NAME_LABEL; ++i) {
            _field_length_mapping.insert(std::make_pair(_residue_name_labels[i], std::make_pair(5, 5)));
       }

       for (int i = 0; i < NUM_CHAIN_ID_LABEL; ++i) {
            _field_length_mapping.insert(std::make_pair(_chain_id_labels[i], std::make_pair(2, 3)));
       }

       _affected_field_set.clear();
}

void PdbWrite::_write()
{
       if (_filename.empty()) return;
       FILE *fp = fopen(_filename.c_str(), "w");
       _write(fp);
       fclose (fp);
}

void PdbWrite::_write(FILE* fp, const int& header_section_only)
{
       if (_output_format == NDB_FILE_FORMAT_PDB && !header_section_only) NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_PDB);

       const std::vector<std::string>& ndbTokens = NdbToken::getNdbTokens();
       for (std::vector<std::string>::const_iterator tpos = ndbTokens.begin(); tpos != ndbTokens.end(); ++tpos) {
            if (*tpos == "ATOMN") {
                 if (!header_section_only) WriteCoordinates(fp);
                 continue;
            } else if (*tpos == "JRNL") {
                 _writeJrnlRecords(fp);
                 continue;
            }

            if (!_pdb_records) continue;

            std::map<std::string, std::list<std::vector<std::string> > >::iterator pos = _pdb_records->find(*tpos);
            if (pos == _pdb_records->end()) continue;

            const ndb_token_format& ndbformat = NdbToken::getTokenFormat(*tpos);
            writeGeneralRecords(fp, ndbformat, pos->second);
       }

       if (_output_format == NDB_FILE_FORMAT_PDB) {
            // print_master_card(fp, numRemark);
       }

       if (!header_section_only) {
            std::vector<std::string> FieldInfo;
            FieldInfo.clear();
            FieldInfo.push_back("END");
            const ndb_token_format& ndbformat = NdbToken::getTokenFormat("END");
            writeGeneralRecord(fp, ndbformat, FieldInfo);
       }

       if (_output_format == NDB_FILE_FORMAT_PDB && !header_section_only) NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_NDB);
}

void PdbWrite::_writeGeneralRecord(std::vector<std::string>& lines, const ndb_token_format& ndbformat, std::vector<std::string>& record,
                                  std::string& prevSeqNo, int& PdbContNo, const bool& full_record_flag)
{
       if (_output_format == NDB_FILE_FORMAT_PDB)
            _updateUpperCase(ndbformat, record);

       std::vector<std::string> continuedField;
       continuedField.clear();
       if (ndbformat.ContinuedF) {
            int i = ndbformat.ContinuedF - 1;
            if (ndbformat.TokenName == "REMARK")
                 continuedField.push_back(record[i]);
            else if (NdbToken::IsJrnlToken(ndbformat.TokenName))
                 _getJrnlContinuedField(ndbformat.TokenName, ndbformat.FieldList[i].FieldWidth, record[i], continuedField);
            else _getContinuedField(ndbformat.FieldList[i].FieldWidth, PdbContNo, record[i], continuedField);
       } else continuedField.push_back("");

       int start_field = 2;
       if (full_record_flag) start_field = 1;
       int ContNo = 0;
       for (unsigned int i = 0; i < continuedField.size(); ++i) {
            ContNo++;

            std::string temp_line = "";
            if (full_record_flag) temp_line = FormattedString(record[0], 6, true);

            for (int j = start_field; j < ndbformat.NumField; ++j) {
                 if (j >= (int) record.size()) break;
                 if (_output_format == NDB_FILE_FORMAT_PDB) {
                      if (ndbformat.FieldList[j].NdbOnly) break;
                      if (ndbformat.TokenName == "REF" && i > 0 && j == ndbformat.ContinuedF) break;
                 }

                 if (j == start_field && !full_record_flag)
                      temp_line += " ";
                 else temp_line += _printSpaceBetweenField(ndbformat, j);

                 temp_line += _printFieldValue(ndbformat, j, record, continuedField[i], ContNo, PdbContNo);
                 if (ndbformat.TokenName == "REMARK" &&
                     j == ndbformat.SeqField - 1) {
                      if (prevSeqNo != record[j] && !record[ndbformat.ContinuedF - 1].empty()) {
                           std::string temp_line1 = temp_line;
                           if (_output_format == NDB_FILE_FORMAT_PDB && full_record_flag)
                                _setEightyCharacters(temp_line1);
                           _num_remark++;
                           lines.push_back(temp_line1);
                      }
                      prevSeqNo = record[j];
                 }
            }

            if (_output_format == NDB_FILE_FORMAT_PDB && full_record_flag)
                 _setEightyCharacters(temp_line);

            lines.push_back(temp_line);
            if (ndbformat.TokenName == "REMARK") _num_remark++;

            if (ndbformat.PdbConField) PdbContNo++;
       }
}

void PdbWrite::_write_ter_card(std::string& text, std::vector<std::string>& FieldInfo)
{
       FieldInfo[0] = "TER";
       FieldInfo[2] = FieldInfo[4];
       FieldInfo[3] = FieldInfo[5];
       if (_segid_flag) FieldInfo[3].clear();
       FieldInfo[4] = FieldInfo[6];
       FieldInfo[5] = FieldInfo[7];
       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("TER");
       std::vector<std::string> lines;
       writeGeneralRecord(lines, ndbformat, FieldInfo);
       for (std::vector<std::string>::const_iterator vpos = lines.begin(); vpos != lines.end(); ++vpos) {
            text += *vpos + "\n";
       }
}

void PdbWrite::_write_ter_card(FILE* fp, std::vector<std::string>& FieldInfo)
{
       FieldInfo[0] = "TER";
       FieldInfo[2] = FieldInfo[4];
       FieldInfo[3] = FieldInfo[5];
       if (_segid_flag) FieldInfo[3].clear();
       FieldInfo[4] = FieldInfo[6];
       FieldInfo[5] = FieldInfo[7];
       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("TER");
       writeGeneralRecord(fp, ndbformat, FieldInfo);
}

void PdbWrite::_writeJrnlRecords(std::string& text)
{
       std::set<std::string> index_set;
       _getJrnlIndexSet(index_set);
       if (index_set.empty()) return;

       std::vector<std::string> index_array, lines;
       index_array.clear();
       if (index_set.find("primary") != index_set.end())
            index_array.push_back("primary");
       if (_output_format == NDB_FILE_FORMAT_PDB) {
            _writeSelectedJrnlRecords(lines, index_array);
            for (std::vector<std::string>::const_iterator vpos = lines.begin(); vpos != lines.end(); ++vpos) {
                 text += *vpos + "\n";
            }
            index_array.clear();
            _getSecondCitationIndex(index_array, index_set, ORIGINAL_CITATION);
            _writeSelectedJrnlRemarkRecords("0", index_array);
            index_array.clear();
            _getSecondCitationIndex(index_array, index_set, "");
            _writeSelectedJrnlRemarkRecords("1", index_array);
       } else {
            _getSecondCitationIndex(index_array, index_set, ORIGINAL_CITATION);
            _getSecondCitationIndex(index_array, index_set, "");
            _writeSelectedJrnlRecords(lines, index_array);
            for (std::vector<std::string>::const_iterator vpos = lines.begin(); vpos != lines.end(); ++vpos) {
                 text += *vpos + "\n";
            }
       }
}

void PdbWrite::_writeJrnlRecords(FILE *fp)
{
       std::set<std::string> index_set;
       _getJrnlIndexSet(index_set);
       if (index_set.empty()) return;

       std::vector<std::string> index_array;
       index_array.clear();
       if (index_set.find("primary") != index_set.end())
            index_array.push_back("primary");
       if (_output_format == NDB_FILE_FORMAT_PDB) {
             _writeSelectedJrnlRecords(fp, index_array);
             index_array.clear();
             _getSecondCitationIndex(index_array, index_set, ORIGINAL_CITATION);
             _writeSelectedJrnlRemarkRecords("0", index_array);
             index_array.clear();
             _getSecondCitationIndex(index_array, index_set, "");
             _writeSelectedJrnlRemarkRecords("1", index_array);
       } else {
            _getSecondCitationIndex(index_array, index_set, ORIGINAL_CITATION);
            _getSecondCitationIndex(index_array, index_set, "");
            _writeSelectedJrnlRecords(fp, index_array);
       }
}

void PdbWrite::_updateUpperCase(const ndb_token_format& ndbformat, std::vector<std::string>& record)
{
       for (int i = 1; i < ndbformat.NumField; ++i) {
            if (i >= (int) record.size()) break;
            if (ndbformat.FieldList[i].NdbOnly) break;
            if (record[i].empty()) continue;
            if (!ndbformat.FieldList[i].MixedCase || record[0] == "AUTHOR")
                 String::UpperCase(record[i]);
       }
}

void PdbWrite::_getJrnlContinuedField(const std::string& TokenName, const int max_len, const std::string& value, std::vector<std::string> &field)
{
       field.clear();
       if ((int) value.size() <= max_len) {
            field.push_back(value);
            return;
       }

       std::string cs = value;
       while (!cs.empty()) {
            String::StripAndCompressWs(cs);
            if ((int) cs.size() <= max_len) {
                 field.push_back(cs);
                 break;
            }
            int space_break = 0;
            int comma_break = 0;
            int dash_break = 0;
            int dot_break = 0;
            int right_parenthese = 0;
            for (int i = 0; i < (int) cs.size(); ++i) {
                 if (i == max_len) break;
                 if (cs[i] == ' ') space_break = i;
                 else if (cs[i] == ',') comma_break = i;
                 else if (cs[i] == '-') dash_break = i;
                 else if (cs[i] == '.') dot_break = i;
                 else if (cs[i] == ')') right_parenthese = i;  
            }
            int prev = 0;
            if (TokenName == "AUTH") {
                 if (comma_break > prev) prev = comma_break;
                 if (right_parenthese > prev) prev = right_parenthese;
                 if (prev < max_len / 5 && dash_break > prev) prev = dash_break;
                 if (prev < max_len / 5 && space_break > prev) prev = space_break;
            } else {
                if (space_break > prev) prev = space_break;
                if (prev < max_len / 5) {
                     if (comma_break > prev) prev = comma_break;
                     if (right_parenthese > prev) prev = right_parenthese;
                     if (TokenName == "REF" && dot_break > prev) prev = dot_break;
                }
            }
            if (prev < max_len / 5) prev = max_len;
            else prev++;
            if (TokenName == "DOI" || (TokenName != "AUTH" &&
               (int) cs.size() &&  cs[max_len] == ' '))
                prev = max_len;
            std::string cs1 = cs.substr(0, prev);
            String::StripAndCompressWs(cs1);
            field.push_back(cs1);
            if (prev < (int) cs.size())
                 cs = cs.substr(prev);
            else cs.clear();
       }
}

void PdbWrite::_getContinuedField(const int default_max_len, const int PdbContNo, const std::string& value, std::vector<std::string> &field)
{
       field.clear();
       int max_len = default_max_len;
       if (PdbContNo > 1) max_len = default_max_len - 1;
       if ((int) value.size() <= max_len) {
            field.push_back(value);
            return;
       }

       bool first = true;
       std::string cs = value;
       while (!cs.empty()) {
            String::StripAndCompressWs(cs);
            if ((int) cs.size() <= max_len) {
                 field.push_back(cs);
                 break;
            }
            int prev = 0;
            for (int i = 0; i < (int) cs.size(); ++i) {
                 if (i == max_len) break;
                 if (cs[i] == ' ' || cs[i] == ',' || cs[i] == '-') {
                      prev = i;
                 } else if (cs[i] == ')' && (i + 1) < (int) cs.size() && (i + 1) < max_len && cs[i+1] != ';') {
                      prev = i;
                 }
            }
            if (prev < max_len / 5) prev = max_len;
            else prev++;
            if (max_len < (int) cs.size() && cs[max_len] == ' ') prev = max_len;
            std::string cs1 = cs.substr(0, prev);
            String::StripAndCompressWs(cs1);
            field.push_back(cs1);
            if (prev < (int) cs.size())
                 cs = cs.substr(prev);
            else cs.clear();
            if (first && max_len == default_max_len) max_len = default_max_len - 1;
            first = false;
       }
}

std::string PdbWrite::_printSpaceBetweenField(const ndb_token_format& ndbformat, const int field_no)
{
       int width = ndbformat.FieldList[field_no].FieldStcol - ndbformat.FieldList[field_no-1].FieldStcol - ndbformat.FieldList[field_no-1].FieldWidth;
       return PrintSpace(width);
}

std::string PdbWrite::_printFieldValue(const ndb_token_format& ndbformat, const int& field_no, const std::vector<std::string>& Fields,
                                       const std::string& continuedField, const int& ContNo, const int& PdbContNo)
{

       std::string value = "";
       if (field_no == ndbformat.ConField - 1) {
            if (ContNo != 1) value = String::IntToString(ContNo);
       } else if (field_no == ndbformat.PdbConField - 1) {
            if (PdbContNo != 1) value = String::IntToString(PdbContNo);
       } else if (field_no == ndbformat.ContinuedF - 1) {
            if (PdbContNo != 1)
                 value = " " + continuedField;
            else value = continuedField;
       } else if (ndbformat.FieldList[field_no].FieldLabel.find("Atom_Name") != std::string::npos
               || ndbformat.FieldList[field_no].FieldLabel.find("atom_id") != std::string::npos) {
            int resField = _getResField(ndbformat.TokenName, field_no);
            if (resField)
                 value = printAtomNameField(_ccDic, "", Fields[field_no],  Fields[resField]);
            else value = Fields[field_no];
       } else {
            value = Fields[field_no];
            if (!value.empty() && (ndbformat.FieldList[field_no].FieldType <= 2)) {
                 if (!String::IsNumber(value)) {
                      value.clear();
                      std::string error = "PdbWrite::_printFieldValue: Token " + ndbformat.TokenName + " Field "
                                        + String::IntToString(field_no + 1) + " suppose to be ";
                      if (ndbformat.FieldList[field_no].FieldType == 1)
                           error += "integer,";
                      else error += "float,";
                      error += " but got " + Fields[field_no] + ".\n";
                      _logIo->message(error.c_str());
                 }
            }
       }
       return _FormattedFieldValue(value, ndbformat.FieldList[field_no]);
}

std::string PdbWrite::_FormattedFieldValue(const std::string& value, const field_format& FieldList)
{
       bool left_adjust = false;
       if (FieldList.FieldType > 2) left_adjust = (FieldList.FieldType % 2 == 1);
       return FormattedFieldValue(value, FieldList.FieldType, FieldList.FieldWidth, FieldList.FieldPrec, left_adjust, true);
}

void PdbWrite::_setEightyCharacters(std::string& line)
{
       int str_len = 80 - line.size();
       if (str_len > 0) line += PrintSpace(str_len); 
}

void PdbWrite::_getJrnlIndexSet(std::set<std::string>& index_set)
{
       index_set.clear();
       if (!_pdb_records) return;

       const std::vector<std::string>& jrnltokens = NdbToken::getJrnlTokens();
       for (unsigned int i = 0; i < jrnltokens.size(); ++i) {
            std::map<std::string, std::list<std::vector<std::string> > >::iterator pos = _pdb_records->find(jrnltokens[i]);
            if (pos == _pdb_records->end()) continue;

            const ndb_token_format& jrnlformat = NdbToken::getTokenFormat(jrnltokens[i]);
            int field_no = jrnlformat.SeqField - 1;
            for (std::list<std::vector<std::string> >::iterator lpos = pos->second.begin(); lpos != pos->second.end(); ++lpos) {
                 index_set.insert((*lpos)[field_no]);
            }
       }
}

void PdbWrite::_getSecondCitationIndex(std::vector<std::string>& index_array, std::set<std::string>& index_set, const std::string& prefix)
{
       for (unsigned int i = 0; i < index_set.size(); ++i) {
            std::string id = prefix + String::IntToString(i + 1);
            if (index_set.find(id) == index_set.end()) continue;
            index_array.push_back(id);
       }
}

void PdbWrite::_writeSelectedJrnlRecords(std::vector<std::string>& lines, std::vector<std::string>& index_array)
{
       lines.clear();
       if (index_array.empty()) return;

       const std::vector<std::string>& jrnltokens = NdbToken::getJrnlTokens();

       std::string prevSeqNo = "-1";
       int PdbContNo = 1;
       for (std::vector<std::string>::iterator vpos = index_array.begin(); vpos != index_array.end(); ++vpos) {
            for (unsigned int i = 0; i < jrnltokens.size(); ++i) {
                 std::map<std::string, std::list<std::vector<std::string> > >::iterator pos = _pdb_records->find(jrnltokens[i]);
                 if (pos == _pdb_records->end()) continue;
     
                 const ndb_token_format& jrnlformat = NdbToken::getTokenFormat(jrnltokens[i]);
                 int field_no = jrnlformat.SeqField - 1;
                 for (std::list<std::vector<std::string> >::iterator lpos = pos->second.begin(); lpos != pos->second.end(); ++lpos) {
                      if ((*lpos)[field_no] != *vpos) continue;
                      if (_output_format == NDB_FILE_FORMAT_PDB && (*lpos)[field_no] == "primary") (*lpos)[field_no].clear();
                      _writeGeneralRecord(lines, jrnlformat, *lpos, prevSeqNo, PdbContNo);
                      break;
                 }
            }
       }
}

void PdbWrite::_writeSelectedJrnlRecords(FILE *fp, std::vector<std::string>& index_array)
{
       std::vector<std::string> lines;
       _writeSelectedJrnlRecords(lines, index_array);

       for (std::vector<std::string>::const_iterator vpos = lines.begin(); vpos != lines.end(); ++vpos) {
            fprintf(fp, "%s\n", vpos->c_str());
       }
}

void PdbWrite::_writeSelectedJrnlRemarkRecords(const std::string& remark_no, std::vector<std::string>& index_array)
{
       if (index_array.empty()) return;

       const std::vector<std::string>& jrnltokens = NdbToken::getJrnlTokens();

       std::vector<std::string> lines, data;
       data.clear();
       data.push_back("REMARK");
       data.push_back(remark_no);
       data.push_back("");
       std::list<std::vector<std::string> > remark_list;
       remark_list.clear();

       std::string prevSeqNo = "-1";
       int PdbContNo = 1;
       int count = 0;
       for (std::vector<std::string>::iterator vpos = index_array.begin(); vpos != index_array.end(); ++vpos) {
            lines.clear();
            for (unsigned int i = 0; i < jrnltokens.size(); ++i) {
                 std::map<std::string, std::list<std::vector<std::string> > >::iterator pos = _pdb_records->find(jrnltokens[i]);
                 if (pos == _pdb_records->end()) continue;
     
                 const ndb_token_format& jrnlformat = NdbToken::getTokenFormat(jrnltokens[i]);
                 int field_no = jrnlformat.SeqField - 1;
                 for (std::list<std::vector<std::string> >::iterator lpos = pos->second.begin(); lpos != pos->second.end(); ++lpos) {
                      if ((*lpos)[field_no] != *vpos) continue;
                      _writeGeneralRecord(lines, jrnlformat, *lpos, prevSeqNo, PdbContNo, false);
                      break;
                 }
            }
            if (lines.empty()) continue;

            if (remark_no == "0") {
                 data[2] = "ORIGINAL DATA REFERENCE " + String::IntToString(count + 1);
                 remark_list.push_back(data);
                 if (!_original_entry_ids.empty()) {
                      if (count < (int) _original_entry_ids.size())
                           data[2] = " PDB ID: " + _original_entry_ids[count];
                      else data[2] = " PDB ID: " + _original_entry_ids[0];
                      remark_list.push_back(data);
                 }
            } else {
                 data[2] = "REFERENCE " + String::IntToString(count + 1);
                 remark_list.push_back(data);
            }
            for (std::vector<std::string>::const_iterator lpos = lines.begin(); lpos != lines.end(); ++lpos) {
                 data[2] = *lpos;
                 remark_list.push_back(data);
            }
            count++;
       }

       if (remark_list.empty()) return;

       std::map<std::string, std::list<std::vector<std::string> > >::iterator pos = _pdb_records->find("REMARK");
       if (pos == _pdb_records->end())
            _pdb_records->insert(std::make_pair("REMARK", remark_list));
       else if (pos->second.empty())
            pos->second = remark_list;
       else {
            std::list<std::vector<std::string> >::iterator lpos; for (lpos = pos->second.begin(); lpos != pos->second.end(); ++lpos) {
                 if (atoi((*lpos)[1].c_str()) > atoi(remark_no.c_str())) break;
            }
            pos->second.splice(lpos, remark_list);
       }
}

int PdbWrite::_getResField(const std::string& TokenName, const int& field_no)
{
       std::map<std::string, std::map<int, int> >::const_iterator mpos = _atomField_mapping.find(TokenName);
       if (mpos == _atomField_mapping.end()) return 0;

       std::map<int, int>::const_iterator ipos = mpos->second.find(field_no);
       if (ipos != mpos->second.end()) return ipos->second;

       return 0;
}

void PdbWrite::_getConnectRecord()
{
       if (_atom_connects.empty()) return;

       if (!_pdb_records) return;

       std::map<int, std::set<int> > _connect;
       _connect.clear();
       for (std::list<std::pair<RCSB::Atom*, RCSB::Atom*> >::const_iterator pos = _atom_connects.begin(); pos != _atom_connects.end(); ++pos) {
            int atom_id1 = atoi(pos->first->atnum().c_str());
            if (atom_id1 > 99999) continue;
            int atom_id2 = atoi(pos->second->atnum().c_str());
            if (atom_id2 > 99999) continue;
            _insert_connect(_connect, atom_id1, atom_id2);
            _insert_connect(_connect, atom_id2, atom_id1);
       }

       if (_connect.empty()) return;

       std::vector<std::string> data;
       data.clear();
       for (int i = 0; i < 12; ++i) data.push_back("");
       data[0] = "CONECT";

       std::list<std::vector<std::string> > connect_list;
       connect_list.clear();
       for (std::map<int, std::set<int> >::const_iterator mpos = _connect.begin(); mpos != _connect.end(); ++mpos) {
            for (int i = 1; i < 12; ++i) data[i].clear();
            data[1] = String::IntToString(mpos->first);
            int field = 0;
            for (std::set<int>::const_iterator spos = mpos->second.begin(); spos != mpos->second.end(); ++spos) {
                 if (field == 4) {
                      connect_list.push_back(data);
                      for (int i = 2; i < 12; ++i) data[i].clear();
                      field = 0;
                 }
                 data[2 + field] = String::IntToString(*spos);
                 field++;
            }
            if (field) connect_list.push_back(data);
       }

       if (connect_list.empty()) return;

       std::map<std::string, std::list<std::vector<std::string> > >::iterator mlpos = _pdb_records->find("CONECT");
       if (mlpos != _pdb_records->end())
            mlpos->second = connect_list;
       else _pdb_records->insert(std::make_pair("CONECT", connect_list));
}

void PdbWrite::_insert_connect(std::map<int, std::set<int> >& connect, const int& atom_id1, const int& atom_id2)
{
       std::map<int, std::set<int> >::iterator mpos = connect.find(atom_id1);
       if (mpos != connect.end()) {
            mpos->second.insert(atom_id2);
            return;
       }

       std::set<int> t_set;
       t_set.clear();
       t_set.insert(atom_id2);
       connect.insert(std::make_pair(atom_id1, t_set));
}

void PdbWrite::_getMasterRecord()
{
       std::vector<std::string> data, tokens;
       data.clear();
       for (int i = 0; i < 13; ++i) data.push_back("");
       data[0] = "MASTER";

       tokens.clear();
       tokens.push_back("");
       tokens.push_back("REMARK");
       tokens.push_back("FTNOTE");
       tokens.push_back("HET");
       tokens.push_back("HELIX");
       tokens.push_back("SHEET");
       tokens.push_back("TURN");
       tokens.push_back("SITE");
       tokens.push_back("");
       tokens.push_back("");
       tokens.push_back("");
       tokens.push_back("CONECT");
       tokens.push_back("SEQRES");
       for (unsigned int i = 1; i < tokens.size(); ++i) {
            int count = 0;
            if (!tokens[i].empty()) {
                 std::map<std::string, std::list<std::vector<std::string> > >::const_iterator mpos = _pdb_records->find(tokens[i]);
                 if (mpos != _pdb_records->end()) count = mpos->second.size();
            }
            data[i] = String::IntToString(count);
       }

       data[1] = String::IntToString(_num_remark);

       tokens.clear();
       tokens.push_back("ORIGX");
       tokens.push_back("SCALE");
       tokens.push_back("MTRIX");
       int count = 0;
       for (std::vector<std::string>::const_iterator pos = tokens.begin(); pos != tokens.end(); ++pos) {
            std::map<std::string, std::list<std::vector<std::string> > >::const_iterator mpos = _pdb_records->find(*pos);
            if (mpos == _pdb_records->end()) continue;
            count += mpos->second.size();
       }
       data[8] = String::IntToString(count);

/*
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator mpos = _pdb_records->find("ATNUMS");
       if (mpos != _pdb_records->end()) {
            const std::vector<std::string>& Field = mpos->second.front();
            data[9] = Field[6];
       }
*/
       data[9] = String::IntToString(_num_atom);

       data[10] = String::IntToString(_num_ter);

       _pdb_records->erase("MASTER");
       std::list<std::vector<std::string> > m_list;
       m_list.clear();
       m_list.push_back(data);
       _pdb_records->insert(std::make_pair("MASTER", m_list));
}

std::string printAtomNameField(ConnectDic *ccdic, const std::string& atomtype, const std::string& atomName, const std::string& residueName)
{
       if (atomName.empty()) return "    ";
       else if (atomName.size() == 4) return atomName;
       if (residueName == "HOH" && atomName[0] == 'O') return " O  ";

       std::string atom_type = atomtype;
       if (atom_type.empty() && !residueName.empty())
            atom_type = ccdic->getAtomType(residueName, atomName);
       if (atom_type.empty())
            atom_type = Element::getAtomSymbol(atomName);

       int align = Element::getAtomAlignPosition(atomName, atom_type);

       std::string cs = atomName;
       if (align == 0) {
            if (atomName.size() == 1)
                 cs = atomName + "   ";
            else if (atomName.size() == 2)
                 cs = atomName + "  ";
            else if (atomName.size() == 3)
                 cs = atomName + " ";
       } else {
            if (atomName.size() == 1)
                 cs = " " + atomName + "  ";
            else if (atomName.size() == 2)
                 cs = " " + atomName + " ";
            else if (atomName.size() == 3)
                 cs = " " + atomName;
       }
       return cs;
}
