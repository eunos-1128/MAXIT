/*
FILE:     NdbToken.C
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

#include "CifFile.h"
#include "NdbToken.h"
#include "NdbToken_global.h"
#include "TypeDef.h"
#include "utillib.h"

void NdbToken::initialize()
{
       _all_tokens.clear();
       _token_mapping.clear();
       _ndb_tokens.clear();
       _ndb_token_set.clear();
       _jrnl_tokens.clear();
       _jrnl_token_set.clear();
       _atom_info.clear();
       _atom_token_set.clear();
       _maxNumField = 0;
}

bool NdbToken::Read(LogUtil& logIo, const std::string& path)
{
       std::string binaryfile = path + "/data/binary/";
       binaryfile += tokenbinfile;
       struct stat statbuf;
       if (stat(binaryfile.c_str(), &statbuf) != 0) {
            std::string error = "NdbToken::Read: Can't find "
                              + binaryfile + "\n";
            logIo.message(error.c_str());
            return false;
       }

       CifFile* fobjR = new CifFile(READ_MODE, binaryfile);
       if (!fobjR) {
            std::string error = "NdbToken::Read: Can't open "
                              + binaryfile + "\n";
            logIo.message(error.c_str());
            return false;
       }

       Block &block = fobjR->GetBlock(fobjR->GetFirstBlockName());

       bool found_error = false;
       if (!block.IsTablePresent("ndb_token")) {
            std::string error = "NdbToken::Read: Can't find ndb_token category in "
                              + binaryfile + "\n";
            logIo.message(error.c_str());
            found_error = true;
       }
       if (!block.IsTablePresent("ndb_field")) {
            std::string error = "NdbToken::Read: Can't find ndb_field category in "
                              + binaryfile + "\n";
            logIo.message(error.c_str());
            found_error = true;
       }
       if (!block.IsTablePresent("jrnl_token")) {
            std::string error = "NdbToken::Read: Can't find jrnl_token category in "
                              + binaryfile + "\n";
            logIo.message(error.c_str());
            found_error = true;
       }
       if (!block.IsTablePresent("jrnl_field")) {
            std::string error = "NdbToken::Read: Can't find jrnl_field category in "
                              + binaryfile + "\n";
            logIo.message(error.c_str());
            found_error = true;
       }
       if (!block.IsTablePresent("atom_token")) {
            std::string error = "NdbToken::Read: Can't find atom_token category in "
                              + binaryfile + "\n";
            logIo.message(error.c_str());
            found_error = true;
       }
       if (found_error) {
            delete fobjR; return false;
       }

       ISTable *Table = getTablePtr(block, "ndb_token");
       ISTable *t = getTablePtr(block, "jrnl_token");
       unsigned int num = Table->GetNumRows() + t->GetNumRows();
       _all_tokens.reserve(num);

       read_token_table(Table, _all_tokens, _ndb_tokens, _ndb_token_set);
       read_token_table(t, _all_tokens, _jrnl_tokens, _jrnl_token_set);
       
       num = 0;
       for (std::vector<ndb_token_format>::iterator
            pos = _all_tokens.begin(); pos != _all_tokens.end(); ++pos) {
            _token_mapping.insert(std::make_pair(pos->TokenName, num));
            num++;
       }

       Table = getTablePtr(block, "ndb_field");
       read_field_table(Table, _all_tokens, _token_mapping);

       Table = getTablePtr(block, "jrnl_field");
       read_field_table(Table, _all_tokens, _token_mapping);

       for (std::vector<ndb_token_format>::iterator
            pos = _all_tokens.begin(); pos != _all_tokens.end(); ++pos) {
            pos->NumField = pos->FieldList.size();
            if (pos->NumField > _maxNumField) _maxNumField = pos->NumField;
            num = 0;
            for (std::vector<field_format>::iterator ppos =
                 pos->FieldList.begin(); ppos != pos->FieldList.end(); ++ppos) {
                 num++;
                 if (ppos->FieldLabel == "Refine_ID") {
                      pos->RefineField = num;
                      break;
                 }
            }
       }

       Table = getTablePtr(block, "atom_token");
       read_atom_table(Table, _atom_info);

       delete fobjR;

       return true;
}

const int& NdbToken::MaxNumField() { return _maxNumField; }

void NdbToken::resetAtomTokenField(const int& format, const bool& strict_pdb_format_chainid_field_flag)
{
       std::set<std::string> unique_atom_token_set;
       unique_atom_token_set.clear();

       for (int i = 0; i < NUM_ATOM_CARD; ++i) {
            std::map<std::string, unsigned int>::iterator pos = _token_mapping.find(_atom_tokens[i].token);
            if (pos == _token_mapping.end()) continue;
            int field_no = _atom_tokens[i].field_no - 1;
            if (format == NDB_FILE_FORMAT_PDB) {
                 _all_tokens[pos->second].FieldList[field_no].FieldStcol = _atom_tokens[i].new_field_start;
                 _all_tokens[pos->second].FieldList[field_no].FieldType  = _atom_tokens[i].new_field_type;
                 _all_tokens[pos->second].FieldList[field_no].FieldWidth = _atom_tokens[i].new_field_width;
                 _all_tokens[pos->second].FieldList[field_no].NdbOnly    = _atom_tokens[i].new_ndb_flag;
                 if (strict_pdb_format_chainid_field_flag && (unique_atom_token_set.find(_atom_tokens[i].token) == unique_atom_token_set.end())) {
                      unique_atom_token_set.insert(_atom_tokens[i].token);
                      // set strict PDB format definition for chainID field.
                      if ((_all_tokens[pos->second].FieldList[5].FieldStcol == 21) && (_all_tokens[pos->second].FieldList[5].FieldWidth == 2)) {
                           _all_tokens[pos->second].FieldList[5].FieldStcol = 22;
                           _all_tokens[pos->second].FieldList[5].FieldWidth = 1;
                      }
                 }
            } else if (format == NDB_FILE_FORMAT_NDB || format == NDB_FILE_FORMAT_MIX) {
                 _all_tokens[pos->second].FieldList[field_no].FieldStcol = _atom_tokens[i].old_field_start;
                 _all_tokens[pos->second].FieldList[field_no].FieldType  = _atom_tokens[i].old_field_type;
                 _all_tokens[pos->second].FieldList[field_no].FieldWidth = _atom_tokens[i].old_field_width;
                 _all_tokens[pos->second].FieldList[field_no].NdbOnly    = _atom_tokens[i].old_ndb_flag;
                 // set loose PDB format definition for chainID field ( defined in ndb_tokens.cif file )
                 if ((_all_tokens[pos->second].FieldList[5].FieldStcol == 22) && (_all_tokens[pos->second].FieldList[5].FieldWidth == 1)) {
                      _all_tokens[pos->second].FieldList[5].FieldStcol = 21;
                      _all_tokens[pos->second].FieldList[5].FieldWidth = 2;
                 }
            }
       }
}

bool NdbToken::IsValidToken(const std::string& token)
{
       if (token.empty()) return false;

       std::string cs = token;
       if (token == "ATOM" || token == "HETATM") cs = "ATOMN";
       if (_token_mapping.find(cs) != _token_mapping.end()) return true;
       int size = cs.size() - 1;
       if (isdigit(cs[size])) {
            if (_token_mapping.find(cs.substr(0, size)) != _token_mapping.end())
                 return true;
       }
       return false;
}

bool NdbToken::IsNdbToken(const std::string& token)
{
       if (token.empty()) return false;

       std::string cs = token;
       if (token == "ATOM" || token == "HETATM") cs = "ATOMN";
       if (_ndb_token_set.find(cs) != _ndb_token_set.end()) return true;
       int size = cs.size() - 1;
       if (isdigit(cs[size])) {
            if (_ndb_token_set.find(cs.substr(0, size)) != _ndb_token_set.end())
                 return true;
       }
       return false;
}

bool NdbToken::IsJrnlToken(const std::string& token)
{
       if (_jrnl_token_set.find(token) != _jrnl_token_set.end()) return true;
       return false;
}

bool NdbToken::IsAtomToken(const std::string& token)
{
       if (_atom_token_set.find(token) != _atom_token_set.end()) return true;
       return false;
}

bool NdbToken::IsMatrixToken(const std::string& token)
{
       if (token.empty()) return false;

       std::map<std::string, unsigned int>::const_iterator
           pos = _token_mapping.find(token);
       if (pos == _token_mapping.end()) {
            int size = token.size() - 1;
            if (isdigit(token[size])) pos = _token_mapping.find(token.substr(0, size));
       }
       if (pos != _token_mapping.end()) {
            if (_all_tokens[pos->second].MatrixType == 3) return true;
       }
       return false;
}

const ndb_token_format& NdbToken::getTokenFormat(const std::string& token)
{
       std::string cs = token;
       if (token == "ATOM" || token == "HETATM") cs = "ATOMN";
       std::map<std::string, unsigned int>::const_iterator
           pos = _token_mapping.find(cs);
       if (pos == _token_mapping.end()) {
            int size = cs.size() - 1;
            if (isdigit(cs[size])) pos = _token_mapping.find(cs.substr(0, size));
       }
       if (pos == _token_mapping.end()) {
            throw std::out_of_range("Token " + token + " is unrecognized.\n");
       }
       return _all_tokens[pos->second];
}

const std::vector<std::string>& NdbToken::getNdbTokens()
{
       return _ndb_tokens;
}

const std::vector<std::string>& NdbToken::getJrnlTokens()
{
       return _jrnl_tokens;
}

void NdbToken::read_token_table(ISTable *Table, std::vector<ndb_token_format> &token_list,
                                std::vector<std::string>& token_names,
                                std::set<std::string>& token_name_set)
{
       int rowNo = Table->GetNumRows();
       std::string cs;
       ndb_token_format ndbtoken;
       for (int row = 0; row < rowNo; ++row) {
            ndbtoken.TokenName.clear();
            ndbtoken.SeqField = 0;
            ndbtoken.NumField = 0;
            ndbtoken.ConField = 0;
            ndbtoken.PdbConField = 0;
            ndbtoken.ContinuedF = 0;
            ndbtoken.MatrixType = 0;
            ndbtoken.NdbOnly = 0;
            ndbtoken.DisplayOrder = row;
            ndbtoken.RefineField = 0;
            ndbtoken.FieldList.clear();

            get_value(cs, Table, row, "name");
            if (!cs.empty()) ndbtoken.TokenName = cs;
            get_value(cs, Table, row, "serial_no");
            if (!cs.empty()) ndbtoken.SeqField = atoi(cs.c_str());
            get_value(cs, Table, row, "continuation_no");
            if (!cs.empty()) ndbtoken.ConField = atoi(cs.c_str());
            get_value(cs, Table, row, "continuation");
            if (!cs.empty()) ndbtoken.ContinuedF = atoi(cs.c_str());
            get_value(cs, Table, row, "matrix_field");
            if (!cs.empty()) ndbtoken.MatrixType = atoi(cs.c_str());
            get_value(cs, Table, row, "ndb_only");
            if (!cs.empty()) ndbtoken.NdbOnly = atoi(cs.c_str());
            get_value(cs, Table, row, "pdb_continued_no_field");
            if (!cs.empty()) ndbtoken.PdbConField = atoi(cs.c_str());
            token_list.push_back(ndbtoken);
            token_names.push_back(ndbtoken.TokenName);
            token_name_set.insert(ndbtoken.TokenName);
       }
}

void NdbToken::read_field_table(ISTable *Table, std::vector<ndb_token_format> &token_list,
                std::map<std::string, unsigned int>& token_mapping)
{
       int rowNo = Table->GetNumRows();
       field_format field;
       std::string cs;
       for (int row = 0; row < rowNo; ++row) {
            get_value(cs, Table, row, "token_name");
            std::map<std::string, unsigned int>::iterator pos = token_mapping.find(cs);
            if (pos == token_mapping.end()) continue;

            field.FieldLabel.clear();
            field.FieldType = 0;
            field.FieldWidth = 0;
            field.FieldPrec = 0;
            field.FieldStcol = 0;
            field.NdbOnly = 0;
            field.MixedCase = 0;
            field.Mandatory = 0;

            get_value(cs, Table, row, "name");
            if (!cs.empty()) field.FieldLabel = cs;
            get_value(cs, Table, row, "type");
            if (!cs.empty()) field.FieldType = atoi(cs.c_str());
            get_value(cs, Table, row, "width");
            if (!cs.empty()) field.FieldWidth = atoi(cs.c_str());
            get_value(cs, Table, row, "precision");
            if (!cs.empty()) field.FieldPrec = atoi(cs.c_str());
            get_value(cs, Table, row, "start_column");
            if (!cs.empty()) field.FieldStcol = atoi(cs.c_str());
            get_value(cs, Table, row, "ndb_only");
            if (!cs.empty()) field.NdbOnly = atoi(cs.c_str());
            get_value(cs, Table, row, "mixed_cases");
            if (!cs.empty()) field.MixedCase = atoi(cs.c_str());
            get_value(cs, Table, row, "mandatory");
            if (!cs.empty()) field.Mandatory = atoi(cs.c_str());
            token_list[pos->second].FieldList.push_back(field);
       }
}

void NdbToken::read_atom_table(ISTable *Table, std::map<std::string,
                               atom_format>& _atom_info) 
{
       int rowNo = Table->GetNumRows();
       atom_format atom;
       std::string cs;
       for (int row = 0; row < rowNo; ++row) {
            atom.AuxAtom = 0;
            atom.Esd = 0;
            atom.SeqAtom = 0;
            atom.PdbTitle.clear();
            get_value(cs, Table, row, "regular_or_aux");
            if (!cs.empty()) atom.AuxAtom = atoi(cs.c_str());
            get_value(cs, Table, row, "esd");
            if (!cs.empty()) atom.Esd = atoi(cs.c_str());
            get_value(cs, Table, row, "serial_atom");
            if (!cs.empty()) atom.SeqAtom = atoi(cs.c_str());
            get_value(cs, Table, row, "pdb_token_name");
            if (!cs.empty()) atom.PdbTitle = cs;
            // get_value(cs, Table, row, "table_name");
            // if (!cs.empty()) atom.TableName = cs;
            get_value(cs, Table, row, "token_name");
            if (!cs.empty()) {
                 _atom_info.insert(std::make_pair(cs, atom));
                 _atom_token_set.insert(cs);
            }
       }
       _atom_token_set.insert("ATOM");
       _atom_token_set.insert("HETATM");
       _atom_token_set.insert("MODEL");
       _atom_token_set.insert("ENDMDL");
}
