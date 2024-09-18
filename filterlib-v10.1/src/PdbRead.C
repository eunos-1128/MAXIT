/*
FILE:     PdbRead.C
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
#include <ctype.h>

#include "PdbRead.h"
#include "PdbRead_global.h"
#include "TypeDef.h"
#include "utillib.h"

#define EXPERIMENT_TYPE_XRAY        2
#define EXPERIMENT_TYPE_NEUTRON     4
#define EXPERIMENT_TYPE_FIBER       8
#define EXPERIMENT_TYPE_ELECTRON   16

#define MINIMUM_ATOM_TOKEN_FIELD   11
#define MAXIMUM_ATOM_TOKEN_FIELD   13

PdbRead::PdbRead()
{
       _init();
       _init_token_set();
}

PdbRead::~PdbRead()
{
       _init();
}

void PdbRead::_init()
{
       _logIo = NULL;
       _messageIo = NULL;
       _pdb_records = NULL;
       _molecules = NULL;
       _filename.clear();
       _input_format = 0;
       _first_atom = true;
       _format_checking = false;
       _strict_pdb_format_chainid_field_flag = false;
       _no_pdb_chain_id_flag = false;
       _atom_token_field = MAXIMUM_ATOM_TOKEN_FIELD;
       _model_set.clear();
       _skip_token_set.clear();
       _check_token_set.clear();
       _atom_token_set.clear();
       _read_token_set.clear();
       _wrong_PDB_chain_ID.clear();
       _error_count.clear();
}

void PdbRead::_init_token_set()
{
       for (int i = 0; i < NUM_SKIP_TOKEN; ++i) { _skip_token_set.insert(_skip_tokens[i]); }
       for (int i = 0; i < NUM_CHECK_TOKEN; ++i) { _check_token_set.insert(_check_tokens[i]); }
       for (int i = 0; i < NUM_ATOM_TOKEN; ++i) { _atom_token_set.insert(_atom_tokens[i]); }
}

void PdbRead::setLog(LogUtil *logPt)
{
       _logIo = logPt;
}

void PdbRead::setMessage(MessageUtil *message)
{
       _messageIo = message;
}

void PdbRead::setRecord(std::map<std::string, std::list<std::vector<std::string> > >* records)
{
       _pdb_records = records;
}

void PdbRead::setMolecule(std::vector<RCSB::Molecule*>* mols)
{
       _molecules = mols;
}

void PdbRead::setFileName(const std::string& filename)
{
       _filename = filename;
}

void PdbRead::setInputFormat(const int& format)
{
       _input_format = format;
}

void PdbRead::setExperimentType(const int& type)
{
       if (type & EXPERIMENT_TYPE_XRAY || type & EXPERIMENT_TYPE_NEUTRON || type & EXPERIMENT_TYPE_FIBER || type & EXPERIMENT_TYPE_ELECTRON)
            _atom_token_field = MAXIMUM_ATOM_TOKEN_FIELD;
       else _atom_token_field = MINIMUM_ATOM_TOKEN_FIELD;
}

void PdbRead::setFormatChecking()
{
       _format_checking = true;
}

void PdbRead::setStrictPdbFormatFlag()
{
       _strict_pdb_format_chainid_field_flag = true;
}

bool PdbRead::getNoPdbRecordFlag()
{
       return _read_token_set.empty();
}

void PdbRead::ReadRecord()
{
       _read();
}

void PdbRead::ReadRecord(const std::string& filename, const int format, std::map<std::string, std::list<std::vector<std::string> > >* records,
                         std::vector<RCSB::Molecule*>* mols)
{
       _filename = filename;
       _input_format = format;
       _pdb_records = records;
       _molecules = mols;
       _read();
}

void PdbRead::ReadRecord(const std::string& coord)
{
       if (coord.empty()) return;

       std::vector<std::string> data;
       get_wordarray_with_space(data, coord, "\n");

       _pdb_records->clear();
       _molecules->clear();

       int MxField = NdbToken::MaxNumField() + 2;

       std::vector<std::string> FieldInfo;
       FieldInfo.reserve(MxField);
       FieldInfo.clear();

       if (_input_format == NDB_FILE_FORMAT_PDB) NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_PDB, _strict_pdb_format_chainid_field_flag);

       std::string token;
       int lineno = 0;
       for (std::vector<std::string>::const_iterator pos = data.begin(); pos != data.end(); ++pos) {
            lineno++;
            if (pos->empty()) continue;

            _get_token_from_line(token, *pos);

            if (token == "END") continue;

            if (!NdbToken::IsValidToken(token)) {
                 if (_skip_token_set.find(token) == _skip_token_set.end()) {
                      std::string error = "PdbRead::_read: The record " + token + " at line '" + *pos + "' is unrecognized.\n";
                      _logIo->message(error.c_str());
                 }
                 continue;
            }

            if (_get_one_line_info(FieldInfo, token, *pos)) {
                 if (_format_checking && (token != "REMARK")) _read_token_set.insert(token);
            }

            if (NdbToken::IsAtomToken(token)) 
                 _processAtomRecord(FieldInfo, lineno);
            else _processTextRecord(FieldInfo);
       }

       _removeLastSemicolon();
       _mergeConcatenateRecords();

       if (_molecules->size() == 1) (*_molecules)[0]->set_Mol_ID(1);

       if (_input_format == NDB_FILE_FORMAT_PDB) NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_NDB);
}

void PdbRead::_read()
{
       _pdb_records->clear();
       _molecules->clear();

       struct stat statbuf;
       if (stat(_filename.c_str(), &statbuf) != 0) {
            std::string error = "PdbRead::_read: File " + _filename + " does not exist.\n";
            _logIo->message(error.c_str());
            if (_messageIo && _format_checking) _messageIo->insertMessage("error", "model", "File " + _filename + " does not exist.", true);
            return;
       }

       FILE *infile = fopen(_filename.c_str(), "r");
       if (infile == (FILE *) NULL) {
            std::string error = "PdbRead::_read: Can not read file " + _filename + ".\n";
            _logIo->message(error.c_str());
            if (_messageIo && _format_checking) _messageIo->insertMessage("error", "model", "Can not read file " + _filename + ".", true);
            return;
       }

       int MxField = NdbToken::MaxNumField() + 2;

       std::vector<std::string> FieldInfo;
       FieldInfo.reserve(MxField);
       FieldInfo.clear();

       if (_input_format == NDB_FILE_FORMAT_PDB) NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_PDB, _strict_pdb_format_chainid_field_flag);

       std::string line, token, model_record, last_line;
       std::string prefix = "";
       int lineno = 0;
       int num_ter = 0;
       int num_end = 0;
       bool found_model = false;
       bool end_in_middle = false;
       bool first_remark_0 = true;
       bool first_remark_1 = true;
       bool start_remark_citation = false;
       bool found_junk_remark = false;
       bool not_found_TER_card_before_ENDMDL = false;
       int serial_no = 0;
       while (!feof(infile)) {
            // get_one_line(infile, line);
            get_line_from_file(infile, line);
            lineno++;
            if (line.empty()) continue;
            _get_token_from_line(token, line);

            if (token.empty() && _messageIo && _format_checking) {
                 _check_line_not_start_from_first_column(line, lineno);
                 continue;
            }

            if (token == "REMARK" && line.find("FILE IS FOR PDB DEPOSITION: REMOVE ALL FROM THIS LINE UP") != std::string::npos) {
                 found_junk_remark = true;
                 _pdb_records->erase("REMARK");
                 continue;
            }
            if (token == "END") {
                 num_end++;
                 _no_pdb_chain_id_flag = false;
                 continue;
            }

            if (!NdbToken::IsValidToken(token)) {
                 if (_skip_token_set.find(token) == _skip_token_set.end()) {
                      std::string error = "PdbRead::_read: The record " + token + " at line " + String::IntToString(lineno)
                                        + " is unrecognized.\n";
                      _logIo->message(error.c_str());
                      if (_messageIo && _format_checking) {
                           error = "Token '" + token + "' is unrecognized in line " + String::IntToString(lineno)
                                 + ". Please refer to the <a href=\"http://www.wwpdb.org/documentation/format33/v3.3.html\">";
                           error += "PDB V3.3 format</a> for the proper token name.";
                           _messageIo->insertMessage("error", "model", error, true);
                      }
                 }
                 continue;
            }

            if (token == "ATOM" || token == "HETATM") {
                 if (_messageIo && _format_checking) last_line = line;
            } else if (token == "MODEL") {
                 if (found_model && _messageIo && _format_checking) _write_no_ENDMDL_error(model_record);
                 found_model = true;
                 num_ter = 0;
                 model_record = line;
                 _no_pdb_chain_id_flag = false;
            } else if (token == "ENDMDL") {
                 found_model = false;
                 if (_no_pdb_chain_id_flag && _messageIo && _format_checking) {
                      std::string error = "Polymer atom records are missing the chain ID field before line " + String::IntToString(lineno)
                                        + ". Please add the chain ID to column 22 of the PDB file and start again.";
                      _messageIo->insertMessage("error", "model", error, true);
                 }
                 _no_pdb_chain_id_flag = false;
                 if (num_ter == 0) not_found_TER_card_before_ENDMDL = true;
            } else if (token == "TER") {
                 num_ter++;
                 if (_no_pdb_chain_id_flag && _messageIo && _format_checking) {
                      std::string error = "Polymer atom records are missing the chain ID field before line " + String::IntToString(lineno)
                                        + ". Please add the chain ID to column 22 of the PDB file and start again.";
                      _messageIo->insertMessage("error", "model", error, true);
                 }
                 _no_pdb_chain_id_flag = false;
            }

            if (token != "END" && num_end) end_in_middle = true;

            bool format_error = _get_one_line_info(FieldInfo, token, line);
            if (!format_error && _format_checking && (token != "REMARK")) _read_token_set.insert(token);

            if (format_error && _messageIo && _format_checking) _write_record_error(token, line, lineno);

            if (token == "REMARK" && _input_format == NDB_FILE_FORMAT_PDB && found_junk_remark &&
               (FieldInfo[1] == "0" || FieldInfo[1] == "1" || FieldInfo[1] == "2")) continue;

            // For REMARK 0 & 1 citation parts, put them into JNRL tokens
            // instead of REMARK tokens
            if (token == "REMARK" && _input_format == NDB_FILE_FORMAT_PDB && (FieldInfo[1] == "0" || FieldInfo[1] == "1")) {
                 _processTextRecord(FieldInfo, false);
                 if (FieldInfo[1] == "0") {
                      if (first_remark_0) {
                           serial_no = 0;
                           start_remark_citation = false;
                      }
                      first_remark_0 = false;
                      prefix = ORIGINAL_CITATION;
                 } else if (FieldInfo[1] == "1") {
                      if (first_remark_1) {
                           serial_no = 0;
                           start_remark_citation = false;
                      }
                      first_remark_1 = false;
                      prefix = "";
                 }
                 bool is_REFERENCE_line = false; 
                 if (FieldInfo[2].find("REFERENCE") != std::string::npos) {
                      std::string cs = FieldInfo[2];
                      String::StripAndCompressWs(cs);
                      if (cs.substr(0, 9) == "REFERENCE") is_REFERENCE_line = true;
                 }
                 if (is_REFERENCE_line) {
                      start_remark_citation = true;
                      serial_no++;
                      FieldInfo.clear();
                 } else if (FieldInfo[2].find("PDB ID:") != std::string::npos) {
                      // _processTextRecord(FieldInfo);
                      FieldInfo.clear();
                      continue;
                 } else if (start_remark_citation) {
                      token = "JRNL";
                      _get_one_line_info(FieldInfo, token, line);
                      if (!FieldInfo.empty()) {
                           FieldInfo[1] = prefix + String::IntToString(serial_no);
                           _processTextRecord(FieldInfo);
                      }
                 } else {
                      // _processTextRecord(FieldInfo);
                      FieldInfo.clear();
                 }
            } else if (NdbToken::IsAtomToken(token)) { 
                 bool status = _processAtomRecord(FieldInfo, lineno);
                 if (!status && _messageIo && _format_checking) _write_aux_atom_record_error(token, line, lineno, last_line);
            } else _processTextRecord(FieldInfo);
       }
       fclose (infile);

       _removeLastSemicolon();
       _mergeConcatenateRecords();

       if (_molecules->size() == 1) (*_molecules)[0]->set_Mol_ID(1);

       if (_input_format == NDB_FILE_FORMAT_PDB) NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_NDB);

       if (!_messageIo || !_format_checking) return;

       if (found_model) _write_no_ENDMDL_error(model_record);

       if (_molecules->size() > 0 && (num_ter == 0 || not_found_TER_card_before_ENDMDL)) {
            std::string error = "All polymer chains should have a TER card at the end of each chain. ";
            error += "No 'TER' records after each polymer. Please add TER card and restart again.";
            _messageIo->insertMessage("error", "model", error, true);
       }

       if (!_wrong_PDB_chain_ID.empty()) {
            std::string chain_ids;
            chain_ids.clear();
            for (std::set<std::string>::const_iterator pos = _wrong_PDB_chain_ID.begin(); pos != _wrong_PDB_chain_ID.end(); ++pos) {
                 if (!chain_ids.empty()) chain_ids += ", ";
                 chain_ids += "'" + *pos + "'";
            }

            std::string error = "Each polymer must be assigned a unique chain ID and only alphanumeric chain ";
            error += "IDs (A-Z, 0-9, a-z) in column 22 are allowed. ";
            if (_wrong_PDB_chain_ID.size() > 1)
                 error += "The PDB chain IDs ";
            else error += "The PDB chain ID ";
            error += chain_ids;
            if (_wrong_PDB_chain_ID.size() > 1)
                 error += " are not allowed.";
            else error += " is not allowed.";
            _messageIo->insertMessage("error", "model", error, true);
       }

       if (num_end > 1 && end_in_middle) {
            std::string error = "File has " + String::IntToString(num_end) + " 'END' tokens and in the middle of file.";
            error += " There should only be one 'END' token at the end of the file. Please remove extra END token(s) and restart again.";
            _messageIo->insertMessage("error", "model", error, true);
       } else if (end_in_middle) {
            _messageIo->insertMessage("error", "model",
                "File has 'END' token in the middle of file. Please move 'END' token to the end of file and restart again.", true);
       }
}

void PdbRead::_get_token_from_line(std::string &token, const std::string &line)
{
       token.clear();
       std::string::size_type index = line.find_first_of(" \t\n", 0);
       if (index == std::string::npos) index = line.length();
       if (index > 6) index = 6;

       if (index == 6 && line.substr(0, 4) == "ATOM" && isdigit(line[4])) index = 4;

       token = line.substr(0, index);
       String::UpperCase(token);
}

bool PdbRead::_get_one_line_info(std::vector<std::string>& FieldInfo, std::string& token, const std::string& line)
{
       FieldInfo.clear();
       if (token == "JRNL") {
            if (line.size() < 17) return false;
            std::string jrnl_token = line.substr(12, 4);
            String::StripAndCompressWs(jrnl_token);
            if (!NdbToken::IsJrnlToken(jrnl_token)) return false;
            token = jrnl_token;
       }

       const ndb_token_format& ndbformat = NdbToken::getTokenFormat(token);
       if (token == "DBREF1") {
            const ndb_token_format& ndbformat1 = NdbToken::getTokenFormat("DBREF");
            for (int i = 0; i < ndbformat1.NumField; ++i) FieldInfo.push_back("");
       } else for (int i = 0; i < ndbformat.NumField; ++i) FieldInfo.push_back("");

       bool checking_flag = false;
       if (_check_token_set.find(token) != _check_token_set.end() && _messageIo && _format_checking) checking_flag = true;
       bool format_error = false;

       int lastColumn = 0;
       std::string not_a_number_error = "";
       int count = 0;
       for (int i = 0; i < ndbformat.NumField; ++i) {
            if (_input_format == NDB_FILE_FORMAT_PDB && ndbformat.FieldList[i].NdbOnly) break;
            int last = ndbformat.FieldList[i].FieldStcol + ndbformat.FieldList[i].FieldWidth - 1;
            if (last > lastColumn) lastColumn = last;
            _get_field_value(i, FieldInfo[i], token, ndbformat.FieldList[i], line);
            if (checking_flag) {
                 if (FieldInfo[i].empty() && ndbformat.FieldList[i].Mandatory) {
                      if ((token == "ATOM") || (token == "HETATM")) {
                           if (i < _atom_token_field) format_error = true;
                      } else format_error = true;
                 }
                 if (!FieldInfo[i].empty() && _atom_token_set.find(token) != _atom_token_set.end() &&
                     ndbformat.FieldList[i].FieldType < 3 && !String::IsNumber(FieldInfo[i])) {
                     if (!not_a_number_error.empty()) not_a_number_error += " ";
                     not_a_number_error += FieldInfo[i];
                     count++;
                 }
            }
       }
       if (!not_a_number_error.empty()) {
            std::string error = "'" + not_a_number_error + "'";
            if (count > 1)
                 error += " are not numbers.";
            else error += " is not a number.";
            _messageIo->insertMessage("error", "model", error, true, line, true);
       }

       if (checking_flag) {
            std::string tmp_line = line;
            if (_atom_token_set.find(token) != _atom_token_set.end()) {
                 if (FieldInfo[5].empty()) _no_pdb_chain_id_flag = true;
                 else if (!IsAlnum(FieldInfo[5])) _wrong_PDB_chain_ID.insert(FieldInfo[5]);
 
                 // remove segID -- Segment identifier 
                 int start = 72;
                 int end = 76;
                 if ((int) tmp_line.size() < end) end = tmp_line.size();
                 for (int i = start; i < end; ++i) tmp_line[i] = ' ';
            }

            std::string value;
            int last_column = tmp_line.size();
            if (lastColumn < last_column) last_column = lastColumn;
            for (int i = 1; i <= ndbformat.NumField; ++i) {
                 int start_column = ndbformat.FieldList[i-1].FieldStcol + ndbformat.FieldList[i-1].FieldWidth - 1;
                 if (start_column >= last_column) break;

                 int length = last_column - start_column;
                 if (i < ndbformat.NumField) {
                      int len = ndbformat.FieldList[i].FieldStcol - ndbformat.FieldList[i-1].FieldStcol - ndbformat.FieldList[i-1].FieldWidth;
                      if (len < length) length = len;
                 }

                 value.clear();
                 if (length > 0) {
                      value = tmp_line.substr(start_column, length);
                      String::StripAndCompressWs(value);
                 }
                 if (!value.empty()) format_error = true;
            }
       }

       if (token == "DBREF1") 
            FieldInfo[0] = "DBREF";
       else if (NdbToken::IsJrnlToken(token)) {
            int field_no = ndbformat.SeqField - 1;
            FieldInfo[0] = "JRNL";
            if (FieldInfo[field_no].empty()) FieldInfo[field_no] = "primary";
       }
       if (!FieldInfo.empty()) {
            String::UpperCase(FieldInfo[0]);
            std::string::size_type pos = FieldInfo[0].find(" ");
            if (pos != std::string::npos) FieldInfo[0].erase(pos);
       }

       return format_error;
}

void PdbRead::_get_field_value(const int field_no, std::string& FieldAnswer, const std::string& token,
                               const field_format& FieldList, const std::string& line)
{
       FieldAnswer.clear();
       int start_column = FieldList.FieldStcol - 1;
       int length = FieldList.FieldWidth;
       if (start_column >= (int) line.size()) return;

       FieldAnswer = line.substr(start_column, length);
       bool strip_leading_flag = true;
       if (((token == "REMARK" || token == "PREMRK" || token == "PCOMPN" ||
             token == "PSOURC") && field_no == 2) ||
            (token == "PDBRMK" && field_no == 3))
            strip_leading_flag = false;
       if (strip_leading_flag) String::StripLeadingWs(FieldAnswer);
       String::StripTrailingWs(FieldAnswer);
}

void PdbRead::_processTextRecord(std::vector<std::string>& FieldInfo, const bool& clear_flag)
{
       if (FieldInfo.empty()) return;

       if (FieldInfo[0] == "DBREF2") {
             std::map<std::string, std::list<std::vector<std::string> > >::iterator pos = _pdb_records->find("DBREF");
             if (pos != _pdb_records->end() && !pos->second.empty()) {
                  std::vector<std::string>& data = pos->second.back();
                  data[9] = data[8];
                  data[8] = FieldInfo[3]; 
                  for (int i = 0; i < 4; ++i) {
                       data[i + 10] = FieldInfo[i + 4];
                  }
             }
             if (clear_flag) FieldInfo.clear();
             return;
       }

       std::string token = FieldInfo[0];
       if (token == "JRNL") token = FieldInfo[2];
       const ndb_token_format& ndbformat = NdbToken::getTokenFormat(token);
       token = ndbformat.TokenName;

       bool is_continued_value = false;
       int field_no = ndbformat.ContinuedF - 1;
       if (token != "REMARK" && token != "PREMRK" && token != "PCOMPN" && token != "PSOURC" && field_no >= 0) {
           int field_no1 = ndbformat.ConField - 1;
           if (field_no1 >= 0 && !FieldInfo[field_no1].empty())
                is_continued_value = true;
           field_no1 = ndbformat.PdbConField - 1;
           if (field_no1 >= 0 && !FieldInfo[field_no1].empty())
                is_continued_value = true;
           if ((token == "COMPND" || token == "SOURCE") && _is_new_field(token, FieldInfo[field_no]))
                is_continued_value = false;
       }

       if (is_continued_value) {
             std::map<std::string, std::list<std::vector<std::string> > >::iterator pos = _pdb_records->find(token);
             if (pos != _pdb_records->end() && !pos->second.empty()) {
                  std::vector<std::string>& data = pos->second.back();
/*
                  if ((token == "COMPND" || token == "SOURCE") &&
                      !data[field_no].empty()) {
                       int size = data[field_no].size() - 1;
                       if (data[field_no][size] == ';') is_continued_value = false;
                  }
                  if (is_continued_value) {
*/
                       _concatenate_continuation_field(token, data[field_no], FieldInfo[field_no]);
                       if (clear_flag) FieldInfo.clear();
                       return;
                  // }
             }
       }

       std::map<std::string, std::list<std::vector<std::string> > >::iterator pos = _pdb_records->find(token);
       if (pos != _pdb_records->end()) pos->second.push_back(FieldInfo);
       else {
            std::list<std::vector<std::string> > t_list;
            t_list.clear();
            t_list.push_back(FieldInfo);
            _pdb_records->insert(std::make_pair(token, t_list));
       }

       if (clear_flag) FieldInfo.clear();
}

void PdbRead::_removeLastSemicolon()
{
       std::vector<std::string> tokens;
       tokens.clear();
       tokens.push_back("COMPND");
       tokens.push_back("SOURCE");
       for (std::vector<std::string>::iterator vpos = tokens.begin(); vpos != tokens.end(); ++vpos) {
             std::map<std::string, std::list<std::vector<std::string> > >::iterator
                 pos = _pdb_records->find(*vpos);
             if (pos == _pdb_records->end()) continue;
             if (pos->second.empty()) continue;

             const ndb_token_format& ndbformat = NdbToken::getTokenFormat(*vpos);
             int field_no = ndbformat.ContinuedF - 1;
             for (std::list<std::vector<std::string> >::iterator lpos = pos->second.begin(); lpos != pos->second.end(); ++lpos) {
                  if (!(*lpos)[field_no].empty()) {
                       int size = (*lpos)[field_no].size() - 1;
                       if ((*lpos)[field_no][size] == ';') (*lpos)[field_no].erase(size);
                  }
             }
       }
}

void PdbRead::_mergeConcatenateRecords()
{
       std::vector<std::string> tokens;
       tokens.clear();
       tokens.push_back("TITLE");
       tokens.push_back("KEYWDS");
       tokens.push_back("AUTHOR");

       std::string text;
       std::list<std::vector<std::string> > t_list;

       for (std::vector<std::string>::iterator vpos = tokens.begin(); vpos != tokens.end(); ++vpos) {
             std::map<std::string, std::list<std::vector<std::string> > >::iterator
                 pos = _pdb_records->find(*vpos);
             if (pos == _pdb_records->end() || pos->second.size() < 2) continue;

             text.clear();
             for (std::list<std::vector<std::string> >::iterator lpos = pos->second.begin(); lpos != pos->second.end(); ++lpos) {
                  if ((*lpos)[2].empty()) continue;
                  if (!text.empty()) text += " ";
                  text += (*lpos)[2];
             }

             std::vector<std::string> t_vec = pos->second.front();
             t_vec[2] = text;

             t_list.clear();
             t_list.push_back(t_vec);
             pos->second.clear();
             pos->second = t_list;
       }
}

bool PdbRead::_is_new_field(const std::string& token, const std::string &value)
{
       if (token == "COMPND") {
           for (int i = 0; i < NUM_COMPND; ++i) {
                if (value.find(__compnd_token[i]) != std::string::npos) return true;
           }
       } else if (token == "SOURCE") {
           for (int i = 0; i < NUM_SOURCE; ++i) {
                if (value.find(__source_token[i]) != std::string::npos) return true;
           }
       }
       return false;
}

void PdbRead::_concatenate_continuation_field(const std::string& token, std::string& FieldInfo, const std::string& NewFieldInfo)
{
       if (NewFieldInfo.empty()) return;
       if (FieldInfo.empty()) FieldInfo = NewFieldInfo;
       else {
            int size = FieldInfo.size() - 1;
            if (token == "PDBRMK")
                 FieldInfo += "\n" + NewFieldInfo;
            if (FieldInfo[size] == '-')
                 FieldInfo += NewFieldInfo;
            else FieldInfo += " " + NewFieldInfo;
       }
}

bool PdbRead::_processAtomRecord(std::vector<std::string>& FieldInfo, const int& lineno)
{
       if (FieldInfo[0] == "ENDMDL") {
            FieldInfo.clear();
            return true;
       }

       if (FieldInfo[0] == "MODEL" || _first_atom) {
            _first_atom = false;
            RCSB::Molecule* mol = new RCSB::Molecule; 
            if (FieldInfo[0] == "MODEL") {
                 int mol_id = atoi(FieldInfo[1].c_str());
                 mol->set_Mol_ID(mol_id);
                 if (_model_set.find(mol_id) != _model_set.end()) {
                      std::string error = "PdbRead::_processAtomRecord: duplicate MODEL number " + String::IntToString(mol_id) + "\n";
                      _logIo->message(error.c_str());
                      if (_messageIo && _format_checking) {
                           std::string error = "Duplicate MODEL number " + String::IntToString(mol_id)
                                             + ". The model number must be between columns 11 and 14. Please make sure ";
                           error += "the MODEL numbers are in correct column and are in sequential order.";
                           _messageIo->insertMessage("error", "model", error, true);
                      }
                 }
                 _model_set.insert(mol_id);
            } else mol->set_Mol_ID(1);
            mol->set_index(_molecules->size());
            _molecules->push_back(mol);
            if (FieldInfo[0] == "MODEL") {
                 FieldInfo.clear();
                 return true;
            }
       }

       unsigned int last = _molecules->size() - 1;
       bool ok = (*_molecules)[last]->insert_a_atom(FieldInfo, true, lineno);

       FieldInfo.clear();

       return ok;
}

void PdbRead::_check_line_not_start_from_first_column(const std::string& line, const int& lineno)
{
       std::string tmp_line = line;
       String::StripAndCompressWs(tmp_line);
       if (tmp_line.empty()) return;

       std::string error = "Line " + String::IntToString(lineno) + " does not start in first column.";
       error += " All tokens must start in the first column. Please modify the file so that the record in line "
              + String::IntToString(lineno) + " start in column 1.";

       _messageIo->insertMessage("error", "model", error, true, line, true);
}

void PdbRead::_write_no_ENDMDL_error(const std::string& model_record)
{
       std::string error = "No 'ENDMDL' found for '" + model_record + "'. Please add the ";
       error += "'ENDMDL' card at the end of the coordinates for each model.";

       _messageIo->insertMessage("error", "model", error, true);
}

void PdbRead::_write_aux_atom_record_error(const std::string& token, const std::string& line, const int& lineno, const std::string& last_line)
{
       std::map<std::string, int>::iterator pos = _error_count.find(token);
       if (pos == _error_count.end()) _error_count.insert(std::make_pair(token, 1));
       else {
            pos->second++;
            // if (pos->second > MAXIMUM_ERROR) return;
       }

       std::string error = "In line " + String::IntToString(lineno) + ", record '" + token
                         + "' does not have corresponding primary ATOM/HETATM coordinate record.\n\n" + line + "\n\nLast coordinate record:\n"
                         + last_line; 

       _messageIo->insertMessage("error", "model", error);
}

void PdbRead::_write_record_error(const std::string& token, const std::string& line, const int& lineno)
{
       if (_atom_token_set.find(token) != _atom_token_set.end()) {
            std::map<std::string, int>::iterator pos = _error_count.find(token);
            if (pos == _error_count.end()) _error_count.insert(std::make_pair(token, 1));
            else {
                 pos->second++;
                 // if (pos->second > MAXIMUM_ERROR) return;
            }
       }

       std::string error = "In line " + String::IntToString(lineno) + ", the format of " +  token
                         + " record is not in compliance with PDB V3.3 format. Please refer to the ";
       error += "<a href=\"http://www.wwpdb.org/documentation/format33/v3.3.html\"> PDB V3.3 format</a> and make sure that all values in line "
              + String::IntToString(lineno) + " are in proper columns.";
       _messageIo->insertMessage("error", "model", error, true, line, true);
}
