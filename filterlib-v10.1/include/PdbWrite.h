/*
FILE:     PdbWrite.h
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
#ifndef _H_PDB_WRITE_H_
#define _H_PDB_WRITE_H_

#include <list>
#include <map>
#include <vector>
#include <set>
#include <string>

#include "ConnectDic.h"
#include "LogUtil.h"
#include "Molecule.h"
#include "NdbToken.h"

class PdbWrite
{
   private:
       LogUtil *_logIo;
       ConnectDic *_ccDic;
       std::map<std::string, std::list<std::vector<std::string> > >* _pdb_records;
       std::vector<RCSB::Molecule*>* _molecules;
       std::vector<std::string> _original_entry_ids;
       std::list<std::pair<RCSB::Atom*, RCSB::Atom*> > _atom_connects;
       std::string _filename;
       bool _segid_flag;
       bool _extended_flag;
       int _output_format;
       int _num_remark;
       int _num_atom;
       int _num_ter;
       std::map<std::string, std::map<int, int> > _atomField_mapping;
       std::map<std::string, std::pair<int, int> > _field_length_mapping;
       std::set<std::string> _affected_field_set;

       void _init();
       void _write();
       void _write(FILE* fp, const int& header_section_only = false);
       void _writeGeneralRecord(std::vector<std::string>& lines, const ndb_token_format& ndbformat, std::vector<std::string>& record,
                               std::string& prevSeqNo, int& PdbContNo, const bool& full_record_flag = true);
       void _write_ter_card(std::string& text, std::vector<std::string>& FieldInfo);
       void _write_ter_card(FILE* fp, std::vector<std::string>& FieldInfo);
       void _writeJrnlRecords(std::string& text);
       void _writeJrnlRecords(FILE* fp);
       void _updateUpperCase(const ndb_token_format& ndbformat, std::vector<std::string>& record);
       void _getJrnlContinuedField(const std::string& TokenName, const int max_len, const std::string& value, std::vector<std::string> &field);
       void _getContinuedField(const int default_max_len, const int PdbContNo, const std::string& value, std::vector<std::string> &field);
       std::string _printSpaceBetweenField(const ndb_token_format& ndbformat, const int field_no);
       std::string _printFieldValue(const ndb_token_format& ndbformat, const int& field_no, const std::vector<std::string>& Fields,
                                    const std::string& continuedField, const int& ContNo, const int& PdbContNo);
       std::string _FormattedFieldValue(const std::string& value, const field_format& FieldList);
       void _setEightyCharacters(std::string& line);
       void _getJrnlIndexSet(std::set<std::string>& index_set);
       void _getSecondCitationIndex(std::vector<std::string>& index_array, std::set<std::string>& index_set, const std::string& prefix);
       void _writeSelectedJrnlRecords(std::vector<std::string>& lines, std::vector<std::string>& index_array);
       void _writeSelectedJrnlRecords(FILE* fp, std::vector<std::string>& index_array);
       void _writeSelectedJrnlRemarkRecords(const std::string& remark_no, std::vector<std::string>& index_array);
       int _getResField(const std::string& TokenName, const int& field_no);
       void _getConnectRecord();
       void _insert_connect(std::map<int, std::set<int> >& connect, const int& atom_id1, const int& atom_id2);
       void _getMasterRecord();
   public:
       PdbWrite();
       ~PdbWrite();
       void setLog(LogUtil *logPt);
       void setCCDic(ConnectDic *ccdic);
       void setRecord(std::map<std::string, std::list<std::vector<std::string> > >* records);
       void setMolecule(std::vector<RCSB::Molecule*>* mols);
       void setSEGIDFlag(const bool& segid_flag);
       void setOriginalEntryIds(const std::vector<std::string>& ids);
       void setConnect(const std::list<std::pair<RCSB::Atom*, RCSB::Atom*> >& connects);
       void setFileName(const std::string& filename);
       void setOutputFormat(const int& format);
       void setExtendedFlag();
       void setFieldLabel(const std::string& fieldLabel);
       void WriteRecord();
       void WriteRecord(FILE* fp);
       void WriteRecord(std::string& text);
       void WriteRecord(const std::string& filename, const int format, std::map<std::string, std::list<std::vector<std::string> > >* records,
                        std::vector<RCSB::Molecule*>* mols, const std::list<std::pair<RCSB::Atom*, RCSB::Atom*> >& connects);
       void WriteRecord(const std::string& filename, const std::vector<RCSB::Residue*>& residues);
       void WriteRecord(const std::string& filename, const std::vector<RCSB::Atom*>& atoms, const std::list<std::pair<RCSB::Atom*, RCSB::Atom*> >& connects = 
                        std::list<std::pair<RCSB::Atom*, RCSB::Atom*> >());
       void WriteRecord(std::string& text, const std::list<std::pair<std::string, RCSB::Atom*> >& atom_list);
       void WriteRecord(const std::string& filename, const std::list<std::pair<std::string, RCSB::Atom*> >& atom_list,
                        const std::list<std::pair<RCSB::Atom*, RCSB::Atom*> >& connects);
       void writeGeneralRecords(FILE *fp, std::map<std::string, std::list<std::vector<std::string> > >* records);
       void writeGeneralRecords(std::vector<std::string>& lines, const ndb_token_format& ndbformat, std::list<std::vector<std::string> >& records,
                                const bool& full_record_flag = true);
       void writeGeneralRecords(FILE* fp, const ndb_token_format& ndbformat, std::list<std::vector<std::string> >& records,
                                const bool& full_record_flag = true);
       void writeGeneralRecord(std::vector<std::string>& lines, const ndb_token_format& ndbformat, std::vector<std::string>& record,
                               const bool& full_record_flag = true);
       void writeGeneralRecord(FILE* fp, const ndb_token_format& ndbformat, std::vector<std::string>& record, const bool& full_record_flag = true);
       void WriteCoordinates(std::string& text);
       void WriteCoordinates(FILE *fp);
       void WriteMolecule(std::string& text, RCSB::Molecule* molecule);
       void WriteMolecule(FILE *fp, RCSB::Molecule* molecule);
       void WriteAtom(std::vector<std::string>& records, const std::string& chain_type, RCSB::Atom* atom, std::vector<std::string>& FieldInfo,
                      const bool include_auxiliary_atoms = true);
       void WriteAtom(FILE* fp, const std::string& chain_type, RCSB::Atom* atom, std::vector<std::string>& FieldInfo, const bool include_auxiliary_atoms = true);
       void ResetAtomCounting();
       void writeMasterRecord(FILE* fp);
       void getAffectedFieldSet(const ndb_token_format& ndbformat, const std::list<std::vector<std::string> >& records, std::set<std::string>& fieldSet);
       void updateNdbTokenFormat(const std::set<std::string>& fieldSet, ndb_token_format& ndbformat); 
};

extern std::string printAtomNameField(ConnectDic *ccdic, const std::string& atomtype, const std::string& atomName, const std::string& residueName);

#endif
