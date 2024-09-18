/*
FILE:     cif-util.C
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
#include <string.h>
#include <sys/stat.h>
#include <math.h>

#include <exception>
#include <string>

#include "utillib.h"
#include "CifFileReadDef.h"
#include "CifParserBase.h"

static std::string __getIDs_from_database_2(CifFile *fobj, const std::string& database_id, const bool& is_upper_case = true);
static std::string __getIDs_from_database_2(Block &block, const std::string& database_id, const bool& is_upper_case = true);

void find_rcsb_directory(const std::string &id, std::string& dirname)
{
       struct stat stat1;

       dirname = "/annotation/prot/" + id;
       if (stat(dirname.c_str(), &stat1) == 0) return;

       dirname = "/annotation/nmr/" + id;
       if (stat(dirname.c_str(), &stat1) == 0) return;

       dirname = "/annotation/ndb/" + id;
       if (stat(dirname.c_str(), &stat1) == 0) return;

       dirname.clear();
}

int file_corres_file(const std::string &id, std::string& filename)
{
       struct stat stat1;

       filename = "/annotation/prot/" + id + "/" + id + ".cor.cif";
       if (stat(filename.c_str(), &stat1) == 0) return 1;

       filename = "/annotation/nmr/" + id + "/" + id + ".cor.cif";
       if (stat(filename.c_str(), &stat1) == 0) return 1;

       filename = "/annotation/ndb/" + id + "/" + id + ".cor.cif";
       if (stat(filename.c_str(), &stat1) == 0) return 1;

       filename.clear();
       return 0;
}

int find_rcsb_file(const std::string &id, std::string& filename)
{
       struct stat stat1;

       filename = "/annotation/prot/" + id + "/" + id + ".cif";
       if (stat(filename.c_str(), &stat1) == 0) return 1;

       filename = "/annotation/nmr/" + id + "/" + id + ".cif";
       if (stat(filename.c_str(), &stat1) == 0) return 1;

       filename = "/annotation/ndb/" + id + "/" + id + ".cif";
       if (stat(filename.c_str(), &stat1) == 0) return 1;

       if (id.size() == 4) {
            std::string hash = id.substr(1, 2);
            filename = "/annotation/new-legacy/" + hash + "/"
                     + id + "/" + id + ".cif";
            if (stat(filename.c_str(), &stat1) == 0) return 1;
       }

       filename.clear();
       return 0;
}

int find_rcsb_file(const std::string &id, std::string& filename, std::string &expfile, std::string& csfile)
{
       struct stat stat1;

       filename = "/annotation/prot/" + id + "/" + id + ".cif";
       expfile = "/annotation/prot/" + id + "/" + id +  "-sf.cif";
       if (stat(filename.c_str(), &stat1) == 0) {
            if (stat(expfile.c_str(), &stat1) != 0) expfile.clear();
            return 1;
       }

       filename = "/annotation/nmr/" + id + "/" + id + ".cif";
       expfile  = "/annotation/nmr/" + id + "/" + id + ".mr";
       csfile   = "/annotation/nmr/" + id + "/" + id + CS_EXTENSION;
       if (stat(filename.c_str(), &stat1) == 0) {
            if (stat(expfile.c_str(), &stat1) != 0) expfile.clear();
            if (stat(csfile.c_str(), &stat1) != 0) csfile.clear();
            return 1;
       }

       filename = "/annotation/ndb/" + id + "/" + id + ".cif";
       expfile =  "/annotation/ndb/" + id + "/" + id + "-sf.cif";
       if (stat(filename.c_str(), &stat1) == 0) {
            if (stat(expfile.c_str(), &stat1) != 0) expfile.clear();
            return 1;
       }

       if (id.size() == 4) {
            std::string hash = id.substr(1, 2);
            filename = "/annotation/new-legacy/" + hash + "/" + id + "/" + id + ".cif";
            if (stat(filename.c_str(), &stat1) == 0) {
                 expfile = "/annotation/new-legacy/" + hash + "/" + id + "/" + id + "-sf.cif";
                 if (stat(expfile.c_str(), &stat1) != 0) {
                      expfile = "/annotation/new-legacy/" + hash + "/" + id + "/" + id + ".mr";
                      if (stat(expfile.c_str(), &stat1) != 0) expfile.clear();
                 }
                 return 1;
            }
       }

       filename.clear();
       expfile.clear();

       return 0;
}

void get_unique_values(ISTable *Table, const std::string& colName, std::set<std::string>& value_set)
{
       value_set.clear();
       if (!Table) return;

       std::string cs; 
       for (unsigned int i = 0; i < Table->GetNumRows(); ++i) {
            get_value_clean(cs, Table, i, colName);
            if (!cs.empty()) value_set.insert(cs);
       }
}

int get_value(std::string &cifstring, ISTable *Table, const unsigned int& row, const std::string& colName)
{
       cifstring.clear();
       try {
            cifstring = (*Table)(row, colName);
            if (cifstring == "." || cifstring == "?" /* || cifstring == "@" */ ) {
                 cifstring.clear();
                 return 0;
            }
            return 1;
       } catch (std::exception &exc) {
            return -1;
       }
}

int get_value_clean(std::string &cifstring, ISTable *Table, const unsigned int& row, const std::string& colName)
{
       int ret = get_value(cifstring, Table, row, colName);
       if (ret == 1) String::StripAndCompressWs(cifstring);
       return ret;
}

int get_value_clean_lower(std::string &cifstring, ISTable *Table, const unsigned int& row, const std::string& colName)
{
       std::string tmpstring;
       cifstring.clear();
       int ret = get_value(tmpstring, Table, row, colName);
       if (ret == 1) {
            String::LowerCase(tmpstring, cifstring);
            String::StripAndCompressWs(cifstring);
       }
       return ret;
}

int get_value_clean_upper(std::string &cifstring, ISTable *Table, const unsigned int& row, const std::string& colName)
{
       std::string tmpstring;
       cifstring.clear();
       int ret = get_value(tmpstring, Table, row, colName);
       if (ret == 1) {
            String::UpperCase(tmpstring, cifstring);
            String::StripAndCompressWs(cifstring);
       }
       return ret;
}

int get_value_lower(std::string &cifstring, ISTable *Table, const unsigned int& row, const std::string& colName)
{
       std::string tmpstring;
       cifstring.clear();
       int ret = get_value(tmpstring, Table, row, colName);
       if (ret == 1) String::LowerCase(tmpstring, cifstring);
       return ret;
}

int get_value_upper(std::string &cifstring, ISTable *Table, const unsigned int& row, const std::string& colName)
{
       std::string tmpstring;
       cifstring.clear();
       int ret = get_value(tmpstring, Table, row, colName);
       if (ret == 1) String::UpperCase(tmpstring, cifstring);
       return ret;
}

void get_values(ISTable *Table, std::vector<std::vector<std::string> >& values)
{
       values.clear();
       if (!Table) return;

       const std::vector<std::string>& itemNames = Table->GetColumnNames();

       get_values(Table, itemNames, values);
}

void get_values(ISTable *Table, std::vector<std::map<std::string, std::string> >& values)
{
       values.clear();
       if (!Table) return;

       const std::vector<std::string>& itemNames = Table->GetColumnNames();

       std::map<std::string, std::string> mapping;
       std::string cs; 

       for (unsigned int i = 0; i < Table->GetNumRows(); ++i) {
            mapping.clear();
            for (std::vector<std::string>::const_iterator pos = itemNames.begin(); pos != itemNames.end(); ++pos) {
                 get_value(cs, Table, i, *pos);
                 if (!cs.empty()) mapping.insert(std::make_pair(*pos, cs));
            }
            if (!mapping.empty()) values.push_back(mapping);
       }
}

void get_values(ISTable *Table, std::list<std::map<std::string, std::string> >& values)
{
       values.clear();
       if (!Table) return;

       const std::vector<std::string>& itemNames = Table->GetColumnNames();

       std::map<std::string, std::string> mapping;
       std::string cs; 

       for (unsigned int i = 0; i < Table->GetNumRows(); ++i) {
            mapping.clear();
            for (std::vector<std::string>::const_iterator pos = itemNames.begin(); pos != itemNames.end(); ++pos) {
                 get_value(cs, Table, i, *pos);
                 if (!cs.empty()) mapping.insert(std::make_pair(*pos, cs));
            }
            if (!mapping.empty()) values.push_back(mapping);
       }
}

void get_values(ISTable *Table, const std::vector<std::string>& itemNames, std::vector<std::vector<std::string> >& values)
{
       values.clear();
       if (!Table) return;

       std::vector<std::string> data;
       std::string cs; 

       for (unsigned int i = 0; i < Table->GetNumRows(); ++i) {
            data.clear();
            for (std::vector<std::string>::const_iterator pos = itemNames.begin(); pos != itemNames.end(); ++pos) {
                 get_value(cs, Table, i, *pos);
                 data.push_back(cs);
            }
            values.push_back(data);
       }
}

void get_values_clean(ISTable *Table, const std::vector<std::string>& itemNames, std::vector<std::vector<std::string> >& values)
{
       values.clear();
       if (!Table) return;

       std::vector<std::string> data;
       std::string cs; 

       for (unsigned int i = 0; i < Table->GetNumRows(); ++i) {
            data.clear();
            for (std::vector<std::string>::const_iterator pos = itemNames.begin(); pos != itemNames.end(); ++pos) {
                 get_value_clean(cs, Table, i, *pos);
                 data.push_back(cs);
            }
            values.push_back(data);
       }
}

void get_conditional_value(Block& block, std::string& value, const std::string& category, const std::string& item,
                           const std::string& condition_item, const std::string& condition_value)
{
       value.clear();

       ISTable *t = getTablePtr(block, category);
       if (!t) return;

       std::string cs;
       for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
            get_value_clean(cs, t, i, condition_item);
            if (cs == condition_value) {
                 get_value_clean(value, t, i, item);
                 return;
            }
       }
}

void update_table(ISTable* t, const std::vector<std::map<std::string, std::string> >& data)
{
       if (data.empty()) return;

       for (std::vector<std::map<std::string, std::string> >::const_iterator vpos = data.begin(); vpos != data.end(); ++vpos) {
            unsigned int row = t->GetNumRows();
            t->AddRow();
            for (std::map<std::string, std::string>::const_iterator mpos = vpos->begin(); mpos != vpos->end(); ++mpos) {
                 t->UpdateCell(row, mpos->first, mpos->second);
            }
       }
}

void update_table(ISTable* t, const std::vector<std::string>& items, const std::vector<std::vector<std::string> >& values)
{
       if (values.empty()) return;

       for (std::vector<std::vector<std::string> >::const_iterator vpos = values.begin(); vpos != values.end(); ++vpos) {
            unsigned int row = t->GetNumRows();
            t->AddRow();
            for (unsigned int i = 0; i < items.size(); ++i) {
                 t->UpdateCell(row, items[i], (*vpos)[i]);
            }
       }

}

bool update_value(ISTable* t, const std::string& item, const unsigned int& row, const std::string& value, const bool& force_update)
{
       if (!t || value.empty()) return false;

       if (!t->IsColumnPresent(item) || (row >= t->GetNumRows())) return false;

       std::string cs;
       get_value_clean(cs, t, row, item);
       if ((cs == value) || (cs.empty() && (value == "?"))) return false;

       if (!cs.empty() && !force_update) return false;

       t->UpdateCell(row, item, value);
       return true;
}

int is_empty_row(ISTable *Table, const std::vector<std::string>& ColumnNames, const unsigned int& rowindex, const std::string& except)
{
       std::string cifstring;
       for (unsigned int i = 0; i < ColumnNames.size(); i++) {
            get_value(cifstring, Table, rowindex, ColumnNames[i]);
            if (!cifstring.empty()) {
                 if (ColumnNames[i] == except) continue;
                 return 0;
            }
       }
       return 1;
}

int is_empty_row(ISTable *Table, const std::vector<std::string>& ColumnNames, const unsigned int& rowindex, const std::set<std::string>& except)
{
       std::string cifstring;
       for (unsigned int i = 0; i < ColumnNames.size(); i++) {
            get_value(cifstring, Table, rowindex, ColumnNames[i]);
            if (!cifstring.empty()) {
                 if (except.find(ColumnNames[i]) != except.end()) continue;
                 return 0;
            }
       }
       return 1;
}

int is_empty_row_loose(ISTable *Table, const std::vector<std::string>& ColumnNames, const unsigned int& rowindex)
{
       std::string cifstring;

       for (unsigned int i = 0; i < ColumnNames.size(); i++) {
            get_value(cifstring, Table, rowindex, ColumnNames[i]);
            if (!cifstring.empty()) {
                 if (ColumnNames[i] == "entry_id" || ColumnNames[i] == "diffrn_id" || ColumnNames[i] == "entity_id" || ColumnNames[i] == "id") continue;
                 return 0;
            }
       }
       return 1;
}

int is_empty_table(ISTable *Table)
{
       if (Table == NULL) return 1;

       if (Table->GetNumRows() == 0) return 1;

       std::set<std::string> except_items;
       except_items.clear();
       except_items.insert("entry_id");
       except_items.insert("ordinal");
       except_items.insert("pdbx_ordinal");

       std::vector<std::string> ColumnNames = Table->GetColumnNames();
       for (unsigned int i = 0; i < Table->GetNumRows(); i++) {
            // if (!is_empty_row(Table, ColumnNames, i, "entry_id")) break;
            if (!is_empty_row(Table, ColumnNames, i, except_items)) return 0;
       }
       return 1;
}

int is_empty_table(ISTable *Table, const std::set<std::string>& except)
{
       if (Table == NULL) return 1;


       if (Table->GetNumRows() == 0) return 1;

       std::vector<std::string> ColumnNames = Table->GetColumnNames();
       for (unsigned int i = 0; i < Table->GetNumRows(); i++) {
            if (!is_empty_row(Table, ColumnNames, i, except)) return 0;
       }
       return 1;
}

int is_empty_table_loose(ISTable *Table)
{
       if (Table == NULL) return 1;

       if (Table->GetNumRows() == 0) return 1;

       std::vector<std::string> ColumnNames = Table->GetColumnNames();
       for (unsigned int i = 0; i < Table->GetNumRows(); i++) {
            if (!is_empty_row_loose(Table, ColumnNames, i)) return 0;
       }
       return 1;
}

int is_empty_table(ISTable *Table, const char *except)
{
       if (Table == NULL) return 1;

       if (Table->GetNumRows() == 0) return 1;

       std::vector<std::string> ColumnNames = Table->GetColumnNames();
       for (unsigned int i = 0; i < Table->GetNumRows(); i++) {
            if (!is_empty_row(Table, ColumnNames, i, except)) return 0;
       }
       return 1;
}

CifFile* read_cif_file(const std::string &odbfile, const std::string &textfile, std::string &diags, const bool& caseSensitive, const int& length)
{
       struct stat statbuf;
       if (stat(textfile.c_str(), &statbuf) != 0) {
            diags = "Can't find " + textfile + " file.";
            return NULL;
       }

       Char::eCompareType caseSense = Char::eCASE_INSENSITIVE;
       if (caseSensitive) caseSense = Char::eCASE_SENSITIVE;
       diags.clear();
       CifFileReadDef readDef;
       CifFile *fobj = NULL;
       if (!odbfile.empty())
            fobj = new CifFile(CREATE_MODE, odbfile, false, caseSense, length, CifString::UnknownValue);
       else fobj = new CifFile(false, caseSense, length, CifString::UnknownValue);
       CifParser *CifParserR = new CifParser(fobj, readDef, false);
       CifParserR->LogFileOff();
       CifParserR->Parse(textfile, diags);
       delete CifParserR;
       if (!diags.empty()) {
            delete fobj;
            return NULL;
       }
       return fobj;
}

CifFile *get_fobj(std::string& error, const std::string& ciffile, const bool& caseSensitive, const int& length)
{
       CifFile *fobj = read_cif_file("", ciffile, error, caseSensitive, length);
       if (!error.empty()) {
            if (fobj) delete fobj;
            return NULL;
       }
       return fobj;
}

CifFile *get_fobj(FILE *log, const std::string& odbfile, const std::string& ciffile, const bool& caseSensitive, const int& length)
{
       std::string cs;
       CifFile *fobj = read_cif_file(odbfile, ciffile, cs, caseSensitive, length);
       if (!cs.empty()) {
            fprintf(log, "%s: read cif file error:\n%s\n", ciffile.c_str(), cs.c_str());
            if (fobj) delete fobj;
            return NULL;
       }
       return fobj;
}

CifFile *get_fobj(LogUtil *log, const std::string& odbfile, const std::string& ciffile, const bool& caseSensitive, const int& length)
{
       std::string cs;
       CifFile *fobj = read_cif_file(odbfile, ciffile, cs, caseSensitive, length);
       if (!cs.empty()) {
            std::string error = ciffile + ": Read cif file failed\n" + cs + "\n";
            log->message(error.c_str());
            if (fobj) delete fobj;
            return NULL;
       }
       return fobj;
}

CifFile *get_fobj(FILE *log, const std::string& cifString, const bool& caseSensitive, const int& length)
{
       std::string diags;
       diags.clear();

       Char::eCompareType caseSense = Char::eCASE_INSENSITIVE;
       if (caseSensitive) caseSense = Char::eCASE_SENSITIVE;
       CifFile *fobj = new CifFile(false, caseSense, length, CifString::UnknownValue);

       CifParser *CifParserR = new CifParser(fobj, false);
       CifParserR->LogFileOff();
       CifParserR->ParseString(cifString, diags);
       delete CifParserR;

       if (!diags.empty()) {
            fprintf(log, "Read Cif String error:\n%s\n", diags.c_str());
            delete fobj;
            return NULL;
       }
       return fobj;
}

CifFile *create_fobj(const std::string& odbfile, const std::string& blockId, const bool& caseSensitive, const int& length)
{
       Char::eCompareType caseSense = Char::eCASE_INSENSITIVE;
       if (caseSensitive) caseSense = Char::eCASE_SENSITIVE;
       CifFile *fobj = NULL;
       if (!odbfile.empty())
            fobj = new CifFile(CREATE_MODE, odbfile, false, caseSense, length, CifString::UnknownValue);
       else fobj = new CifFile(false, caseSense, length, CifString::UnknownValue);
       fobj->AddBlock(blockId);
       return fobj;
}

void write_cif_file(const std::string& blockID, const std::string& filename, ISTable* table)
{
       std::string BlockID = blockID;
       if (BlockID.empty()) BlockID = "XXXX";
       CifFile* fobj = create_fobj("", BlockID);
       Block& block = fobj->GetBlock(BlockID);
       block.WriteTable(table);
       fobj->SetQuoting(CifFile::eDOUBLE);
       fobj->Write(filename);
       delete fobj;
}

void write_cif_file(const std::string& blockID, const std::string& filename, const std::vector<ISTable*>& tables)
{
       std::string BlockID = blockID;
       if (BlockID.empty()) BlockID = "XXXX";
       CifFile* fobj = create_fobj("", BlockID);
       Block& block = fobj->GetBlock(BlockID);
       for (std::vector<ISTable*>::const_iterator pos = tables.begin(); pos != tables.end(); ++pos) {
            block.WriteTable(*pos);
       }
       fobj->SetQuoting(CifFile::eDOUBLE);
       fobj->Write(filename);
       delete fobj;
}

void write_cif_file(const std::string& filename, const std::list<std::pair<std::string, std::vector<ISTable*> > >& bt_lists)
{
       CifFile *fobj = new CifFile(false, Char::eCASE_SENSITIVE, 250, CifString::UnknownValue);
       for (std::list<std::pair<std::string, std::vector<ISTable*> > >::const_iterator lpos = bt_lists.begin(); lpos != bt_lists.end(); ++lpos) {
            fobj->AddBlock(lpos->first);
            Block& block = fobj->GetBlock(lpos->first);
            for (std::vector<ISTable*>::const_iterator pos = lpos->second.begin(); pos != lpos->second.end(); ++pos) {
                 block.WriteTable(*pos);
            }
       }
       fobj->SetQuoting(CifFile::eDOUBLE);
       fobj->Write(filename);
       delete fobj;
}

std::string getSingleBestValue(Block& block, const std::string& category, const std::string& item, const int& number_compare_flag)
{
       // number_compare_flag = 0 for smaller, number_compare_flag = 1 for larger
       std::string value = "";
       ISTable* t = getTablePtr(block, category);
       if (!t) return value;

       std::string cs;
       for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
            get_value_clean(cs, t, i, item);
            if (cs.empty() || !String::IsNumber(cs)) continue;
            if (value.empty()) value = cs;
            else if (number_compare_flag) { // for larger
                 if (atof(cs.c_str()) > atof(value.c_str())) value = cs;
            } else { // for small 
                 if (atof(cs.c_str()) < atof(value.c_str())) value = cs;
            }
       }
       return value;
}

ISTable *getTableCopy(Block &block, const std::string& table_name)
{
       if (!block.IsTablePresent(table_name)) return NULL;

       ISTable *t = block.GetTablePtr(table_name);
       if (is_empty_table(t)) return NULL;

       ISTable *t1 = new ISTable(*t);
       return t1;
}

ISTable *getTablePtr(Block &block, const std::string& table_name)
{
       if (!block.IsTablePresent(table_name)) return NULL;

       ISTable *t = block.GetTablePtr(table_name);
       if (is_empty_table(t)) return NULL;
       return t;
}

void deleteTable(Block &block, const std::string& table_name)
{
       if (block.IsTablePresent(table_name)) block.DeleteTable(table_name);
}

ISTable *add_new_table(const std::string &tablename, const int& num_items, const char** items)
{
       ISTable *table = new ISTable(tablename);
       for (int i = 0; i < num_items; i++) {
            table->AddColumn(items[i]);
       }
       return table;
}

ISTable *add_new_table(const std::string &tablename, const std::vector<std::string> &items)
{
       ISTable *table = new ISTable(tablename);
       for (unsigned int i = 0; i < items.size(); i++) {
            table->AddColumn(items[i]);
       }
       return table;
}

ISTable *add_new_table(const std::string &tablename, const std::vector<std::string> &items,
                       const std::vector<std::vector<std::string> >& values)
{
       ISTable *t = add_new_table(tablename, items);
       int row = 0;
       for (std::vector<std::vector<std::string> >::const_iterator pos = values.begin(); pos != values.end(); ++pos) {
            t->AddRow();
            for (unsigned int i = 0; i < items.size(); ++i) {
                 t->UpdateCell(row, items[i], (*pos)[i]);
            }
            row++;
       }
       return t;
}

void check_missing_item(ISTable *table, const int& num_items, const char** items, const bool mandatary)
{
       std::vector<std::string> data;

       data.clear();
       if (mandatary)
            for (unsigned int i = 0; i < table->GetNumRows(); i++) data.push_back(".");
       else for (unsigned int i = 0; i < table->GetNumRows(); i++) data.push_back("");
       for (int i = 0; i < num_items; i++) {
            if (!table->IsColumnPresent(items[i])) {
                 table->AddColumn(items[i], data);
            }
       }
}

void check_missing_item(ISTable *table, const std::vector<std::string> &items)
{
       std::vector<std::string> data;

       data.clear();
       for (unsigned int i = 0; i < table->GetNumRows(); i++) data.push_back("");
       for (unsigned int i = 0; i < items.size(); i++) {
            if (!table->IsColumnPresent(items[i])) {
                 table->AddColumn(items[i], data);
            }
       }
}

void copy_item(ISTable *table, const std::string& existing_item, const std::string& new_item)
{
       if (table->IsColumnPresent(new_item)) return;

       vector<string> data;
       try {
            table->GetColumn(data, existing_item);
       } catch (std::exception &exc) {
            data.clear();
            for (unsigned int i = 0; i < table->GetNumRows(); ++i) data.push_back("");
       }
       table->AddColumn(new_item, data);
}

void rename_item(ISTable *table, const std::string& existing_item, const std::string& new_item)
{
       if (table->IsColumnPresent(existing_item) && !table->IsColumnPresent(new_item))
            table->RenameColumn(existing_item, new_item);
}

void get_values(ISTable *t, const unsigned int& row, const std::vector<std::vector<std::string> >& items, std::vector<std::vector<std::string> >& values)
{
       values.clear();
       if (!t) return;

       std::vector<std::string> data;
       std::string cs;
       for (std::vector<std::vector<std::string> >::const_iterator vvpos = items.begin(); vvpos != items.end(); ++vvpos) {
            data.clear();
            for (std::vector<std::string>::const_iterator vpos = vvpos->begin(); vpos != vvpos->end(); ++vpos) {
                 get_value(cs, t, row, *vpos);
                 data.push_back(cs);
            }
            values.push_back(data);
       }
}

void get_values_clean(ISTable *t, const unsigned int& row, const std::vector<std::vector<std::string> >& items, std::vector<std::vector<std::string> >& values)
{
       values.clear();
       if (!t) return;

       std::vector<std::string> data;
       std::string cs;
       for (std::vector<std::vector<std::string> >::const_iterator vvpos = items.begin(); vvpos != items.end(); ++vvpos) {
            data.clear();
            for (std::vector<std::string>::const_iterator vpos = vvpos->begin(); vpos != vvpos->end(); ++vpos) {
                 get_value_clean(cs, t, row, *vpos);
                 data.push_back(cs);
            }
            values.push_back(data);
       }
}

std::string getLowerCaseDEPID(CifFile *fobj)
{
       return __getIDs_from_database_2(fobj, "WWPDB", false);
}

std::string getUpperCaseDEPID(CifFile *fobj)
{
       return __getIDs_from_database_2(fobj, "WWPDB");
}

std::string getLowerCasePDBID(CifFile *fobj)
{
       return __getIDs_from_database_2(fobj, "PDB", false);
}

std::string getUpperCasePDBID(CifFile *fobj)
{
       return __getIDs_from_database_2(fobj, "PDB");
}

std::string getLowerCaseRCSBID(CifFile *fobj)
{
       return __getIDs_from_database_2(fobj, "RCSB", false);
}

std::string getUpperCaseRCSBID(CifFile *fobj)
{
       return __getIDs_from_database_2(fobj, "RCSB");
}

std::string getLowerCaseNDBID(CifFile *fobj)
{
       return __getIDs_from_database_2(fobj, "NDB", false);
}

std::string getUpperCaseNDBID(CifFile *fobj)
{
       return __getIDs_from_database_2(fobj, "NDB");
}

std::string getLowerCaseBMRBID(CifFile *fobj)
{
       return __getIDs_from_database_2(fobj, "BMRB", false);
}

std::string getUpperCaseBMRBID(CifFile *fobj)
{
       return __getIDs_from_database_2(fobj, "BMRB");
}

std::string getLowerCaseEMDBID(CifFile *fobj)
{
       return __getIDs_from_database_2(fobj, "EMDB", false);
}

std::string getUpperCaseEMDBID(CifFile *fobj)
{
       return __getIDs_from_database_2(fobj, "EMDB");
}

std::string getLowerCaseDEPID(Block& block)
{
       return __getIDs_from_database_2(block, "WWPDB", false);
}

std::string getUpperCaseDEPID(Block& block)
{
       return __getIDs_from_database_2(block, "WWPDB");
}

std::string getLowerCasePDBID(Block& block)
{
       return __getIDs_from_database_2(block, "PDB", false);
}

std::string getUpperCasePDBID(Block& block)
{
       return __getIDs_from_database_2(block, "PDB");
}

std::string getLowerCaseRCSBID(Block& block)
{
       return __getIDs_from_database_2(block, "RCSB", false);
}

std::string getUpperCaseRCSBID(Block& block)
{
       return __getIDs_from_database_2(block, "RCSB");
}

std::string getLowerCaseNDBID(Block& block)
{
       return __getIDs_from_database_2(block, "NDB", false);
}

std::string getUpperCaseNDBID(Block& block)
{
       return __getIDs_from_database_2(block, "NDB");
}

std::string getLowerCaseBMRBID(Block& block)
{
       return __getIDs_from_database_2(block, "BMRB", false);
}

std::string getUpperCaseBMRBID(Block& block)
{
       return __getIDs_from_database_2(block, "BMRB");
}

std::string getLowerCaseEMDBID(Block& block)
{
       return __getIDs_from_database_2(block, "EMDB", false);
}

std::string getUpperCaseEMDBID(Block& block)
{
       return __getIDs_from_database_2(block, "EMDB");
}

void get_rcsb_and_pdb_ids(CifFile *fobj, std::string &RCSB_ID, std::string &rcsb_id, std::string &PDB_ID, std::string &pdb_id, std::string &NDB_ID,
                          std::string &ndb_id, std::string &BMRB_ID, std::string &bmrb_id)
{
       Block &block = fobj->GetBlock(fobj->GetFirstBlockName());
       get_rcsb_and_pdb_ids(block, RCSB_ID, rcsb_id, PDB_ID, pdb_id, NDB_ID, ndb_id, BMRB_ID, bmrb_id);
}

void get_rcsb_and_pdb_ids(Block &block, std::string &RCSB_ID, std::string &rcsb_id, std::string &PDB_ID, std::string &pdb_id, std::string &NDB_ID,
                          std::string &ndb_id, std::string &BMRB_ID, std::string &bmrb_id)
{
       RCSB_ID.clear(); rcsb_id.clear();
       PDB_ID.clear(); pdb_id.clear();
       NDB_ID.clear(); ndb_id.clear();
       BMRB_ID.clear(); bmrb_id.clear();

       ISTable *Table = block.GetTablePtr("database_2");
       if (!is_empty_table(Table)) {
            std::string cs;
            for (unsigned int i = 0; i < Table->GetNumRows(); i++) {
                 get_value_clean_upper(cs, Table, i, "database_id");
                 if (cs == "RCSB") {
                      get_value_clean_upper(RCSB_ID, Table, i, "database_code");
                 } else if (cs == "PDB") {
                      get_value_clean_upper(PDB_ID, Table, i, "database_code");
                 } else if (cs == "NDB") {
                      get_value_clean_upper(NDB_ID, Table, i, "database_code");
                 } else if (cs == "BMRB") {
                      get_value_clean_upper(BMRB_ID, Table, i, "database_code");
                 }
            }
       }

       if (RCSB_ID.empty() && PDB_ID.empty() && NDB_ID.empty() && !BMRB_ID.empty())
           RCSB_ID = BMRB_ID;

       String::LowerCase(RCSB_ID, rcsb_id);
       String::LowerCase(PDB_ID, pdb_id);
       String::LowerCase(NDB_ID, ndb_id);
       String::LowerCase(BMRB_ID, bmrb_id);
}

void get_exp_and_mol_type(CifFile *fobj, std::string &exp_type, std::string &mol_type)
{
       exp_type = "X-ray";

       std::string cs;
  
       Block &block = fobj->GetBlock(fobj->GetFirstBlockName());
       ISTable *Table = block.GetTablePtr("exptl");
       if (!is_empty_table(Table)) {
            for (unsigned int i = 0; i < Table->GetNumRows(); i++) {
                 get_value_clean_upper(cs, Table, i, "method");
                 if (cs.find("NMR") != std::string::npos || cs.find("INFRARED SPECTROSCOPY") != std::string::npos) {
                      exp_type = "NMR";
                      break;
                 } else if (cs.find("THEORETICAL MODEL") != std::string::npos) {
                      exp_type = "model";
                      break;
                 } else if (cs.find("CRYO-ELECTRON MICROSCOPY") != std::string::npos || cs.find("ELECTRON DIFFRACTION") != std::string::npos ||
                            cs.find("ELECTRON TOMOGRAPHY") != std::string::npos || cs.find("ELECTRON MICROSCOPY") != std::string::npos) {
                      exp_type = "EM";
                      break;
                 }
            }
       }

       mol_type.clear();
       Table = block.GetTablePtr("entity_poly");
       if (!is_empty_table(Table)) {
            int has_na = 0;
            int has_prot = 0;
            for (unsigned int i = 0; i < Table->GetNumRows(); i++) {
                 get_value_clean_lower(cs, Table, i, "type");
                 if (cs.empty()) {
                      get_value_clean_upper(cs, Table, i, "ndb_seq_one_letter_code");
                      if (!cs.empty()) {
                           if (cs.find("R") != std::string::npos || cs.find("D") != std::string::npos || cs.find("Q") != std::string::npos ||
                               cs.find("E") != std::string::npos || cs.find("H") != std::string::npos || cs.find("L") != std::string::npos ||
                               cs.find("K") != std::string::npos || cs.find("M") != std::string::npos || cs.find("F") != std::string::npos ||
                               cs.find("P") != std::string::npos || cs.find("S") != std::string::npos || cs.find("T") != std::string::npos ||
                               cs.find("W") != std::string::npos || cs.find("Y") != std::string::npos || cs.find("V") != std::string::npos)
                                has_prot = 1;
                           else has_na = 1;
                      }
                 } else if (cs.find("polypeptide") != std::string::npos) {
                      has_prot = 1;
                 } else if (cs.find("ribonucleotide") != std::string::npos) {
                      has_na = 1;
                 }
            }
            if (has_prot && has_na)
                 mol_type = "naprot";
            else if (has_na)
                 mol_type = "na";
            else if (has_prot)
                 mol_type = "prot";
       }
       if (mol_type.empty()) mol_type = "prot";
}

static std::string __getIDs_from_database_2(CifFile *fobj, const std::string& database_id, const bool& is_upper_case)
{
       Block &block = fobj->GetBlock(fobj->GetFirstBlockName());
       return __getIDs_from_database_2(block, database_id, is_upper_case);
}

static std::string __getIDs_from_database_2(Block &block, const std::string& database_id, const bool& is_upper_case)
{
       ISTable *t = getTablePtr(block, "database_2");
       if (!t) return "";

       std::map<std::string, std::string> ids_mapping;
       ids_mapping.clear();

       std::string db_name, db_code;
       for (unsigned int i = 0; i < t->GetNumRows(); i++) {
            get_value_clean_upper(db_name, t, i, "database_id");
            get_value_clean_upper(db_code, t, i, "database_code");
            if (!db_name.empty() && !db_code.empty()) ids_mapping.insert(make_pair(db_name, db_code));
       }
/*
       db_code.clear();
       std::map<std::string, std::string>::const_iterator mpos = ids_mapping.find("BMRB");
       if (mpos != ids_mapping.end()) db_code = mpos->second;
       if (ids_mapping.find("RCSB") == ids_mapping.end() && ids_mapping.find("PDB") == ids_mapping.end() &&
           ids_mapping.find("NDB") == ids_mapping.end() && !db_code.empty())
            ids_mapping.insert(make_pair("RCSB", db_code));

       db_code.clear();
       mpos = ids_mapping.find("RCSB");
       if (mpos != ids_mapping.end()) db_code = mpos->second;
       if (ids_mapping.find("WWPDB") == ids_mapping.end() && !db_code.empty())
            ids_mapping.insert(make_pair("WWPDB", db_code));
*/
       db_code.clear();
       std::map<std::string, std::string>::const_iterator mpos = ids_mapping.find(database_id);
       if (mpos != ids_mapping.end()) {
            db_code = mpos->second;
            if (!is_upper_case) String::LowerCase(db_code);
       }

       return db_code;
}
