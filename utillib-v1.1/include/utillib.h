/*
FILE:     utillib.h
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
#ifndef _H_UTIL_TOOL_H_
#define _H_UTIL_TOOL_H_

#include <stdio.h>

#include <list>
#include <set>
#include <string>
#include <map>
#include <vector>

#include "CifFile.h"
#include "LogUtil.h"

#define MASTER_TO_INTERNAL 1
#define INTERNAL_TO_MASTER 2

#define CS_EXTENSION  "-cs.cif"


extern void add_missing_zero_before_point(std::string& f);
extern ISTable *add_new_table(const std::string &tablename, const int& num_items, const char** items);
extern ISTable *add_new_table(const std::string &tablename, const std::vector<std::string> &items);
extern ISTable *add_new_table(const std::string &tablename, const std::vector<std::string> &items,
                              const std::vector<std::vector<std::string> >& values);
extern void add_type(std::map<std::string, int>& all_types, const std::string& type);
extern void capture_upload_file(const std::string &fileid, const std::string &outfilename,
         const std::string &outextension, const std::string &dirPath,
         const std::string &HT_ZCAT_COMMAND, const std::string &HT_GZCAT_COMMAND,
         const std::string &HT_ASCII_FILTER_COMMAND);
extern void check_missing_item(ISTable *table, const int& num_items, const char** items, const bool mandatary = false);
extern void check_missing_item(ISTable *table, const std::vector<std::string> &items);
extern void check_path_variable(std::string &path);
extern bool checking_regular_expression_match_ok(const std::string& pattern, const std::string& value);
/*
extern void clustering_with_insertion(std::vector<std::set<int> >& groups,
                       std::vector<std::pair<int, int> >& links);
extern void clustering_with_merging(std::vector<std::set<int> >& groups,
                       std::vector<std::pair<int, int> >& links);
*/
extern void combine_values_delimit_by_comma(std::string &value1, const std::string &value2);
extern void convert_time_format(const int& length, std::string& datestring);
extern void copy_file(const std::string &old_file, const std::string &new_file);
extern void copy_item(ISTable *table, const std::string& existing_item, const std::string& new_item);
extern int count_pattern(const std::string& line, const std::string& pattern);
extern void create_directory(const std::string& top_dirname, const std::string& sub_dirname);
extern CifFile *create_fobj(const std::string& odbfile, const std::string& blockId, const bool& caseSensitive = false, const int& length = 250);
extern void delete_space_between_names(std::string& namestring);
extern void deleteTable(Block &block, const std::string& table_name);
extern int display_raw_file(FILE *io, const std::string &fName);
extern int file_corres_file(const std::string &id, std::string& filename);
extern  void find_rcsb_directory(const std::string &id, std::string& dirname);
extern int find_rcsb_file(const std::string &id, std::string &filename);
extern int find_rcsb_file(const std::string &id, std::string& filename, std::string &expfile, std::string& csfile);
extern void find_rings(const unsigned int& size, const std::map<unsigned int, std::vector<unsigned int> >& neighbor_map,
                       std::vector<std::vector<unsigned int> >& ringLists);
extern std::string FloatToString(const double&);
extern std::string FloatToString(const double&, const int&, const int&, const bool& left_adjust = false, const bool& forcedWidthLimit = false);
extern std::string FormattedString(const std::string&, const int&, const bool& left_adjust = false, const bool& forcedWidthLimit = false);
extern std::string FormattedFieldValue(const std::string& value, const int& FieldType, const int& FieldWidth, const int& FieldPrec,
                                       const bool& left_adjust, const bool& forcedWidthLimit = false);
extern void get_address_book(const std::string&, std::map<string, std::vector<std::string> >&);
extern void get_chain_start_and_end(const std::set<unsigned int>& chain_break_set, unsigned int size, 
                                    std::vector<std::pair<unsigned int, unsigned int> >& pair_array);
extern void get_conditional_value(Block& block, std::string& value, const std::string& category, const std::string& item,
                                  const std::string& condition_item, const std::string& condition_value);
extern void get_current_directory_path(std::string &dirname);
extern std::string getLowerCaseDEPID(CifFile *fobj);
extern std::string getUpperCaseDEPID(CifFile *fobj);
extern std::string getLowerCasePDBID(CifFile *fobj);
extern std::string getUpperCasePDBID(CifFile *fobj);
extern std::string getLowerCaseRCSBID(CifFile *fobj);
extern std::string getUpperCaseRCSBID(CifFile *fobj);
extern std::string getLowerCaseNDBID(CifFile *fobj);
extern std::string getUpperCaseNDBID(CifFile *fobj);
extern std::string getLowerCaseBMRBID(CifFile *fobj);
extern std::string getUpperCaseBMRBID(CifFile *fobj);
extern std::string getLowerCaseEMDBID(CifFile *fobj);
extern std::string getUpperCaseEMDBID(CifFile *fobj);
extern std::string getLowerCaseDEPID(Block &block);
extern std::string getUpperCaseDEPID(Block &block);
extern std::string getLowerCasePDBID(Block &block);
extern std::string getUpperCasePDBID(Block &block);
extern std::string getLowerCaseRCSBID(Block &block);
extern std::string getUpperCaseRCSBID(Block &block);
extern std::string getLowerCaseNDBID(Block &block);
extern std::string getUpperCaseNDBID(Block &block);
extern std::string getLowerCaseBMRBID(Block &block);
extern std::string getUpperCaseBMRBID(Block &block);
extern std::string getLowerCaseEMDBID(Block &block);
extern std::string getUpperCaseEMDBID(Block &block);
extern void get_date(std::string &date_string);
extern void get_date_from_time(std::string &date_string, const time_t time_reading);
extern void get_exp_and_mol_type(CifFile *fobj, std::string &exp_type, std::string &mol_type);
extern void get_file_content(const std::string &fName, std::string &content);
extern void get_file_list_from_directory(std::vector<std::string> &file_list, const std::string &dir_name);
extern void get_full_date(std::string &full_date, const std::string pdb_date);
extern int  get_line_from_file(FILE *fp, std::string &str);
extern void get_line_from_file(FILE *fp, std::string &str, std::vector<std::string>& non_ascii_positions);
extern void get_local_file_copy(const std::string &input_filename, std::string &filename_without_path, std::string &filename_of_local_copy);
extern void get_max_length_words(std::vector<std::string> &words, const std::string &line, const unsigned int max_length);
extern std::string get_next_available_asym_id(std::set<std::string>& used_ids);
extern std::string get_next_available_pdb_chain_id(std::set<std::string>& used_ids);
extern void get_next_friday(std::string &date_string);
extern std::string get_next_id(const std::string& firstcolumn, const std::string& secondcolumn, const std::string& thirdcolumn, std::set<std::string>& id_set);
extern void get_next_thursday(std::string &date_string);
extern void get_next_tuesday(std::string &date_string);
extern void get_next_wednesday(std::string &next_wednesday);
extern void get_no_loop_category_set(const std::string& filePath, std::set<std::string>& no_loop_categories);
extern void get_one_line(FILE *fp, std::string &cs);
extern void get_order_list(const std::string& filename, Block& block, std::vector<std::string>& order_list);
extern std::string getPDBSymmetry(const std::string& SymOP);
extern void get_rcsb_and_pdb_ids(CifFile *fobj, std::string &RCSB_ID, std::string &rcsb_id, std::string &PDB_ID, std::string &pdb_id, std::string &NDB_ID,
                                 std::string &ndb_id, std::string &BMRB_ID, std::string &bmrb_id);
extern void get_rcsb_and_pdb_ids(Block &block, std::string &RCSB_ID, std::string &rcsb_id, std::string &PDB_ID, std::string &pdb_id, std::string &NDB_ID,
                                 std::string &ndb_id, std::string &BMRB_ID, std::string &bmrb_id);
extern std::set<std::string> get_set_common_element(const std::set<std::string>& set1, const std::set<std::string>& set2);
extern std::string getSingleBestValue(Block& block, const std::string& category, const std::string& item, const int& number_compare_flag);
extern std::set<std::string> get_stop_words(const std::string& filename);
extern ISTable *getTableCopy(Block &block, const std::string& table_name);
extern ISTable *getTablePtr(Block &block, const std::string& table_name);
extern void get_temp_directory(std::string &temp_dirname, const std::string& dirPath);
extern void get_temp_filename(std::string &temp_filename);
extern void get_temp_filename(std::string &temp_filename, const std::string& dirPath);
extern std::string get_type(const std::map<std::string, int>& all_types, const std::string& d_value = "HETAIN");
extern std::string get_reformated_text_with_width(const std::string& input_text, const unsigned int& width, const std::vector<std::string>& delimList,
                                                  const std::string& next_line_pre_empty_space);
extern void get_text_array_from_block(std::vector<std::string>& text_array, const std::string& block_text, const int& width);
extern void get_text_array_from_block_with_prefix(std::vector<std::string>& text_array, const std::string& block_text, const std::string& prefix,
                                                  const int& width);
extern std::string get_unique_atom_name(std::set<std::string>& existing_name_set, const std::string& name, const std::string& type);
extern void get_unique_values(ISTable *Table, const std::string& colName, std::set<std::string>& value_set);
extern void get_validation_letter_template(const std::string&, const bool, std::vector<std::vector<std::string> >&, std::map<string, std::string>&,
                                           std::vector<std::vector<std::string> >&); 
extern int get_value(std::string &cifstring, ISTable *Table, const unsigned int& row, const std::string& colName);
extern int get_value_clean(std::string &cifstring, ISTable *Table, const unsigned int& row, const std::string& colName);
extern int get_value_clean_lower(std::string &cifstring, ISTable *Table, const unsigned int& row, const std::string& colName);
extern int get_value_clean_upper(std::string &cifstring, ISTable *Table, const unsigned int& row, const std::string& colName);
extern int get_value_lower(std::string &cifstring, ISTable *Table, const unsigned int& row, const std::string& colName);
extern int get_value_upper(std::string &cifstring, ISTable *Table, const unsigned int& row, const std::string& colName);
extern void get_values(ISTable *Table, std::vector<std::vector<std::string> >& values);
extern void get_values(ISTable *Table, std::vector<std::map<std::string, std::string> >& values);
extern void get_values(ISTable *Table, std::list<std::map<std::string, std::string> >& values);
extern void get_values(ISTable *Table, const std::vector<std::string>& itemNames, std::vector<std::vector<std::string> >& values);
extern void get_values(ISTable *t, const unsigned int& row, const std::vector<std::vector<std::string> >& items, std::vector<std::vector<std::string> >& values);
extern void get_values_clean(ISTable *Table, const std::vector<std::string>& itemNames, std::vector<std::vector<std::string> >& values);
extern void get_values_clean(ISTable *t, const unsigned int& row, const std::vector<std::vector<std::string> >& items, std::vector<std::vector<std::string> >& values);
extern string::size_type get_word(const std::string &line, const std::string& delims, string::size_type begIdx, std::string &word);
extern void get_wordarray(std::vector<std::string> &wordarray, const std::string &line, const std::string& delims);
extern void get_wordarray_delimit_by_string(std::vector<std::string> &wordarray, const std::string &line, const std::string &delimit);
extern void get_wordlist_from_string_separated_by_delimits(std::vector<std::string>& wordlist, const std::string& string_buff, const std::string& delimit,
                                                           const std::string& erase_string = "", const bool& with_empty = false);
extern void get_wordarray_with_space(std::vector<std::string> &wordarray, const std::string &line, const std::string& delims);
extern void get_wordset(std::set<std::string> &wordset, const std::string &line, const std::string& delims);
extern int is_empty_row(ISTable *Table, const std::vector<std::string>& ColumnNames, const unsigned int& rowindex, const std::string& except);
extern int is_empty_row(ISTable *Table, const std::vector<std::string>& ColumnNames, const unsigned int& rowindex, const std::set<std::string>& except);
extern int is_empty_row_loose(ISTable *Table, const std::vector<std::string>& ColumnNames, const unsigned int& rowindex);
extern int is_empty_table(ISTable *Table);
extern int is_empty_table(ISTable *Table, const std::set<std::string>& except);
extern int is_empty_table_loose(ISTable *Table);
extern int is_empty_table(ISTable *Table, const char *except);
extern bool is_same_map(const std::map<std::string, std::string>& map1, const std::map<std::string, std::string>& map2);
extern bool is_same_set(const std::set<std::string>& set1, const std::set<std::string>& set2);
extern int is_wrong_date(const std::string& date_string);
extern bool IsAlnum(const std::string& value);
extern bool IsInteger(const std::string& sValue, int& iValue);
extern bool IsResidueName(const std::string& value);
extern bool IsSameStringCombination(const std::string& string1, const std::string& string2, const std::string& delims);
extern bool IsNotValidatedValue(const std::string& value);
extern std::string join_string(const std::vector<std::string>& data, const std::string& delimiter);
extern std::string join_string(const std::set<std::string>& data, const std::string& delimiter);
extern void master_internal_conversion(const std::string&, Block&, const int&);
extern void _MultipleSort(const std::vector<std::vector<int> >& score_matrix, std::vector<int>& index);
extern std::string numerical_multiplier(const unsigned int& number);
extern std::string other_numerical_multiplier(const unsigned int& number);
extern CifFile* read_cif_file(const std::string &odbfile, const std::string &textfile, std::string &diags, const bool& caseSensitive = false,
                              const int& length = 250);
extern CifFile *get_fobj(std::string& error, const std::string& ciffile, const bool& caseSensitive = false, const int& length = 250);
extern CifFile *get_fobj(FILE *log, const std::string& odbfile, const std::string& ciffile, const bool& caseSensitive = false, const int& length = 250);
extern CifFile *get_fobj(LogUtil *log, const std::string& odbfile, const std::string& ciffile, const bool& caseSensitive = false, const int& length = 250);
extern CifFile *get_fobj(FILE *log, const std::string& cifString, const bool& caseSensitive = false, const int& length = 250);
extern int parseString(const std::string &orginal_string, std::vector<std::string> &list);
extern void preprocess_file(const std::string& infile, const std::string& outfile);
extern std::string PrintSpace(const int& width);
extern std::string read_last_line_from_file(const std::string& filename);
extern void reformatting_float_number_without_extra_zero(std::string& f);
extern void reformatting_charge_to_prefix(std::string& charge);
extern void reformatting_charge_to_suffix(std::string& charge);
extern void replace_semicolons(std::string &str);
extern std::string replace_string(const std::string& in_string, const std::vector<std::string>& data);
extern void replace_string(std::string &line, const std::string &removed_word, const std::string &added_word);
extern int replace_string_new_to_old(std::string &line, const std::string &added_word, const std::string &removed_word);
extern void RemoveNonAlnum(std::string& value, const bool& sign_flag = false);
extern void remove_directory(const std::string &dir_name);
extern void remove_redundant_keywords(const std::string& filename, std::string& keyword, std::string& header);
extern void rename_item(ISTable *table, const std::string& existing_item, const std::string& new_item);
extern bool reset_entry_id(Block& block, const std::string& entry_id, const bool& include_entry = true);
extern void update_table(ISTable* t, const std::vector<std::map<std::string, std::string> >& data);
extern void update_table(ISTable* t, const std::vector<std::string>& items, const std::vector<std::vector<std::string> >& values);
extern bool update_value(ISTable* t, const std::string& item, const unsigned int& row, const std::string& value, const bool& force_update=false);
extern void write_cif_file(const std::string& blockID, const std::string& filename, ISTable* table);
extern void write_cif_file(const std::string& blockID, const std::string& filename, const std::vector<ISTable*>& tables);
extern void write_cif_file(const std::string& filename, const std::list<std::pair<std::string, std::vector<ISTable*> > >& bt_lists);
extern void write_json_map_list(FILE* fp, const std::list<std::map<std::string, std::string> >& map_list, const bool& has_value_item_only = false);
extern void write_json_map_list(FILE* fp, const std::vector<std::map<std::string, std::string> >& map_list, const bool& has_value_item_only = false);
// extern void write_json_map(FILE* fp, const std::map<std::string, std::map<std::string, std::string> >& map_value, const bool& has_value_item_only = false);
extern void write_json_map(FILE* fp, const std::map<std::string, std::string>& map_value, const bool& has_value_item_only = false);
extern void write_json_open_dictionary(FILE* fp);
extern void write_json_close_dictionary(FILE* fp);
extern void write_json_open_list(FILE* fp);
extern void write_json_close_list(FILE* fp);
extern void write_json_dictionry_key(FILE* fp, const std::string& key);
extern void write_json_dictionry_key_value_pair(FILE* fp, const std::string& key, const std::string& value);
extern void write_json_comma_delimit(FILE* fp, const bool& need_delimit);
extern void write_json_string_vector(FILE* fp, const std::vector<std::string>& values);
extern void write_json_string_vectors(FILE* fp, const std::vector<std::vector<std::string> >& values);
extern void write_json_string_value(FILE* fp, const std::string& value);
extern void write_json_int_vector(FILE* fp, const std::vector<int>& values);
extern void write_json_int_vectors(FILE* fp, const std::vector<std::vector<int> >& values);
extern int wsystem(const char *cmdstring, const char *path, const int time_limit = 0);

#endif
