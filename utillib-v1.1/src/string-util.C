/*
FILE:     string-util.C
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
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <ctype.h>

#include "regex.h"
#include "utillib.h"

#define  NS               10

static bool IsDashValue(const std::string& value);
static std::string get_next_available_id(const std::string& letters, std::set<std::string>& used_ids);
// static std::string get_next_available_id(const std::string& letters, const std::set<std::string>& used_ids, std::vector<std::string>& input_list);

std::string read_last_line_from_file(const std::string& filename)
{
       std::string lastLine = "";

       std::ifstream fin(filename.c_str());
       if (fin) {
            fin.seekg(0, std::ios_base::end);         // Start at end of file
            char ch = ' ';                            // Init ch not equal to '\n'
            while (ch != '\n') {
                 fin.seekg(-2, std::ios_base::cur);   // Two steps back, this means we will NOT check the last character
                 if ((int)fin.tellg() <= 0) {         // If passed the start of the file,
                      fin.seekg(0);                   // this is the start of the line
                      break;
                 }
                 fin.get(ch);                         // Check the next character
            }
            getline(fin, lastLine);                   // Read the current line
            fin.close();
       }

       return lastLine;
}

bool checking_regular_expression_match_ok(const std::string& pattern, const std::string& value)
{
       regex_t preg;
       regmatch_t pmatch[NS];
       rcsb_regcomp(&preg, pattern.c_str(), REG_EXTENDED);

       bool wrong_pattern = false;
       int ret = rcsb_regexec(&preg, value.c_str(), NS, pmatch, 0);
       if (ret != 0) wrong_pattern = true;
       else {
            int len = pmatch[0].rm_eo - pmatch[0].rm_so;
            if (len != (int) value.size()) wrong_pattern = true;
       }

       rcsb_regfree(&preg);

       return (!wrong_pattern);
}

bool is_same_map(const std::map<std::string, std::string>& map1, const std::map<std::string, std::string>& map2)
{
       if (map1.size() != map2.size()) return false;

       for (std::map<std::string, std::string>::const_iterator mpos1 = map1.begin(); mpos1 != map1.end(); ++mpos1) {
            std::map<std::string, std::string>::const_iterator mpos2 = map2.find(mpos1->first);
            if (mpos2 == map2.end()) return false;
            if (mpos1->second != mpos2->second) return false;
       }

       for (std::map<std::string, std::string>::const_iterator mpos2 = map2.begin(); mpos2 != map2.end(); ++mpos2) {
            std::map<std::string, std::string>::const_iterator mpos1 = map1.find(mpos2->first);
            if (mpos1 == map1.end()) return false;
            if (mpos1->second != mpos2->second) return false;
       }
       return true;
}

bool is_same_set(const std::set<std::string>& set1, const std::set<std::string>& set2)
{
       if (set1.size() != set2.size()) return false;

       for (std::set<std::string>::const_iterator spos1 = set1.begin(); spos1 != set1.end(); ++spos1) {
            std::set<std::string>::const_iterator spos2 = set2.find(*spos1);
            if (spos2 == set2.end()) return false;
       }

       for (std::set<std::string>::const_iterator spos2 = set2.begin(); spos2 != set2.end(); ++spos2) {
            std::set<std::string>::const_iterator spos1 = set1.find(*spos2);
            if (spos1 == set1.end()) return false;
       }
       return true;
}

void get_max_length_words(std::vector<std::string> &words, const std::string &line, const unsigned int max_length)
{
       words.clear();
       if (line.empty()) return;

       if (line.length() <= max_length) {
            words.push_back(line);
            return;
       }

       unsigned int start = 0;
       unsigned int end = 0;
       for (unsigned int i = 0; i < line.size(); ++i) {
            if (!isalnum(line[i]) && line[i] != '.' &&
                 line[i] != '(' && line[i] != '[') {
                 if ((i - start + 1) > max_length) {
                      if (end < start) {
                           end = i;
                           if ((end - start + 1) > max_length) end = start + max_length - 1;
                      }

                      std::string cs = line.substr(start, end - start + 1);
                      String::StripLeadingWs(cs);
                      String::StripTrailingWs(cs);
                      words.push_back(cs);
                      start = end + 1;
                 } else end = i; 
            }
       }
       if ((line.size() - start) > max_length && end > 0 &&
           (end - start + 1) <= max_length) {
            std::string cs = line.substr(start, end - start + 1);
            String::StripLeadingWs(cs);
            String::StripTrailingWs(cs);
            words.push_back(cs);
            start = end + 1;
       }

       if (start < line.size()) {
            std::string cs = line.substr(start, line.size() - start);
            String::StripLeadingWs(cs);
            String::StripTrailingWs(cs);
            words.push_back(cs);
       }
}

std::string::size_type get_word(const std::string &line, const std::string& delims, std::string::size_type begIdx, std::string &word)
{
       std::string::size_type endIdx;

       word.clear();
       if (begIdx != std::string::npos) {
            endIdx = line.find_first_of(delims, begIdx);
            if (endIdx == std::string::npos) {
                 endIdx = line.length();
            }
            word = line.substr(begIdx, endIdx - begIdx);
            begIdx = line.find_first_not_of(delims, endIdx);
       }
       return begIdx;
}

void get_wordarray(std::vector<std::string> &wordarray, const std::string &line, const std::string& delims)
{
       wordarray.clear();
       if (line.empty()) return;

       std::string word;
       std::string::size_type begIdx, endIdx;

       begIdx = line.find_first_not_of(delims);
       while (begIdx != std::string::npos) {
            endIdx = line.find_first_of(delims, begIdx);
            if (endIdx == std::string::npos) {
                 endIdx = line.length();
            }
            word = line.substr(begIdx, endIdx - begIdx);
            if (word == "NULL") word.clear();
            String::StripAndCompressWs(word);
            wordarray.push_back(word);
            begIdx = line.find_first_not_of(delims, endIdx);
       }
}

void get_wordarray_with_space(std::vector<std::string> &wordarray, const std::string &line, const std::string& delims)
{
       wordarray.clear();
       if (line.empty()) return;

       std::string::size_type begIdx = line.find(delims);
       if (begIdx != std::string::npos) {
            std::string word = line.substr(0, begIdx);
            wordarray.push_back(word);
            while (begIdx != std::string::npos) {
                 std::string::size_type endIdx = line.find(delims, begIdx + delims.size());
                 if (endIdx == std::string::npos) endIdx = line.length();
                 word = line.substr(begIdx + delims.size(), endIdx - begIdx -
                                delims.size());
                 wordarray.push_back(word);
                 endIdx = line.find(delims, begIdx + delims.size());
                 begIdx = endIdx;
            }
       } else wordarray.push_back(line);
}

void get_wordarray_delimit_by_string(std::vector<std::string> &wordarray, const std::string &line, const std::string &delimit)
{
       wordarray.clear();
       if (line.empty()) return;

       std::string::size_type begIdx = line.find(delimit);
       if (begIdx != std::string::npos) {
            std::string word = line.substr(0, begIdx);
            String::StripAndCompressWs(word);
            wordarray.push_back(word);
            while (begIdx != std::string::npos) {
                 std::string::size_type endIdx = line.find(delimit, begIdx + delimit.size());
                 if (endIdx == std::string::npos) endIdx = line.length();
                 word = line.substr(begIdx + delimit.size(), endIdx - begIdx - delimit.size());
                 String::StripAndCompressWs(word);
                 wordarray.push_back(word);
                 endIdx = line.find(delimit, begIdx + delimit.size());
                 begIdx = endIdx;
            }
       } else wordarray.push_back(line);
}

void get_wordlist_from_string_separated_by_delimits(std::vector<std::string>& wordlist, const std::string& string_buff, const std::string& delimit,
                                                    const std::string& erase_string, const bool& with_empty)
{
       wordlist.clear();
       std::vector<std::string> data;
       get_wordarray_delimit_by_string(data, string_buff, delimit);
       for (unsigned int i = 0; i < data.size(); ++i) {
            if (!erase_string.empty()) {
                 std::string::size_type p = data[i].find(erase_string);
                 if (p != std::string::npos) data[i].erase(p);
            }
            String::StripAndCompressWs(data[i]);
            if (data[i].empty() || data[i] == "NULL" || data[i] == "N/A") {
                 if (with_empty) wordlist.push_back("");
                 continue;
            }
            wordlist.push_back(data[i]);
       }
}

void get_wordset(std::set<std::string> &wordset, const std::string &line, const std::string& delims)
{
       wordset.clear();
       if (line.empty()) return;

       std::string word;
       std::string::size_type begIdx, endIdx;

       begIdx = line.find_first_not_of(delims);
       while (begIdx != std::string::npos) {
            endIdx = line.find_first_of(delims, begIdx);
            if (endIdx == std::string::npos) {
                 endIdx = line.length();
            }
            word = line.substr(begIdx, endIdx - begIdx);
            String::StripAndCompressWs(word);
            wordset.insert(word);
            begIdx = line.find_first_not_of(delims, endIdx);
       }
}

std::string join_string(const std::vector<std::string>& data, const std::string& delimiter)
{
       std::string join_cs = "";
       for (std::vector<std::string>::const_iterator pos = data.begin(); pos != data.end(); ++pos) {
            if (!join_cs.empty()) join_cs += delimiter;
            join_cs += *pos; 
       }
       return join_cs;
}

std::string join_string(const std::set<std::string>& data, const std::string& delimiter)
{
       std::string join_cs = "";
       for (std::set<std::string>::const_iterator pos = data.begin(); pos != data.end(); ++pos) {
            if (!join_cs.empty()) join_cs += delimiter;
            join_cs += *pos; 
       }
       return join_cs;
}

void replace_semicolons(std::string &str)
{
       if (str.empty()) return;

       for (unsigned int i = 0; i < str.size() - 1; i++) {
            if (str[i] == '\n' && str[i + 1] == ';') {
                str[i] = ';'; str[i + 1] = '\n';
            }
       }
}

std::string replace_string(const std::string& in_string, const std::vector<std::string>& data)
{
       std::string line = in_string;
       for (std::vector<std::string>::const_iterator pos = data.begin();pos != data.end(); ++pos) {
            std::string::size_type size_pos = line.find("%s");
            if (size_pos != std::string::npos) line.replace(size_pos, 2, *pos);
       }
       return line;
}

void replace_string(std::string &line, const std::string &removed_word, const std::string &added_word)
{
       std::string::size_type pos = line.find(removed_word);
       while (pos != std::string::npos) {
            line.replace(pos, removed_word.size(), added_word);
            pos = line.find(removed_word, pos + added_word.size());
       }
}


int replace_string_new_to_old(std::string &line, const std::string &added_word, const std::string &removed_word)
{
       int count = 0;
       std::string::size_type pos = line.find(removed_word);
       while (pos != std::string::npos) {
            count++;
            line.replace(pos, removed_word.size(), added_word);
            pos = line.find(removed_word, pos + added_word.size());
       }
       return count;
}

int count_pattern(const std::string& line, const std::string& pattern)
{
       int count = 0;
       std::string::size_type pos = line.find(pattern);
       while (pos != std::string::npos) {
            count++;
            pos = line.find(pattern, pos + pattern.size());
       }
       return count;
} 

void combine_values_delimit_by_comma(std::string &value1, const std::string &value2)
{
       if (value1.empty()) {
            value1 = value2;
            return;
       }

       std::map<std::string, std::string> mapping;
       mapping.clear();
       std::vector<std::string> data;
       std::string cs;
       get_wordarray(data, value1, ",");
       for (std::vector<std::string>::const_iterator pos = data.begin(); pos != data.end(); ++pos) {
            String::UpperCase(*pos, cs);
            mapping.insert(std::make_pair(cs, *pos));
       }
       get_wordarray(data, value2, ",");
       for (std::vector<std::string>::const_iterator pos = data.begin(); pos != data.end(); ++pos) {
            String::UpperCase(*pos, cs);
            mapping.insert(std::make_pair(cs, *pos));
       }
       value1.clear();
       for (std::map<std::string, std::string>::const_iterator mpos = mapping.begin(); mpos != mapping.end(); ++mpos) {
            if (!value1.empty()) value1 += ", ";
            value1 += mpos->second;
       }
}

void delete_space_between_names(std::string& namestring)
{
       if (namestring.empty()) return;

       std::vector<std::string> data;
       get_wordarray(data, namestring, ",");
       namestring.clear();
       for (std::vector<std::string>::const_iterator
            pos = data.begin(); pos != data.end(); ++pos) {
            if (!namestring.empty()) namestring += ",";
            namestring += *pos;
       }
}

std::string FloatToString(const double& val)
{
       char buffer[20];
       sprintf(buffer, "%.10f", val);
       std::string cs = buffer;
       return cs;
}

std::string FloatToString(const double& val, const int& w, const int& p, const bool& left_adjust, const bool& forcedWidthLimit)
{
       double val1 = val;
       double cut_off = 5.0 / pow(10.0, p + 1);
       if (fabs(val) < cut_off) val1 = fabs(val);
       std::ostringstream oss;
       if (left_adjust) {
           if (w)
                oss << std::left << std::setw(w) << std::fixed << std::setprecision(p) << val1;
           else oss << std::left << std::fixed << std::setprecision(p) << val1;
       } else {
           if (w)
                oss << std::setw(w) << std::fixed << std::setprecision(p) << val1;
           else oss << std::fixed << std::setprecision(p) << val1;
       }
       std::string rs = oss.str();
       if (forcedWidthLimit && w && (int) rs.size() > w) {
            if (left_adjust)
                 return rs.substr(0, w);
            else return rs.substr(rs.size() - w);
       }
       return (oss.str());
}

std::string FormattedString(const std::string& cs, const int& w, const bool& left_adjust, const bool& forcedWidthLimit)
{
       std::ostringstream oss;
       if (left_adjust)
            oss << std::left << std::setw(w) << cs;
       else oss << std::setw(w) << cs;
       std::string rs = oss.str();
       if (forcedWidthLimit && w && (int) rs.size() > w) {
            if (left_adjust)
                 return rs.substr(0, w);
            else return rs.substr(rs.size() - w);
       }
       return (oss.str());
}

std::string PrintSpace(const int& width)
{
       std::string cs = "";
       for (int i = 0; i < width; ++i) cs += " ";
       return cs;
}

std::string FormattedFieldValue(const std::string& value, const int& FieldType, const int& FieldWidth, const int& FieldPrec,
                                const bool& left_adjust, const bool& forcedWidthLimit)
{
       if (value.empty()) return PrintSpace(FieldWidth);

       if (FieldType > 2) return FormattedString(value, FieldWidth, left_adjust, forcedWidthLimit);
       else {
            double temp_double = atof(value.c_str());
            if (FieldType == 2) {
                 int len = 1;
                 std::string::size_type p = value.find(".");
                 if (p != std::string::npos) len = value.substr(p).size();
                 if (len > FieldPrec) temp_double += 1.0 / pow(10.0, len);
            }
            return FloatToString(temp_double, FieldWidth, FieldPrec, left_adjust, forcedWidthLimit);
       }
}

bool IsNotValidatedValue(const std::string& value)
{
       std::string cs;
       String::LowerCase(value, cs);

       if (cs == "n/a") return true;
       else if (IsDashValue(cs)) return true;

       return false;
}

bool IsAlnum(const std::string& value)
{
       if (value.empty()) return true;

       for (unsigned int i = 0; i < value.size(); ++i) {
            if (!isalnum(value[i])) return false;
       }

       return true;
}

bool IsResidueName(const std::string& value)
{
       if (IsAlnum(value)) return true;

       if (value[0] != '~' && value[0] != '%') return false;
       for (unsigned int i = 1; i < value.size(); ++i) {
            if (!isdigit(value[i])) return false;
       }
       return true;
}

void RemoveNonAlnum(std::string& value, const bool& sign_flag)
{
       std::string buff;
       buff.clear();
       for (unsigned int i = 0; i < value.size(); ++i) {
            if (isalnum(value[i])) buff.push_back(value[i]);
            else if (sign_flag && (value[i] == '+' || value[i] == '-')) buff.push_back(value[i]);
       }
       value = buff;
}

bool IsSameStringCombination(const std::string& string1, const std::string& string2, const std::string& delims)
{
       if (string1 == string2) return true;

       std::set<std::string> str_set1, str_set2;
       get_wordset(str_set1, string1, delims);
       get_wordset(str_set2, string2, delims);

       if (str_set1 == str_set2) return true;

       for (std::set<std::string>::const_iterator spos = str_set1.begin(); spos != str_set1.end(); ++spos) {
            if (str_set2.find(*spos) != str_set2.end()) return true;
       }

       return false;
}

std::string get_reformated_text_with_width(const std::string& input_text, const unsigned int& width, const std::vector<std::string>& delimList,
                                           const std::string& next_line_pre_empty_space)
{
       if (input_text.size() <= width) return input_text;

       std::map<std::string, unsigned int> delimPositionMap;
       std::map<unsigned int, std::string> positionDelimMap;

       std::string output_text = "";
       std::string text = input_text;
       unsigned int adjust_width = width;

       while (!text.empty()) {
            if (!output_text.empty()) adjust_width = width - next_line_pre_empty_space.size();
            if (text.size() <= width) {
                 output_text += "\n" + next_line_pre_empty_space + text;
                 break;
            } 

            delimPositionMap.clear();
            for (std::vector<std::string>::const_iterator vpos = delimList.begin(); vpos != delimList.end(); ++vpos) {
                 std::string::size_type begIdx = text.find(*vpos);
                 while (begIdx != std::string::npos) {
                      unsigned int position = begIdx;
                      if ((*vpos) != " ") position = begIdx + vpos->size();
                      if (position > adjust_width) break;

                      std::map<std::string, unsigned int>::iterator mpos = delimPositionMap.find(*vpos);
                      if (mpos != delimPositionMap.end()) mpos->second = position;
                      else delimPositionMap.insert(std::make_pair(*vpos, position));

                      begIdx = text.find(*vpos, begIdx + vpos->size());
                 }
            }

            if (delimPositionMap.empty()) {
                 if (!output_text.empty())
                      output_text += "\n" + next_line_pre_empty_space + text.substr(0, adjust_width);
                 else output_text = text.substr(0, adjust_width);
                 text = text.substr(adjust_width);
            } else {
                 positionDelimMap.clear();
                 for (std::map<std::string, unsigned int>::const_iterator mpos = delimPositionMap.begin(); mpos != delimPositionMap.end(); ++mpos) {
                      positionDelimMap.insert(std::make_pair(mpos->second, mpos->first));
                 }
                 std::map<unsigned int, std::string>::const_reverse_iterator rmpos = positionDelimMap.rbegin();
                 if (!output_text.empty())
                      output_text += "\n" + next_line_pre_empty_space + text.substr(0, rmpos->first);
                 else output_text = text.substr(0, rmpos->first);
                 if (rmpos->second == " ")
                      text = text.substr(rmpos->first + 1);
                 else text = text.substr(rmpos->first);
            }
       }

       return output_text;
}

void get_text_array_from_block(std::vector<std::string>& text_array, const std::string& block_text, const int& width)
{
       bool has_exceeded_width_line = false;
       get_wordarray_with_space(text_array, block_text, "\n");
       int start = -1;
       int end = -1;
       for (unsigned int i = 0; i < text_array.size(); ++i) {
            String::StripTrailingWs(text_array[i]);
            if (!text_array[i].empty()) {
                 if (start < 0) start = i;
                 end = i;
            }
            if ((int) text_array[i].size() > width) {
                 has_exceeded_width_line = true;
                 break;
            }
       }
       if (start < 0 || end < 0) {
            text_array.clear();
            return;
       }

       if (!has_exceeded_width_line) {
            if (start > 0 || end < ((int) text_array.size() - 1)) {
                 std::vector<std::string> data;
                 data.clear();
                 for (int i = start; i <= end; ++i) {
                      data.push_back(text_array[i]);
                 }
                 text_array = data;
            }
            return;
       }

       text_array.clear();

       std::string cs = block_text;
       String::StripAndCompressWs(cs);
       get_max_length_words(text_array, cs, width);
}

void get_text_array_from_block_with_prefix(std::vector<std::string>& text_array, const std::string& block_text, const std::string& prefix, const int& width)
{
       text_array.clear();

       std::string space = PrintSpace(prefix.size());
       std::vector<std::string> tmp_array;
       get_text_array_from_block(tmp_array, block_text, width);
       for (unsigned int i = 0; i < tmp_array.size(); ++i) {
            if (i == 0)
                 text_array.push_back(prefix + tmp_array[i]);
            else text_array.push_back(space + tmp_array[i]);
       }
}

std::string getPDBSymmetry(const std::string& SymOP)
{
       if (SymOP.empty()) return "1555";

       std::string cs = SymOP;
       std::string::size_type p = cs.find("_");
       if (p != std::string::npos) cs.erase(p, 1);
       return cs;
}

void reformatting_charge_to_prefix(std::string& charge)
{
       if (charge.empty() || (charge == "0") || (charge == "?") || (charge == ".")) return;

       int i = atoi(charge.c_str());
       if ((i < 0) || (charge[charge.size() - 1] == '-'))
            charge = "-" + String::IntToString(abs(i));
       else charge = String::IntToString(i);
}

void reformatting_charge_to_suffix(std::string& charge)
{
       if (charge.empty() || (charge == "0") || (charge == "?") || (charge == ".")) return;

       reformatting_charge_to_prefix(charge);
       int i = atoi(charge.c_str());
       if (i > 0)
            charge = String::IntToString(i) + "+";
       else if (i < 0)
            charge = String::IntToString(abs(i)) + "-";
}

void add_missing_zero_before_point(std::string& f)
{
       if (f[0] == '.')
            f.insert(0, "0");
       else if (f.substr(0, 2) == "-.")
            f.insert(1, "0");
}

void reformatting_float_number_without_extra_zero(std::string& f)
{
       add_missing_zero_before_point(f);

       std::vector<std::string> data;
       get_wordarray(data, f, ".");
       if (data.size() != 2) return;

       int pos = -1;
       for (int i = (int) data[1].size() - 1; i >= 0; --i) {
            if (data[1][i] != '0') break;
            pos = i;
       }
       if (pos < 0) return;
       if (pos == 0) pos = 1;
       data[1].erase(pos);
       f = data[0] + "." + data[1];
}

std::set<std::string> get_set_common_element(const std::set<std::string>& set1, const std::set<std::string>& set2)
{
       std::set<std::string> common_set;
       common_set.clear();
       if (set1.empty() || set2.empty()) return common_set;
       for (std::set<std::string>::const_iterator pos = set1.begin(); pos != set1.end(); ++pos) {
            if (set2.find(*pos) == set2.end()) continue;
            common_set.insert(*pos);
       }
       return common_set;
}

std::string get_next_available_asym_id(std::set<std::string>& used_ids)
{
       std::string letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
       return get_next_available_id(letters, used_ids);
}

std::string get_next_available_pdb_chain_id(std::set<std::string>& used_ids)
{
       std::string letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
       return get_next_available_id(letters, used_ids);
}

std::string get_unique_atom_name(std::set<std::string>& existing_name_set, const std::string& name, const std::string& type)
{
       if (existing_name_set.find(name) == existing_name_set.end()) {
            existing_name_set.insert(name);
            return name;
       }
       if (name.size() < 3) {
            std::string letters = "123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
            for (unsigned int i = 0; i < letters.size(); ++i) {
                 std::string new_name = name;
                 new_name += letters[i];
                 if (existing_name_set.find(new_name) == existing_name_set.end()) {
                      existing_name_set.insert(new_name);
                      return new_name;
                 }
            }
       }
       int end = 1000;
       if (type.size() == 2) end = 100;
       for (int i = 1; i < end; ++i) {
            std::string new_name = type + String::IntToString(i);
            if (existing_name_set.find(new_name) == existing_name_set.end()) {
                 existing_name_set.insert(new_name);
                 return new_name;
            }
       }
       return name;
}

bool IsInteger(const std::string& sValue, int& iValue)
{
       iValue = -999999999;
       try {
            iValue = String::StringToInt(sValue);
            return true;
       } catch (const std::exception& e) {
            return false;
       }
}

static bool IsDashValue(const std::string& value)
{
       if (value.empty()) return false;

       bool flag = true;
       for (unsigned int i = 0; i < value.size(); ++i) {
            if (value[i] != '-') {
                 flag = false;
                 break;
            }
       }
       return flag;
}

static std::string get_next_available_id(const std::string& letters, std::set<std::string>& used_ids)
{

       std::string id = "";

       std::vector<std::string> letter_list, id_list;
       letter_list.clear();
       id_list.clear();
       id_list.reserve(used_ids.size() + 3);

       for (unsigned int i = 0; i < letters.size(); ++i) {
            id.clear();
            id += letters[i];
            letter_list.push_back(id);
       }

       id_list.push_back("");
       bool found = false;
       unsigned int start = 0;
       while (true) {
            unsigned int end = id_list.size();
            for (std::vector<std::string>::const_iterator vpos = letter_list.begin(); vpos != letter_list.end(); ++vpos) {
                 for (unsigned int i = start; i < end; ++i) {
                      id = id_list[i] + *vpos;
                      id_list.push_back(id);
                      if (id_list.size() > (used_ids.size() + 2)) {
                           found = true;
                           break;
                      }
                 }
                 if (found) break;
            }
            if (found) break;

            start = end;
       }

       id.clear();
       std::vector<std::string>::const_iterator vpos = id_list.begin();
       vpos++;
       for ( ; vpos != id_list.end(); ++vpos) {
            if (used_ids.find(*vpos) == used_ids.end()) {
                 used_ids.insert(*vpos);
                 id = *vpos;
                 break;
            }
       }
       return id;
}
/*
static std::string get_next_available_id(const std::string& letters, std::set<std::string>& used_ids)
{
       std::string id = "";

       std::vector<std::string> id_list;
       id_list.clear();

       while (true) {
            id = get_next_available_id(letters, used_ids, id_list);
            if (!id.empty()) {
                 used_ids.insert(id);
                 break;
            }
       }

       return id;
}

static std::string get_next_available_id(const std::string& letters, const std::set<std::string>& used_ids, std::vector<std::string>& input_list)
{
       std::vector<std::string> output_list;
       output_list.clear();

       for (unsigned int i = 0; i < letters.size(); ++i) {
            if (input_list.empty()) {
                 std::string id = "";
                 id += letters[i];
                 if (used_ids.find(id) == used_ids.end()) return id;
                 output_list.push_back(id);
            } else {
                 for (std::vector<std::string>::const_iterator pos = input_list.begin(); pos != input_list.end(); ++pos) {
                      std::string id = *pos;
                      id += letters[i];
                      if (used_ids.find(id) == used_ids.end()) return id;
                      output_list.push_back(id);
                 }
            }
       }

       input_list = output_list;
       return "";
}
*/
