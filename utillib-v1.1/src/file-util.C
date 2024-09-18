/*
FILE:     file-util.C
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
#include <unistd.h>
#include <limits.h>
#include <string.h>
#include <string>
#include <sys/dir.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>

#include "utillib.h"

void get_temp_filename(std::string &temp_filename)
{
       char name[100], templateFileName[9] = "./XXXXXX";

       strcpy(name, templateFileName);
       int fd = mkstemp(name);
       close (fd);

       std::string cs = name;
       std::string::size_type p = cs.find_last_of('/');
       if (p != std::string::npos)
            temp_filename = cs.substr(p + 1);
       else temp_filename = cs;
}

void get_temp_filename(std::string &temp_filename, const std::string& dirPath)
{
       char name[1000];

       strcpy(name, dirPath.c_str());
       if (dirPath[dirPath.size() - 1] == '/')
            strcat(name, "XXXXXX");
       else strcat(name, "/XXXXXX");

       int fd = mkstemp(name);
       close (fd);

       temp_filename = name;
}

void get_temp_directory(std::string &temp_dirname, const std::string& dirPath)
{
       temp_dirname.clear();

       char name[1000];

       strcpy(name, dirPath.c_str());
       if (dirPath[dirPath.size() - 1] == '/')
            strcat(name, "tmpdir.XXXXXX");
       else strcat(name, "/tmpdir.XXXXXX");

       char *dir_name = mkdtemp(name);
       if (dir_name == NULL) return;

       temp_dirname = dir_name;
}

void get_current_directory_path(std::string &dirname)
{
       dirname.clear();
   
       char buff[PATH_MAX + 1];
       char* cwd = getcwd(buff, PATH_MAX + 1);
       if (cwd != NULL) dirname = cwd;
}

void create_directory(const std::string& top_dirname, const std::string& sub_dirname)
{
       std::string full_dirname = top_dirname;
       if (top_dirname[top_dirname.size() - 1] == '/')
            full_dirname += sub_dirname;
       else full_dirname += "/" + sub_dirname;

       remove_directory(full_dirname);
       mkdir(full_dirname.c_str(), 0755);
}

void get_local_file_copy(const std::string &input_filename, std::string &filename_without_path, std::string &filename_of_local_copy)
{
       get_temp_filename(filename_of_local_copy);

       std::string::size_type p = input_filename.find_last_of('/');
       if (p != std::string::npos)
            filename_without_path = input_filename.substr(p + 1);
       else filename_without_path = input_filename;
 
       std::string tmp_name, command;
       int size = filename_without_path.size();
       if (size > 2 && filename_without_path.substr(size - 2, 2) == ".Z") {
            filename_without_path.erase(size - 2);
            get_temp_filename(tmp_name);
            command = "zcat " + input_filename + " > " + tmp_name;
            wsystem(command.c_str(), ".");
            preprocess_file(tmp_name, filename_of_local_copy);
            // command = "rm " + tmp_name;
            // wsystem(command.c_str(), ".");
            remove(tmp_name.c_str());
       } else if (size > 3 && filename_without_path.substr(size - 3, 3) == ".gz") {
            filename_without_path.erase(size - 3);
            get_temp_filename(tmp_name);
            struct utsname sname;
            uname(&sname);
            if (!strncmp(sname.sysname, "Linux", 5))
                 command = "zcat " + input_filename + " > " + tmp_name;
            else command = "gzcat " + input_filename + " > " + tmp_name;
            wsystem(command.c_str(), ".");
            preprocess_file(tmp_name, filename_of_local_copy);
            // command = "rm " + tmp_name;
            // wsystem(command.c_str(), ".");
            remove(tmp_name.c_str());
       } else preprocess_file(input_filename, filename_of_local_copy);
}

void check_path_variable(std::string &path)
{
       if (path.empty()) return;
       if (path[path.size()-1] == '/') return;
       path += '/';
}

int get_line_from_file(FILE *fp, std::string &str)
{
       int  non_ascii = 0;
       long backward = -1;

       str.clear();
       while (!feof(fp)) {
            char c = fgetc(fp);
            if (c == EOF) break;
            else if (c == 9) str += ' ';
            else if (c == '\n' || c == 13) {
                 char i = fgetc(fp);
                 if (i == EOF) break;
                 if ((c == '\n' && i != 13) || (c == 13 && i != '\n')) fseek(fp, backward, SEEK_CUR);
                 break;
            } else if (c > 31 && c < 127) str += c;
            else non_ascii++;
       }
       return non_ascii;
}

void get_line_from_file(FILE *fp, std::string &str, std::vector<std::string>& non_ascii_positions)
{
       long backward = -1;

       str.clear();
       non_ascii_positions.clear();
       while (!feof(fp)) {
            char c = fgetc(fp);
            if (c == EOF) break;
            else if (c == 9) str += ' ';
            else if (c == '\n' || c == 13) {
                 char i = fgetc(fp);
                 if (i == EOF) break;
                 if ((c == '\n' && i != 13) || (c == 13 && i != '\n')) fseek(fp, backward, SEEK_CUR);
                 break;
            } else {
                 str += c;
                 if ((c <= 31) || (c >= 127)) non_ascii_positions.push_back(String::IntToString(str.size()));
            }
       }
}

void get_one_line(FILE *fp, std::string &cs)
{
       char c;

       cs.clear();
       while (!feof(fp)) {
            c = fgetc(fp);
            if (c == '\n' || c == EOF) break;
            cs += c;
       }
}

void get_file_list_from_directory(std::vector<std::string> &file_list, const std::string &dir_name)
{
       file_list.clear();

       DIR *dfd = opendir(dir_name.c_str());
       struct direct *dp = NULL;
       if (dfd) {
            while ((dp = readdir(dfd)) != NULL) {
                 if (strcmp(dp->d_name, ".") && strcmp(dp->d_name, "..")) {
                      std::string cs = dp->d_name;
                      file_list.push_back(cs);
                 }
            }
            closedir(dfd);
       }
}

void remove_directory(const std::string &dir_name)
{
       if (dir_name.empty()) return;

       struct stat statbuf;
       if (stat(dir_name.c_str(), &statbuf) != 0) return;

       if (S_ISDIR(statbuf.st_mode)) {
            DIR *dfd = opendir(dir_name.c_str());
            struct direct *dp = NULL;
            if (dfd) {
                 while ((dp = readdir(dfd)) != NULL) {
                      if (!strcmp(dp->d_name, ".") || !strcmp(dp->d_name, "..")) continue;
                      
                      std::string subdir_name = dir_name;
                      if (dir_name[dir_name.size() - 1] != '/') subdir_name += "/";
                      subdir_name += dp->d_name;
                      remove_directory(subdir_name);
                 }
                 closedir(dfd);
            }
            rmdir(dir_name.c_str());
       } else unlink(dir_name.c_str());
}

void capture_upload_file(const std::string &fileid, const std::string &outfilename, const std::string &outextension, const std::string &dirPath,
                         const std::string &HT_ZCAT_COMMAND, const std::string &HT_GZCAT_COMMAND, const std::string &HT_ASCII_FILTER_COMMAND)
/*
 * outextension should include ".", like ".cif", "-sf.cif"
 */
{
       int i = fileid.size();
       std::string script = "cd " + dirPath + "; ";
       if (i > 2 && fileid[i-2] == '.' && fileid[i-1] == 'Z') {
            script += HT_ZCAT_COMMAND + " " + fileid + " | " + HT_ASCII_FILTER_COMMAND + " > " + outfilename + outextension;
       } else if  (i > 3 && fileid[i-3]=='.' && fileid[i-2]=='g' && fileid[i-1]=='z') {
            script += HT_GZCAT_COMMAND + " " + fileid + " | " + HT_ASCII_FILTER_COMMAND + " > " + outfilename + outextension;
       } else {
            std::string distPath = dirPath; check_path_variable(distPath);
            std::string orgiPath = "";
            std::string orgiFile = "";
            std::string distFile = outfilename + outextension;
            std::string::size_type pos = fileid.find_last_of('/');
            if (pos != std::string::npos) {
                 orgiPath = fileid.substr(0, pos + 1);
                 if (orgiPath.size() < fileid.size()) {
                      orgiFile = fileid.substr(orgiPath.size());
                 }
            } else orgiFile = fileid;

            if (distPath != orgiPath || distFile != orgiFile) {
                 script += HT_ASCII_FILTER_COMMAND + " < " + fileid + " > " + outfilename + outextension;
            } else {
                 distFile += ".org";
                 std::string filename = fileid + ".org";
                 rename(fileid.c_str(), filename.c_str());
                 script += HT_ASCII_FILTER_COMMAND + " < " + filename + " > " + outfilename + outextension;
            }
       }
       wsystem(script.c_str(), dirPath.c_str());
}

int display_raw_file(FILE *io, const std::string &fName)
{
       int iret = 0;
       if (fName.empty()) return iret;

       FILE *fp = fopen(fName.c_str(), "r");
       if (fp != NULL) {
            struct stat statbuf;
            stat(fName.c_str(), &statbuf);
            char *data = new char[statbuf.st_size+1];
            memset(data, 0, statbuf.st_size+1);
            size_t size = fread(data, sizeof(char), statbuf.st_size+1, fp);
            fclose (fp);

            if ((data != NULL) && (size > 0)) {
                 fprintf(io, "%s", data);
                 delete [] data;
                 iret = 1;
            }
       }
       return(iret);
}

void copy_file(const std::string &old_file, const std::string &new_file)
{
       FILE *fp = fopen(new_file.c_str(), "w");
       int ret = display_raw_file(fp, old_file);
       fclose (fp);
       if (!ret) remove(new_file.c_str());
}

void get_file_content(const std::string &fName, std::string &content)
{
       content.clear();
       if (fName.empty()) return;

       FILE *fp = fopen(fName.c_str(), "r");
       if (fp != NULL) {
            struct stat statbuf;
            stat(fName.c_str(), &statbuf);
            char *data = new char[statbuf.st_size+1];
            memset(data, 0, statbuf.st_size+1);
            size_t size = fread(data, sizeof(char), statbuf.st_size+1, fp);
            fclose (fp);

            if ((data != NULL) && (size > 0)) {
                 content = data;
                 delete [] data;
            }
       }
}

void preprocess_file(const std::string& infile, const std::string& outfile)
{
       FILE *in = fopen(infile.c_str(), "r");
       FILE *out = fopen(outfile.c_str(), "w");
       int c = getc(in);
       while (!feof(in)) {
            int i = getc(in);
            if (c == 13 && i != 10) c = 10;
            if (c > 31 || c == 10)
                 putc(c, out);
            else if (c == 9 || c > 126)
                 putc(' ', out);
            c = i;
       }
       fclose (in);
       fclose (out);
}

void get_address_book(const std::string& annotator_e_mail_file, std::map<std::string, std::vector<std::string> >& address_book)
{
       address_book.clear();

       std::vector<std::string> array;
       std::string cs, initial, address, name;

       FILE *fp = fopen(annotator_e_mail_file.c_str(), "r");
       while (!feof(fp)) {
            get_one_line(fp, cs);
            if (feof(fp)) break;
            get_wordarray(array, cs, " \n");
            String::UpperCase(array[0], initial);
            address = array[1];
            name.clear();
            for (unsigned int i = 2; i < array.size(); i++) {
                 if (!name.empty()) name += " ";
                 name += array[i];
            }
            array.clear();
            array.push_back(address);
            array.push_back(name);
            address_book.insert(std::make_pair(initial, array));
       }
       fclose (fp);
}
