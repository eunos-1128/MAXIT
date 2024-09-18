/*
FILE:     MessageUtil.C
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
#include <stdarg.h>

#include "MessageUtil.h"
#include "utillib.h"

MessageUtil::MessageUtil()
{
       clear();
}

MessageUtil::~MessageUtil()
{
       clear();
}

void MessageUtil::clear()
{
       _formatted_length = 80;
       _messageMap.clear();
}

void MessageUtil::clear_by_type(const std::string& type)
{
       if (_messageMap.find(type) != _messageMap.end()) _messageMap.erase(type);
}

bool MessageUtil::empty() { return _messageMap.empty(); }

bool MessageUtil::empty(const std::string& type) { return (_messageMap.find(type) == _messageMap.end()); }

bool MessageUtil::empty(const std::string& type, const std::string& sub_type)
{
       std::map<std::string, std::map<std::string, std::vector<std::map<std::string, std::string> > > >::const_iterator pos = _messageMap.find(type);
       if (pos == _messageMap.end()) return true;

       return (pos->second.find(sub_type) == pos->second.end());
}

void MessageUtil::set_formatted_length(const unsigned int& default_length)
{
       _formatted_length = default_length;
}

void MessageUtil::insertMessage(const std::string& type, const std::string& sub_type, const std::string& message, const bool& reformat_flag,
                                const std::string& line, const bool& column_flag, const std::string& msg_type)
{
       std::map<std::string, std::string> msg_map;
       msg_map.clear();
       msg_map["message"] = message;
       if (!line.empty()) msg_map["line"] = line;
       if (column_flag) msg_map["column"] = "yes";
       if (reformat_flag) msg_map["reformat"] = "yes";
       if (!msg_type.empty()) msg_map["type"] = msg_type;

       std::vector<std::map<std::string, std::string> > t_data;
       t_data.clear();
       t_data.push_back(msg_map);

       std::map<std::string, std::vector<std::map<std::string, std::string> > > t_map;
       t_map.clear();
       t_map.insert(std::make_pair(sub_type, t_data));

       std::map<std::string, std::map<std::string, std::vector<std::map<std::string, std::string> > > >::iterator pos = _messageMap.find(type);
       if (pos == _messageMap.end()) {
            _messageMap.insert(std::make_pair(type, t_map));
       } else {
            std::map<std::string, std::vector<std::map<std::string, std::string> > >::iterator pos1 = pos->second.find(sub_type);
            if (pos1 != pos->second.end()) {
                 if (sub_type == "label") return;
                 pos1->second.push_back(msg_map);
            } else pos->second.insert(std::make_pair(sub_type, t_data));
       }
}

void MessageUtil::getAllMessage(std::vector<std::string>& messages)
{
       messages.clear();
       if (_messageMap.empty()) return;

       for (std::map<std::string, std::map<std::string, std::vector<std::map<std::string, std::string> > > >::const_iterator
            pos = _messageMap.begin(); pos != _messageMap.end(); ++pos) {
            for (std::map<std::string, std::vector<std::map<std::string, std::string> > >::const_iterator
                 pos1 = pos->second.begin(); pos1 != pos->second.end(); ++pos1) {
                 _getMessage(messages, pos1->second);
            }
       }
}

void MessageUtil::getMessageByType(std::vector<std::string>& messages, const std::string& type)
{
       messages.clear();

       std::map<std::string, std::map<std::string, std::vector<std::map<std::string, std::string> > > >::const_iterator pos = _messageMap.find(type);
       if (pos == _messageMap.end()) return;

       for (std::map<std::string, std::vector<std::map<std::string, std::string> > >::const_iterator
            pos1 = pos->second.begin(); pos1 != pos->second.end(); ++pos1) {
            _getMessage(messages, pos1->second);
       }
}

void MessageUtil::getMessageByType(std::vector<std::string>& messages, const std::string& type, const std::string& sub_type)
{
       messages.clear();

       std::map<std::string, std::map<std::string, std::vector<std::map<std::string, std::string> > > >::const_iterator pos = _messageMap.find(type);
       if (pos == _messageMap.end()) return;

       std::map<std::string, std::vector<std::map<std::string, std::string> > >::const_iterator pos1 = pos->second.find(sub_type);
       if (pos1 == pos->second.end()) return;

       _getMessage(messages, pos1->second);
}

std::string MessageUtil::getMessageByType(const std::string& type, const std::string& sub_type)
{
       std::map<std::string, std::map<std::string, std::vector<std::map<std::string, std::string> > > >::const_iterator mpos = _messageMap.find(type);
       if (mpos == _messageMap.end()) return "";

       std::map<std::string, std::vector<std::map<std::string, std::string> > >::const_iterator mpos1 = mpos->second.find(sub_type);
       if (mpos1 == mpos->second.end()) return "";

       std::string msg = "";
       std::vector<std::string> msg_data;
       for (std::vector<std::map<std::string, std::string> >::const_iterator pos = mpos1->second.begin(); pos != mpos1->second.end(); ++pos) {
            std::map<std::string, std::string>::const_iterator mpos2 = pos->find("message");
            if (mpos2 == pos->end()) continue;

            bool reformat_flag = false;
            std::map<std::string, std::string>::const_iterator mpos3 = pos->find("reformat");
            if (mpos3 != pos->end()) reformat_flag = true;

            if (!reformat_flag) msg += mpos2->second + "\n\n";
            else {
                 _reformatString(false, mpos2->second, msg_data);
                 for (std::vector<std::string>::const_iterator vpos = msg_data.begin(); vpos != msg_data.end(); ++vpos) {
                      msg += *vpos + "\n";
                 }
                 msg += "\n";
            }
       }
       return msg;
}

void MessageUtil::getAllMessage(std::vector<std::pair<std::string, std::string> >& messages)
{
       messages.clear();
       if (_messageMap.empty()) return;

       for (std::map<std::string, std::map<std::string, std::vector<std::map<std::string, std::string> > > >::const_iterator
            pos = _messageMap.begin(); pos != _messageMap.end(); ++pos) {
            for (std::map<std::string, std::vector<std::map<std::string, std::string> > >::const_iterator
                 pos1 = pos->second.begin(); pos1 != pos->second.end(); ++pos1) {
                 _getMessage(messages, pos1->second);
            }
       }
}

void MessageUtil::getMessageByType(std::vector<std::pair<std::string, std::string> >& messages, const std::string& type)
{
       messages.clear();

       std::map<std::string, std::map<std::string, std::vector<std::map<std::string, std::string> > > >::const_iterator pos = _messageMap.find(type);
       if (pos == _messageMap.end()) return;

       for (std::map<std::string, std::vector<std::map<std::string, std::string> > >::const_iterator
            pos1 = pos->second.begin(); pos1 != pos->second.end(); ++pos1) {
            _getMessage(messages, pos1->second);
       }
}

void MessageUtil::getMessageByType(std::vector<std::pair<std::string, std::string> >& messages, const std::string& type, const std::string& sub_type)
{
       messages.clear();

       std::map<std::string, std::map<std::string, std::vector<std::map<std::string, std::string> > > >::const_iterator pos = _messageMap.find(type);
       if (pos == _messageMap.end()) return;

       std::map<std::string, std::vector<std::map<std::string, std::string> > >::const_iterator pos1 = pos->second.find(sub_type);
       if (pos1 == pos->second.end()) return;

       _getMessage(messages, pos1->second);
}

void MessageUtil::printAllMessage(FILE* fp)
{
       if (_messageMap.empty()) return;

       for (std::map<std::string, std::map<std::string, std::vector<std::map<std::string, std::string> > > >::const_iterator
            pos = _messageMap.begin(); pos != _messageMap.end(); ++pos) {
            for (std::map<std::string, std::vector<std::map<std::string, std::string> > >::const_iterator
                 pos1 = pos->second.begin(); pos1 != pos->second.end(); ++pos1) {
                 _printMessage(fp, pos->first, pos1->second);
            }
       }
}

void MessageUtil::printMessageByType(FILE* fp, const std::string& type, const std::map<std::string, std::string>& header)
{
       std::map<std::string, std::map<std::string, std::vector<std::map<std::string, std::string> > > >::const_iterator pos = _messageMap.find(type);
       if (pos == _messageMap.end()) return;

       std::map<std::string, std::string>::const_iterator mpos = header.find(type);
       if (mpos != header.end()) fprintf(fp, "%s\n", mpos->second.c_str());

       for (std::map<std::string, std::vector<std::map<std::string, std::string> > >::const_iterator
            pos1 = pos->second.begin(); pos1 != pos->second.end(); ++pos1) {
            _printMessage(fp, type, pos1->second);
       }
}

void MessageUtil::printMessageByType(FILE* fp, const std::string& type, const std::string& sub_type)
{
       std::map<std::string, std::map<std::string, std::vector<std::map<std::string, std::string> > > >::const_iterator pos = _messageMap.find(type);
       if (pos == _messageMap.end()) return;

       std::map<std::string, std::vector<std::map<std::string, std::string> > >::const_iterator pos1 = pos->second.find(sub_type);
       if (pos1 == pos->second.end()) return;

       _printMessage(fp, type, pos1->second);
}

void MessageUtil::printAllMessage(LogUtil* log)
{
       if (_messageMap.empty()) return;

       for (std::map<std::string, std::map<std::string, std::vector<std::map<std::string, std::string> > > >::const_iterator
            pos = _messageMap.begin(); pos != _messageMap.end(); ++pos) {
            for (std::map<std::string, std::vector<std::map<std::string, std::string> > >::const_iterator
                 pos1 = pos->second.begin(); pos1 != pos->second.end(); ++pos1) {
                 _printMessage(log, pos->first, pos1->second);
            }
       }
}

void MessageUtil::printMessageByType(LogUtil* log, const std::string& type, const std::map<std::string, std::string>& header)
{
       std::map<std::string, std::map<std::string, std::vector<std::map<std::string, std::string> > > >::const_iterator pos = _messageMap.find(type);
       if (pos == _messageMap.end()) return;

       std::map<std::string, std::string>::const_iterator mpos = header.find(type);
       if (mpos != header.end()) {
            log->message(mpos->second.c_str());
            log->message("\n");
       }

       for (std::map<std::string, std::vector<std::map<std::string, std::string> > >::const_iterator
            pos1 = pos->second.begin(); pos1 != pos->second.end(); ++pos1) {
            _printMessage(log, type, pos1->second);
       }
}

void MessageUtil::printMessageByType(LogUtil* log, const std::string& type, const std::string& sub_type)
{
       std::map<std::string, std::map<std::string, std::vector<std::map<std::string, std::string> > > >::const_iterator pos = _messageMap.find(type);
       if (pos == _messageMap.end()) return;

       std::map<std::string, std::vector<std::map<std::string, std::string> > >::const_iterator pos1 = pos->second.find(sub_type);
       if (pos1 == pos->second.end()) return;

       _printMessage(log, type, pos1->second);
}

void MessageUtil::_getMessage(std::vector<std::string>& messages, const std::vector<std::map<std::string, std::string> >& msgMap)
{
       for (std::vector<std::map<std::string, std::string> >::const_iterator pos = msgMap.begin(); pos != msgMap.end(); ++pos) {
            std::map<std::string, std::string>::const_iterator mpos = pos->find("message");
            if (mpos != pos->end()) messages.push_back(mpos->second);
       }
}

void MessageUtil::_getMessage(std::vector<std::pair<std::string, std::string> >& messages, const std::vector<std::map<std::string, std::string> >& msgMap)
{
       for (std::vector<std::map<std::string, std::string> >::const_iterator pos = msgMap.begin(); pos != msgMap.end(); ++pos) {
            std::string type = "";
            std::map<std::string, std::string>::const_iterator mpos = pos->find("type");
            if (mpos != pos->end()) type = mpos->second;
            mpos = pos->find("message");
            if (mpos != pos->end()) messages.push_back(std::make_pair(type, mpos->second));
       }
}

void MessageUtil::_printMessage(FILE* fp, const std::string& type, const std::vector<std::map<std::string, std::string> >& msgMap)
{
       std::string msg, line, upperType;
       bool column_flag, reformat_flag; 
       std::vector<std::string> msg_data;

       bool error_flag = false;
       String::UpperCase(type, upperType);
       if (upperType == "ERROR") error_flag = true;
 
       for (std::vector<std::map<std::string, std::string> >::const_iterator pos = msgMap.begin(); pos != msgMap.end(); ++pos) {
            _getMessageInfo(*pos, msg, line, column_flag, reformat_flag);
            if (msg.empty()) continue;

            if (!reformat_flag) {
                 if (error_flag)
                      fprintf(fp, "ERROR: %s\n\n", msg.c_str());
                 else fprintf(fp, "%s\n\n", msg.c_str());
            } else {
                 _reformatString(error_flag, msg, msg_data);
                 for (std::vector<std::string>::const_iterator vpos = msg_data.begin(); vpos != msg_data.end(); ++vpos) {
                      fprintf(fp, "%s\n", vpos->c_str());
                 }
                 fprintf(fp, "\n");
            }
            if (line.empty()) continue;

            if (column_flag) fprintf(fp, "12345678901234567890123456789012345678901234567890123456789012345678901234567890\n");
            fprintf(fp, "%s\n\n", line.c_str());
       }
}

void MessageUtil::_printMessage(LogUtil* log, const std::string& type, const std::vector<std::map<std::string, std::string> >& msgMap)
{
       std::string msg, line, upperType;
       bool column_flag, reformat_flag; 
       std::vector<std::string> msg_data;

       bool error_flag = false;
       String::UpperCase(type, upperType);
       if (upperType == "ERROR") error_flag = true;
 
       for (std::vector<std::map<std::string, std::string> >::const_iterator pos = msgMap.begin(); pos != msgMap.end(); ++pos) {
            _getMessageInfo(*pos, msg, line, column_flag, reformat_flag);
            if (msg.empty()) continue;

            if (!reformat_flag) {
                 if (error_flag) log->message("ERROR: ");
                 log->message(msg.c_str());
                 log->message("\n\n");
            } else {
                 _reformatString(error_flag, msg, msg_data);
                 for (std::vector<std::string>::const_iterator vpos = msg_data.begin(); vpos != msg_data.end(); ++vpos) {
                      log->message(vpos->c_str());
                      log->message("\n");
                 }
                 log->message("\n");
            }
            if (line.empty()) continue;

            if (column_flag) log->message("12345678901234567890123456789012345678901234567890123456789012345678901234567890\n");
            log->message(line.c_str());
            log->message("\n\n");
       }
}

void MessageUtil::_getMessageInfo(const std::map<std::string, std::string>& msgMap, std::string& msg, std::string& line,
                                  bool& column_flag, bool& reformat_flag)
{
       msg.clear();
       line.clear();
       column_flag = false;
       reformat_flag = false;

       std::map<std::string, std::string>::const_iterator mpos = msgMap.find("message"); if (mpos != msgMap.end()) msg = mpos->second; 
       mpos = msgMap.find("line"); if (mpos != msgMap.end()) line = mpos->second;
       mpos = msgMap.find("column"); if (mpos != msgMap.end()) column_flag = true;
       mpos = msgMap.find("reformat"); if (mpos != msgMap.end()) reformat_flag = true;
}

void MessageUtil::_reformatString(const bool& flag, const std::string& msg, std::vector<std::string>& msg_data)
{
       msg_data.clear();

       std::string line = "";
       if (flag) line = "ERROR:";

       std::vector<std::string> data;
       get_wordarray(data, msg, " ");
       unsigned int len = 0;
       for (std::vector<std::string>::const_iterator pos = data.begin(); pos != data.end(); ++pos) {
            if ((len + pos->size() + 1) > _formatted_length) {
                 msg_data.push_back(line);
                 len = 0;
                 line = "";
                 if (flag) line = "      ";
            }
            if (!line.empty()) {
                 line += " ";
                 len += 1;
            }
            line += *pos;
            len += pos->size();
       }
       if (len) msg_data.push_back(line);
}
