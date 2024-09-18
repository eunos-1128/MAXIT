/*
FILE:     time-util.C
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
#include <string.h>
#include <time.h>

#include "utillib.h"

static time_t check_daylight_saving(const time_t local_reading, const time_t final_reading);

void convert_time_format(const int& length, std::string& datestring)
/* ---------------------------------------------------------------
 * converts different format of time output
 * length == 9:  format = day-mon-year, mon is in JAN, FEB,... and year is two digit
 *               (e.g. 09-JAN-13)
 * length == 10: format = year-mon-day, mon are in 1, 2, 3,... and year is four digit
 *               (e.g. 2013-01-09)
 * length == 11: format = day-mon-year, mon is in JAN, FEB,... and year is two digit
 *               (e.g. 09-JAN-2013) 
 * ---------------------------------------------------------------*/
{
       String::StripAndCompressWs(datestring);
       if (datestring.empty()) return;

       std::string::size_type p = datestring.find_first_of(" ;:");
       if (p != std::string::npos) datestring.erase(p);

       if (datestring.size() == 4 || (datestring.size() == 9 && datestring[4] == '-'))
            return;

       std::vector<std::string> data;
       get_wordarray(data, datestring, "-");

       if (data.size() == 1) return;

       std::string day, mon, year;
       day.clear();
       mon.clear();
       year.clear();
       if (data.size() == 2) {
            if (data[0].size() == 4) {
                 mon  = data[1];
                 year = data[0];
            } else {
                 mon  = data[0];
                 year = data[1];
            }
       } else if (data.size() == 3) {
            mon = data[1];
            if (data[0].size() == 4) {
                 day  = data[2];
                 year = data[0];
            } else {
                 day  = data[0];
                 year = data[2];
            }
       }

       std::map<std::string, std::string> mapping;
       mapping.clear();

       if (length == 9 || length == 11) {
            mapping.insert(std::make_pair("1", "JAN"));
            mapping.insert(std::make_pair("01", "JAN"));
            mapping.insert(std::make_pair("2", "FEB"));
            mapping.insert(std::make_pair("02", "FEB"));
            mapping.insert(std::make_pair("3", "MAR"));
            mapping.insert(std::make_pair("03", "MAR"));
            mapping.insert(std::make_pair("4", "APR"));
            mapping.insert(std::make_pair("04", "APR"));
            mapping.insert(std::make_pair("5", "MAY"));
            mapping.insert(std::make_pair("05", "MAY"));
            mapping.insert(std::make_pair("6", "JUN"));
            mapping.insert(std::make_pair("06", "JUN"));
            mapping.insert(std::make_pair("7", "JUL"));
            mapping.insert(std::make_pair("07", "JUL"));
            mapping.insert(std::make_pair("8", "AUG"));
            mapping.insert(std::make_pair("08", "AUG"));
            mapping.insert(std::make_pair("9", "SEP"));
            mapping.insert(std::make_pair("09", "SEP"));
            mapping.insert(std::make_pair("10", "OCT"));
            mapping.insert(std::make_pair("11", "NOV"));
            mapping.insert(std::make_pair("12", "DEC"));
       } else if (length == 10) {
            mapping.insert(std::make_pair("JAN", "01"));
            mapping.insert(std::make_pair("FEB", "02"));
            mapping.insert(std::make_pair("MAR", "03"));
            mapping.insert(std::make_pair("APR", "04"));
            mapping.insert(std::make_pair("MAY", "05"));
            mapping.insert(std::make_pair("JUN", "06"));
            mapping.insert(std::make_pair("JUL", "07"));
            mapping.insert(std::make_pair("AUG", "08"));
            mapping.insert(std::make_pair("SEP", "09"));
            mapping.insert(std::make_pair("OCT", "10"));
            mapping.insert(std::make_pair("NOV", "11"));
            mapping.insert(std::make_pair("DEC", "12"));
       }

       String::UpperCase(mon);
       std::map<std::string, std::string>::const_iterator
           mpos = mapping.find(mon);
       if (mpos != mapping.end()) mon = mpos->second;

       if (length == 9) {
            if (year.size() != 2) {
                 int iyear = atoi(year.c_str());
                 iyear = iyear % 100;
                 year = String::IntToString(iyear);
                 if (year.size() == 1) year = "0" + year;
            }
       } else if (length == 10 || length == 11) {
            if (year.size() != 4) {
                 int iyear = atoi(year.c_str());
                 iyear = iyear % 100;
                 if (iyear > 60)
                      year = "19" + String::IntToString(iyear);
                 else if (iyear < 10)
                      year = "200" + String::IntToString(iyear);
                 else year = "20" + String::IntToString(iyear);
            }
       }

       if (day.size() == 1) day = "0" + day;

       if (length == 9 || length == 11) {
            if (day.empty())
                 datestring.clear();
            else datestring = day + "-";
            datestring += mon + "-" + year;
       } else if (length == 10) {
            datestring = year + "-" + mon;
            if (!day.empty()) datestring += "-" + day;
       }
}

void get_date(std::string &date_string)
{
       time_t start_reading;
       time(&start_reading);
       get_date_from_time(date_string, start_reading);
}

void get_next_tuesday(std::string &next_tuesday)
{
       char stringbuff[20];
       struct tm *tmptr;
       time_t local_reading, final_reading;

       memset(stringbuff, 0, 20);
       time(&local_reading);
       tmptr = localtime(&local_reading);
       strftime(stringbuff, 4, "%w", tmptr);
       int today = atoi(stringbuff);
       switch (today) {
            case 0:
                 final_reading = local_reading + 9 * 86400;
                 break;
            case 1:
                 final_reading = local_reading + 8 * 86400;
                 break;
            case 2:
                 final_reading = local_reading + 7 * 86400;
                 break;
            case 3:
                 final_reading = local_reading + 6 * 86400;
                 break;
            case 4:
                 final_reading = local_reading + 5 * 86400;
                 break;
            case 5:
                 final_reading = local_reading + 4 * 86400;
                 break;
            case 6:
                 final_reading = local_reading + 10 * 86400;
                 break;
            default:
                 final_reading = local_reading + 8 * 86400;
                 break;
       }
       final_reading = check_daylight_saving(local_reading, final_reading);
       get_date_from_time(next_tuesday, final_reading);
}

void get_next_wednesday(std::string &next_wednesday)
{
       char stringbuff[20];
       struct tm *tmptr;
       time_t local_reading, final_reading;

       memset(stringbuff, 0, 20);
       time(&local_reading);
       tmptr = localtime(&local_reading);
       strftime(stringbuff, 4, "%w", tmptr);
       int today = atoi(stringbuff);
       switch (today) {
            case 0:
                 final_reading = local_reading + 10 * 86400;
                 break;
            case 1:
                 final_reading = local_reading + 9 * 86400;
                 break;
            case 2:
                 final_reading = local_reading + 8 * 86400;
                 break;
            case 3:
                 final_reading = local_reading + 7 * 86400;
                 break;
            case 4:
                 final_reading = local_reading + 6 * 86400;
                 break;
            case 5:
                 final_reading = local_reading + 5 * 86400;
                 break;
            case 6:
                 final_reading = local_reading + 11 * 86400;
                 break;
            default:
                 final_reading = local_reading + 9 * 86400;
                 break;
       }
       final_reading = check_daylight_saving(local_reading, final_reading);
       get_date_from_time(next_wednesday, final_reading);
}

void get_next_thursday(std::string &next_thursday)
{
       char stringbuff[20];
       struct tm *tmptr;
       time_t local_reading, final_reading;

       memset(stringbuff, 0, 20);
       time(&local_reading);
       tmptr = localtime(&local_reading);
       strftime(stringbuff, 4, "%w", tmptr);
       int today = atoi(stringbuff);
       switch (today) {
            case 0:
                 final_reading = local_reading + 4 * 86400;
                 break;
            case 1:
                 final_reading = local_reading + 3 * 86400;
                 break;
            case 2:
                 final_reading = local_reading + 2 * 86400;
                 break;
            case 3:
                 final_reading = local_reading + 1 * 86400;
                 break;
            case 4:
                 final_reading = local_reading;
                 break;
            case 5:
                 final_reading = local_reading + 6 * 86400;
                 break;
            case 6:
                 final_reading = local_reading + 5 * 86400;
                 break;
            default:
                 final_reading = local_reading;
                 break;
       }
       final_reading = check_daylight_saving(local_reading, final_reading);
       get_date_from_time(next_thursday, final_reading);
}

void get_next_friday(std::string &next_friday)
{
       char stringbuff[20];
       struct tm *tmptr;
       time_t local_reading, final_reading;

       memset(stringbuff, 0, 20);
       time(&local_reading);
       tmptr = localtime(&local_reading);
       strftime(stringbuff, 4, "%w", tmptr);
       int today = atoi(stringbuff);
       switch (today) {
            case 0:
                 final_reading = local_reading + 5 * 86400;
                 break;
            case 1:
                 final_reading = local_reading + 4 * 86400;
                 break;
            case 2:
                 final_reading = local_reading + 3 * 86400;
                 break;
            case 3:
                 final_reading = local_reading + 2 * 86400;
                 break;
            case 4:
                 final_reading = local_reading + 1 * 86400; 
                 break;
            case 5:
                 final_reading = local_reading;
                 break;
            case 6:
                 final_reading = local_reading + 6 * 86400;
                 break;
            default:
                 final_reading = local_reading;
                 break;
       }
       final_reading = check_daylight_saving(local_reading, final_reading);
       get_date_from_time(next_friday, final_reading);
}

void get_date_from_time(std::string &date_string, const time_t time_reading)
{
       char stringbuff[20];
       struct tm *tmptr = NULL;

       memset(stringbuff, 0, 20);
       tmptr = localtime(&time_reading);
       strftime(stringbuff, 18, "%Y-%m-%d", tmptr);
       date_string = stringbuff;
       String::StripAndCompressWs(date_string);
}

void get_full_date(std::string &full_date, const std::string pdb_date)
{
       int i;
       std::string day, month, year;

       day.clear();
       month.clear();
       year.clear();
       day = pdb_date.substr(0, 2);
       month = pdb_date.substr(3, 3);
       year = pdb_date.substr(7, 2);
       if (month == "JAN")
            month = "01";
       else if (month == "FEB")
            month = "02";
       else if (month == "MAR")
            month = "03";
       else if (month == "APR")
            month = "04";
       else if (month == "MAY")
            month = "05";
       else if (month == "JUN")
            month = "06";
       else if (month == "JUL")
            month = "07";
       else if (month == "AUG")
            month = "08";
       else if (month == "SEP")
            month = "09";
       else if (month == "OCT")
            month = "10";
       else if (month == "NOV")
            month = "11";
       else if (month == "DEC")
            month = "12";
       i = atoi(year.c_str());
       if (i > 60)
            year.insert(0, "19");
       else year.insert(0, "20");
       full_date = year + "-" + month + "-" + day;
}

int is_wrong_date(const std::string& date_string)
{
       char year[20];
       int end_year, i_year, i_month, i_day, is_leap_year = 0, days_per_month;
       struct tm *tmptr;
       time_t local_reading;
       std::vector<std::string> datearray;

       get_wordarray(datearray, date_string, "-");
       if (datearray.size() != 3) return 1;

       memset(year, 0, 20);
       time(&local_reading);
       tmptr = localtime(&local_reading);
       strftime(year, 5, "%Y", tmptr);
       end_year = atoi(year) + 5;

       i_year = atoi(datearray[0].c_str());

       if (i_year < 1970 || i_year > end_year) return 1;

       i_month = atoi(datearray[1].c_str());
       if (i_month < 1 || i_month > 12) return 1;

       i_day = atoi(datearray[2].c_str());

       if ((i_year % 4 == 0 && i_year % 100 != 0) || i_year % 400 == 0)
            is_leap_year = 1;

       switch (i_month) {
             case 1:
             case 3:
             case 5:
             case 7:
             case 8:
             case 10:
             case 12:
                  days_per_month = 31;
                  break;
             case 2:
                  days_per_month = 28;
                  if (is_leap_year) days_per_month = 29;
                  break;
             case 4:
             case 6:
             case 9:
             case 11:
                  days_per_month = 30;
                  break;
             default:
                  days_per_month = 30;
                  break;
       }

       if (i_day < 1 || i_day > days_per_month) return 1;

       return 0;
}

static time_t check_daylight_saving(const time_t local_reading, const time_t final_reading)
{
       struct tm* tmptr = localtime(&local_reading);
       int current_flag = tmptr->tm_isdst;
       tmptr = localtime(&final_reading);
       int future_flag = tmptr->tm_isdst;

       if (current_flag == 0 && future_flag > 0) {
            return (final_reading - 3600);
       } else if (current_flag > 0 && future_flag == 0) {
            return (final_reading + 3600);
       }
       return final_reading;
}
