/*
FILE:     NdbToken.h
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
#ifndef _H_NDB_TOKEN_H_
#define _H_NDB_TOKEN_H_

#include <map>
#include <set>
#include <string>
#include <vector>

#define NDB_DATE_FIELD            7

#include "ISTable.h"
#include "LogUtil.h"

typedef struct {
        std::string FieldLabel;     /* Field name */
        int   FieldType;            /* See above for type descriptions */
        int   FieldWidth;           /* Width of the field */
        int   FieldPrec;            /* Number of digit after decimal point */
        int   FieldStcol;           /* Starting column of the field */
        int   NdbOnly;              /* 1: for NDB used only, 0: used in both
                                     * NDB /PDB */
        int   MixedCase;            /* 1 consists capital and small cases
                                       0 for captial only */
        int   Mandatory;            /* 1: Mandatory yes, 0: Mandatory no */
} field_format;

typedef struct {
        std::string  TokenName;
        int          SeqField;     /* field of serial number */
        int          NumField;     /* number of fields */
        int          ConField;     /* mark for continuation indication */
        int          PdbConField;  /* Pdb continued field, for example
                                    * title consist paper TITLE and descriptor,
                                    * which is two parts for NDB files,
                                    * however, it is one thing for PDB files*/
        int          ContinuedF;   /* fields could be continued */
        int          MatrixType;   /* 0: non-matrix token, 3: matrix token */
        int          NdbOnly;      /* 1: for NDB used only, 0: used in both
                                    * NDB/PDB */
        int          DisplayOrder;
        int          RefineField;
        std::vector<field_format> FieldList;
} ndb_token_format;

typedef struct {
     // int  NdbTokenId;      /* TokenId corresponding to ndb_token_format */
        int  AuxAtom;         /* 1 to identify auxiliary atom, such as ANISOU, etc */
        int  Esd;             /* 1 to identify esd atom, such as SIGUIJ, SIGATM */
        int  SeqAtom;         /* 1 to identify atomserial number increment, TER is
                               * the only one in AuxAtom group with 1 value */
        std::string PdbTitle; /* TokenTitle used in PDB */
     // std::string TableName; /* table name in database */
} atom_format;


class NdbToken
{
   private:
       static void read_token_table(ISTable *Table, std::vector<ndb_token_format> &token_list,
                                    std::vector<std::string>& token_names,
                                    std::set<std::string>& token_name_set);
       static void read_field_table(ISTable *Table, std::vector<ndb_token_format> &token_list,
                                    std::map<std::string, unsigned int>& token_mapping);
       static void read_atom_table(ISTable *Table, std::map<std::string,
                                   atom_format>& _atom_info);
   public:
       static void initialize();
       static bool Read(LogUtil& logIo, const std::string& path);
       static const int& MaxNumField();
       static void resetAtomTokenField(const int& format, const bool& strict_pdb_format_chainid_field_flag = false);
       static bool IsValidToken(const std::string& token);
       static bool IsNdbToken(const std::string& token);
       static bool IsJrnlToken(const std::string& token);
       static bool IsAtomToken(const std::string& token);
       static bool IsMatrixToken(const std::string& token);
       static const ndb_token_format& getTokenFormat(const std::string& token);
       static const std::vector<std::string>& getNdbTokens();
       static const std::vector<std::string>& getJrnlTokens();
};

#endif
