/*
FILE:     SwitchLabelingUtil_global.h
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
#ifndef _H_SWITCH_LABELING_UTIL_GLOBAL_H_
#define _H_SWITCH_LABELING_UTIL_GLOBAL_H_

#include <map>
#include <string>

#include "SwitchLabelingUtil.h"

static CHIRAL_LABELINGS chiral_labeling_na[3] = {
       { 2, 0, {{"C2'", "C1'", "C3'", "H2'", "H2''",1,  1.0, 0},
                {"C2'", "C1'", "C3'", "H2''","H2'",1, -1.0, 0}} },
       { 2, 0, {{"C5'", "O5'", "C4'", "H5'", "H5''",1,  1.0, 0},
                {"C5'", "O5'", "C4'", "H5''","H5'",1, -1.0, 0}} },
       { 2, 1, {{"P",   "O5'", "O3'", "OP2", "OP1", 0,  1.0, 0},
                {"P",   "O5'", "O3'", "OP1", "OP2", 0, -1.0, 0}} }
};

static CHIRAL_LABELINGS chiral_labeling_ARG[3] = {
       { 2, 0, {{"CB", "CA", "CG",  "HB2", "HB3", 1, -1.0, 0},
                {"CB", "CA", "CG",  "HB3", "HB2", 1,  1.0, 0}} },
       { 2, 0, {{"CG", "CB", "CD",  "HG2", "HG3", 1, -1.0, 0},
                {"CG", "CB", "CD",  "HG3", "HG2", 1,  1.0, 0}} },
       { 2, 0, {{"CD", "CG", "NE",  "HD2", "HD3", 1, -1.0, 0},
                {"CD", "CG", "NE",  "HD3", "HD2", 1,  1.0, 0}} }
};

static CHIRAL_LABELINGS chiral_labeling_ASN = {
         2, 0, {{"CB", "CA", "CG",  "HB2", "HB3", 1, -1.0, 0},
                {"CB", "CA", "CG",  "HB3", "HB2", 1,  1.0, 0}}
};

static CHIRAL_LABELINGS chiral_labeling_CYS = {
         2, 0, {{"CB", "CA", "SG",  "HB2", "HB3", 1, -1.0, 0},
                {"CB", "CA", "SG",  "HB3", "HB2", 1,  1.0, 0}}
};

static CHIRAL_LABELINGS chiral_labeling_GLN[2] = {
       { 2, 0, {{"CB", "CA", "CG",  "HB2", "HB3", 1, -1.0, 0},
                {"CB", "CA", "CG",  "HB3", "HB2", 1,  1.0, 0}} },
       { 2, 0, {{"CG", "CB", "CD",  "HG2", "HG3", 1, -1.0, 0},
                {"CG", "CB", "CD",  "HG3", "HG2", 1,  1.0, 0}} }
};

static CHIRAL_LABELINGS chiral_labeling_GLY = {
         2, 0, {{"CA", "N",  "C",   "HA2", "HA3", 1, -1.0, 0},
                {"CA", "N",  "C",   "HA3", "HA2", 1,  1.0, 0}}
};

static CHIRAL_LABELINGS chiral_labeling_ILE = {
         2, 0, {{"CG1","CB", "CD1", "HG12","HG13",1, -1.0, 0},
                {"CG1","CB", "CD1", "HG13","HG12",1,  1.0, 0}}
};

static CHIRAL_LABELINGS chiral_labeling_LEU[2] = {
       { 2, 0, {{"CB", "CA", "CG",  "HB2", "HB3", 1, -1.0, 0},
                {"CB", "CA", "CG",  "HB3", "HB2", 1,  1.0, 0}} },
       { 2, 0, {{"CG", "CB", "CD1", "CD2", "CD1", 0, -1.0, 3,
               {{"HD21","HD11"},{"HD22","HD12"},{"HD23","HD13"}}}, 
                {"CG", "CB", "CD2", "CD1", "CD2", 0,  1.0, 3,
               {{"HD11","HD21"},{"HD12","HD22"},{"HD13","HD23"}}}}}
};

static CHIRAL_LABELINGS chiral_labeling_LYS[4] = {
       { 2, 0, {{"CB", "CA", "CG",  "HB2", "HB3", 1, -1.0, 0},
                {"CB", "CA", "CG",  "HB3", "HB2", 1,  1.0, 0}} },
       { 2, 0, {{"CG", "CB", "CD",  "HG2", "HG3", 1, -1.0, 0},
                {"CG", "CB", "CD",  "HG3", "HG2", 1,  1.0, 0}} },
       { 2, 0, {{"CD", "CG", "CE",  "HD2", "HD3", 1, -1.0, 0},
                {"CD", "CG", "CE",  "HD3", "HD2", 1,  1.0, 0}} },
       { 2, 0, {{"CE", "CD", "NZ",  "HE2", "HE3", 1, -1.0, 0},
                {"CE", "CD", "NZ",  "HE3", "HE2", 1,  1.0, 0}} }
};

static CHIRAL_LABELINGS chiral_labeling_MET[2] = {
       { 2, 0, {{"CB", "CA", "CG",  "HB2", "HB3", 1, -1.0, 0},
                {"CB", "CA", "CG",  "HB3", "HB2", 1,  1.0, 0}} },
       { 2, 0, {{"CG", "CB", "SD",  "HG2", "HG3", 1, -1.0, 0},
                {"CG", "CB", "SD",  "HG3", "HG2", 1,  1.0, 0}} }
};

static CHIRAL_LABELINGS chiral_labeling_PRO[4] = {
       { 2, 0, {{"N",  "CA", "CD",  "HT1", "HT2", 1, -1.0, 0},
                {"N",  "CA", "CD",  "HT2", "HT1", 1,  1.0, 0}} },
       { 2, 0, {{"CB", "CA", "CG",  "HB2", "HB3", 1, -1.0, 0},
                {"CB", "CA", "CG",  "HB3", "HB2", 1,  1.0, 0}} },
       { 2, 0, {{"CG", "CB", "CD",  "HG2", "HG3", 1, -1.0, 0},
                {"CG", "CB", "CD",  "HG3", "HG2", 1,  1.0, 0}} },
       { 2, 0, {{"CD", "CG", "N",   "HD2", "HD3", 1, -1.0, 0},
                {"CD", "CG", "N",   "HD3", "HD2", 1,  1.0, 0}} }
};

static CHIRAL_LABELINGS chiral_labeling_SER = {
         2, 0, {{"CB", "CA", "OG",  "HB2", "HB3", 1, -1.0, 0},
                {"CB", "CA", "OG",  "HB3", "HB2", 1,  1.0, 0}}
};

static CHIRAL_LABELINGS chiral_labeling_VAL = {
         2, 0, {{"CB", "CA", "CG1", "CG2", "CG1", 0, -1.0, 3,
               {{"HG21","HG11"},{"HG22","HG12"},{"HG23","HG13"}}},
                {"CB", "CA", "CG2", "CG1", "CG2", 0,  1.0, 3,
               {{"HG11","HG21"},{"HG12","HG22"},{"HG13","HG23"}}}} 
};

#define NUM_OF_RESIDUES  30

static CHIRAL_LABELING_POOL total_list[NUM_OF_RESIDUES] = {
       { "A",    3, chiral_labeling_na   },
       { "ARG",  3, chiral_labeling_ARG  },
       { "ASN",  1, &chiral_labeling_ASN },
       { "ASP",  1, &chiral_labeling_ASN },
       { "C",    3, chiral_labeling_na   },
       { "CYS",  1, &chiral_labeling_CYS },
       { "DA",   3, chiral_labeling_na   },
       { "DC",   3, chiral_labeling_na   },
       { "DG",   3, chiral_labeling_na   },
       { "DI",   3, chiral_labeling_na   },
       { "DT",   3, chiral_labeling_na   },
       { "DU",   3, chiral_labeling_na   },
       { "G",    3, chiral_labeling_na   },
       { "GLN",  2, chiral_labeling_GLN  },
       { "GLU",  2, chiral_labeling_GLN  },
       { "GLY",  1, &chiral_labeling_GLY },
       { "HIS",  1, &chiral_labeling_ASN },
       { "I",    3, chiral_labeling_na   },
       { "ILE",  1, &chiral_labeling_ILE },
       { "LEU",  2, chiral_labeling_LEU  },
       { "LYS",  4, chiral_labeling_LYS  },
       { "MET",  2, chiral_labeling_MET  },
       { "PHE",  1, &chiral_labeling_ASN },
       { "PRO",  4, chiral_labeling_PRO  },
       { "SER",  1, &chiral_labeling_SER },
       { "T",    3, chiral_labeling_na   },
       { "TRP",  1, &chiral_labeling_ASN },
       { "TYR",  1, &chiral_labeling_ASN },
       { "U",    3, chiral_labeling_na   },
       { "VAL",  1, &chiral_labeling_VAL }
};

static std::map<std::string, int> pool_index;

#endif
