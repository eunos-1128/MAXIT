/*
FILE:     Remark_3_Phenix.h
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

REMARKS Remark_3_Phenix1[Num_Remark_3_Phenix1] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DATA USED IN REFINEMENT."}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"RESOLUTION RANGE HIGH (ANGSTROMS) :" },
    { "RFACTR",    7, 39, 2, 6, 2, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"RESOLUTION RANGE LOW  (ANGSTROMS) :" },
    { "RFACTR",    6, 39, 2, 7, 2, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,{
    { "",              0,  3, 3, 0, 0,PDB_REMARK_LEFT_JUSTIFIED,
        "MIN(FOBS/SIGMA_FOBS)              :"},
    {"RFACTR",   5,  39, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
       "Sigma_F_Ref"},}
  },
  { 3, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"COMPLETENESS FOR RANGE        (%) :"},
    { "PEROBS",  2,  39, 2, 6, 1, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"NUMBER OF REFLECTIONS             :"},
    { "RFACTR",  3,   39, 1, 10, 0,  PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
};

REMARKS Remark_3_Phenix2[Num_Remark_3_Phenix2] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FIT TO DATA USED IN REFINEMENT."}}
  },
  { 3,  2, PDB_REMARK_R_VALUE, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"R VALUE     (WORKING + TEST SET) :"},
    { "RFACTR",  8,   38, 2, 6, 3,  PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_R_VALUE, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"R VALUE            (WORKING SET) :"},
    { "RFACTR",  11,   38, 2, 6, 3,  PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_R_VALUE, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FREE R VALUE                     :"},
    { "FREERF",  2,   38, 2, 5, 3,  PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FREE R VALUE TEST SET SIZE   (%) :"},
    { "FREERF",  4,   38, 2, 6, 3,  PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FREE R VALUE TEST SET COUNT      :"},
    { "FREERF",  5,   38, 1, 5, 0,  PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_Phenix_BIN[Num_Remark_3_Phenix_BIN] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FIT TO DATA USED IN REFINEMENT (IN BINS)."}}
  },
  { 3,  1, PDB_REMARK_R_VALUE, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BIN  RESOLUTION RANGE  COMPL.    NWORK NFREE   RWORK  RFREE"}}
  }
};

REMARKS Remark_3_PHENIX_BULK[Num_Remark_3_PHENIX_BULK] = {
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BULK SOLVENT MODELLING."}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "METHOD USED        :"},
    { "SOLMOD", 2, 24, 3, 36, 3, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SOLVENT RADIUS     :"},
    { "SOLMOD", 5, 24, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SHRINKAGE RADIUS   :"},
    { "SOLMOD", 7, 24, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "K_SOL              :"},
    { "SOLMOD", 3, 24, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "B_SOL              :"},
    { "SOLMOD", 4, 24, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  }
};

REMARKS Remark_3_PHENIX_ESU[Num_Remark_3_PHENIX_ESU] = {
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "ERROR ESTIMATES."}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "COORDINATE ERROR (MAXIMUM-LIKELIHOOD BASED)     :"},
    { "REFESU", 4, 53, 2, 7, 3, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "PHASE ERROR (DEGREES, MAXIMUM-LIKELIHOOD BASED) :"},
    { "REFESU", 6, 53, 2, 7, 3, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  }
};

REMARKS Remark_3_PHENIX_RMSD[Num_Remark_3_PHENIX_RMSD] = {
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "DEVIATIONS FROM IDEAL VALUES."}}
  },
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 17, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "RMSD          COUNT"}}
  },
  { 3, 3, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BOND      :"},
    { "RMSTN1", 2, 15, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "RMSTN1", 4, 30, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 3, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "ANGLE     :"},
    { "RMSTN1", 5, 15, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "RMSTN1", 7, 30, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 3, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "CHIRALITY :"},
    { "RMSTN1",11, 15, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "RMSTN1",13, 30, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 3, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "PLANARITY :"},
    { "RMSTN2", 5, 15, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "RMSTN2", 7, 30, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 3, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "DIHEDRAL  :"},
    { "RMSTN1", 8, 15, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "RMSTN1",10, 30, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  }
};

REMARKS Remark_3_BValue_PHENIX[Num_Remark_3_BValue_PHENIX] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"B VALUES."}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FROM WILSON PLOT           (A**2) :" },
    { "BVALUE",    3, 39, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"MEAN B VALUE      (OVERALL, A**2) :" },
    { "BVALUE",    2, 39, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"OVERALL ANISOTROPIC B VALUE."}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"B11 (A**2) :" },
    { "BVALUE", 4, 17, 2, 10, 5, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"B11 :" },
    { "BVALUE", 4, 10, 2, 10, 5, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"B22 (A**2) :" },
    { "BVALUE", 5, 17, 2, 10, 5, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"B22 :" },
    { "BVALUE", 5, 10, 2, 10, 5, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"B33 (A**2) :" },
    { "BVALUE", 6, 17, 2, 10, 5, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"B33 :" },
    { "BVALUE", 6, 10, 2, 10, 5, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"B12 (A**2) :" },
    { "BVALUE", 7, 17, 2, 10, 5, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"B12 :" },
    { "BVALUE", 7, 10, 2, 10, 5, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"B13 (A**2) :" },
    { "BVALUE", 8, 17, 2, 10, 5, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"B13 :" },
    { "BVALUE", 8, 10, 2, 10, 5, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"B23 (A**2) :" },
    { "BVALUE", 9, 17, 2, 10, 5, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"B23 :" },
    { "BVALUE", 9, 10, 2, 10, 5, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_PHENIX_TWINNING[Num_Remark_3_PHENIX_TWINNING] = {
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TWINNING INFORMATION."}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "FRACTION:"},
    { "TWIN", 3, 13, 2, 6, 4, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "OPERATOR:"},
    { "TWIN", 2, 13, 3, 47, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  }
};

REMARKS Remark_3_PHENIX_TLS_GROUP[Num_Remark_3_PHENIX_TLS_GROUP] = {
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TLS GROUP :"},
    { "TLSGRO", 2, 15, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SELECTION:"},
    { "TLSGRO", 8, 15, 3, 45, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "ORIGIN FOR THE GROUP (A):"},
    { "TLSGRO", 4, 29, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "TLSGRO", 5, 38, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "TLSGRO", 6, 47, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "T TENSOR"}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "", 0, 6, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "T11:"},
    { "TTENSR", 3, 10, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "T22:"},
    { "TTENSR", 4, 24, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "", 0, 6, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "T33:"},
    { "TTENSR", 5, 10, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "T12:"},
    { "TTENSR", 6, 24, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "", 0, 6, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "T13:"},
    { "TTENSR", 7, 10, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "T23:"},
    { "TTENSR", 8, 24, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "L TENSOR"}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "", 0, 6, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "L11:"},
    { "LTENSR", 3, 10, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "L22:"},
    { "LTENSR", 4, 24, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "", 0, 6, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "L33:"},
    { "LTENSR", 5, 10, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "L12:"},
    { "LTENSR", 6, 24, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "", 0, 6, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "L13:"},
    { "LTENSR", 7, 10, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "L23:"},
    { "LTENSR", 8, 24, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S TENSOR"}}
  },
  { 3, 6, PDB_REMARK_GENERAL, 0, {
    { "", 0, 6, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S11:"},
    { "STENSR", 3, 10, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S12:"},
    { "STENSR", 4, 24, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 34, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S13:"},
    { "STENSR", 5, 38, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 6, PDB_REMARK_GENERAL, 0, {
    { "", 0, 6, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S21:"},
    { "STENSR", 6, 10, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S22:"},
    { "STENSR", 7, 24, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 34, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S23:"},
    { "STENSR", 8, 38, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 6, PDB_REMARK_GENERAL, 0, {
    { "", 0, 6, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S31:"},
    { "STENSR", 9, 10, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S32:"},
    { "STENSR",10, 24, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 34, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S33:"},
    { "STENSR",11, 38, 2, 9, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  }
};

REMARKS Remark_3_PHENIX_NCS[Num_Remark_3_PHENIX_NCS] = {
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NCS DETAILS"}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NUMBER OF NCS GROUPS :"},
    { "NCSTLS", 2, 26, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  }
};

REMARKS Remark_3_PHENIX_NCS_GROUP =
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NCS GROUP :"},
    { "NCSGRO", 5, 15, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  };

REMARKS Remark_3_PHENIX_NCS_GROUP_OP[Num_Remark_3_PHENIX_NCS_GROUP_OP] = {
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NCS OPERATOR :"},
    { "PHENCS", 3, 19, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "REFERENCE SELECTION:"},
    { "PHENCS", 4, 26, 3, 34, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SELECTION          :"},
    { "PHENCS", 5, 26, 3, 34, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "ATOM PAIRS NUMBER  :"},
    { "PHENCS", 6, 26, 3, 34, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "RMSD               :"},
    { "PHENCS", 7, 26, 3, 34, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  }
};

GROUP_REMARKS Remark_Phenix_Write = 
{
  13,
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0 },
  { "", "", "", "", "", "", "", "", "", "", "TLSGRO", "", "" },
  { &Remark_3_REFMAC_TARGET, Remark_3_Phenix1, Remark_3_Phenix2, Remark_3_Phenix_BIN,
    Remark_3_PHENIX_BULK, Remark_3_PHENIX_ESU, Remark_3_BValue, Remark_3_PHENIX_TWINNING,
    Remark_3_PHENIX_RMSD, Remark_3_REFMAC_TLS, Remark_3_PHENIX_TLS_GROUP,
    Remark_3_PHENIX_NCS, &Remark_3_Other
  },
  { 1, Num_Remark_3_Phenix1, Num_Remark_3_Phenix2, Num_Remark_3_Phenix_BIN,
    Num_Remark_3_PHENIX_BULK, Num_Remark_3_PHENIX_ESU, Num_Remark_3_BValue,
    Num_Remark_3_PHENIX_TWINNING, Num_Remark_3_PHENIX_RMSD,
    Num_Remark_3_REFMAC_TLS, Num_Remark_3_PHENIX_TLS_GROUP,
    Num_Remark_3_PHENIX_NCS, 1
  }
};

GROUP_REMARKS Remark_Phenix_Read = 
{
  15,
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 1, 0 },
  { "", "", "", "", "", "", "", "", "", "", "TLSGRO", "", "", "", "" },
  { &Remark_3_REFMAC_TARGET, Remark_3_Phenix1, Remark_3_Phenix2, Remark_3_Phenix_BIN,
    Remark_3_PHENIX_BULK, Remark_3_PHENIX_ESU, Remark_3_BValue_PHENIX, Remark_3_PHENIX_TWINNING,
    Remark_3_PHENIX_RMSD, Remark_3_REFMAC_TLS, Remark_3_PHENIX_TLS_GROUP,
    Remark_3_PHENIX_NCS, &Remark_3_PHENIX_NCS_GROUP, Remark_3_PHENIX_NCS_GROUP_OP,
    &Remark_3_Other
  },
  { 1, Num_Remark_3_Phenix1, Num_Remark_3_Phenix2, Num_Remark_3_Phenix_BIN,
    Num_Remark_3_PHENIX_BULK, Num_Remark_3_PHENIX_ESU, Num_Remark_3_BValue_PHENIX,
    Num_Remark_3_PHENIX_TWINNING, Num_Remark_3_PHENIX_RMSD,
    Num_Remark_3_REFMAC_TLS, Num_Remark_3_PHENIX_TLS_GROUP,
    Num_Remark_3_PHENIX_NCS, 1, Num_Remark_3_PHENIX_NCS_GROUP_OP, 1
  }
};

const char *_method_remark[NUM_METHOD] = {
       "X-RAY DATA.",
       "NEUTRON DATA."
};

const char *_method_name[NUM_METHOD] = {
       "X-RAY DIFFRACTION",
       "NEUTRON DIFFRACTION"
};
