/*
FILE:     Remark_General.h
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
REMARKS Remark_2 = {
     2,  3, PDB_REMARK_GENERAL, 0, {
      { "",                0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED, 
	  "RESOLUTION." }, 
      { "RFACTR",     7, 13, 2, 7, 2, PDB_REMARK_RIGHT_JUSTIFIED,
	  "Up_Lim_Resol_Ref"},
      { "",                 0, 21, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	  "ANGSTROMS."}}
};

REMARKS Remark_2_NMR = {
     2, 1, PDB_REMARK_GENERAL, 0, {
      {"",                0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
         "RESOLUTION. NOT APPLICABLE." }}
};

REMARKS Remark_3_1[Num_Remark_3_1] = {
  { 3,  1,PDB_REMARK_GENERAL, 0,  {
    { "",             0,    1, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"REFINEMENT."}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"PROGRAM     :" },
    { "REFMET",    4, 17, 3, 33, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"AUTHORS     :" },
    { "",    0, 17, 3, 33, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
};

REMARKS Remark_3_2[Num_Remark_3_2] = {
  { 3,  1,PDB_REMARK_GENERAL, 0,  {
    { "",             0,    1, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"REFINEMENT."}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"PROGRAM 1   :" },
    { "REFMET",    4, 17, 3, 33, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"AUTHORS 1   :" },
    { "",    0, 17, 3, 33, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"PROGRAM 2   :" },
    { "REFMET",    4, 17, 3, 33, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"AUTHORS 2   :" },
    { "",    0, 17, 3, 33, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_3[Num_Remark_3_3] = {
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "PROGRAM" },
    { "REFMET",    4, 30, 3, 30, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "AUTHORS" },
    { "",    0, 30, 3, 30, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  }
};

REMARKS Remark_3_NMR[Num_Remark_3_NMR] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,  {
    { "",             0,    1, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "REFINEMENT."}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "PROGRAM     :" },
    { "REFMET",    4, 17, 3, 43, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "AUTHORS     :" },
    { "REFMET",    5, 17, 3, 43, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "OTHER REFINEMENT REMARKS:" },
    { "REFREM",    3, 28, 3, 32, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  }
};

REMARKS Remark_3_General[Num_Remark_3_General] = {
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
        "DATA CUTOFF            (SIGMA(F)) :"},
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

REMARKS Remark_3_Fit_Agreement[Num_Remark_3_Fit_Agreement] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FIT/AGREEMENT OF MODEL WITH ALL DATA."}}
  },
  { 3, 2, PDB_REMARK_R_VALUE, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"R VALUE   (WORKING + TEST SET, NO CUTOFF) :"},
    { "RNOCUT",    2, 47, 2, 6, 4, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_R_VALUE, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"R VALUE          (WORKING SET, NO CUTOFF) :"},
    { "RNOCUT",    3, 47, 2, 6, 4, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_R_VALUE, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FREE R VALUE                  (NO CUTOFF) :"},
    { "RNOCUT",    4, 47, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FREE R VALUE TEST SET SIZE (%, NO CUTOFF) :"},
    { "RNOCUT",    5, 47, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FREE R VALUE TEST SET COUNT   (NO CUTOFF) :"},
    { "RNOCUT",    6, 47, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"TOTAL NUMBER OF REFLECTIONS   (NO CUTOFF) :"},
    { "RNOCUT",    7, 47, 2, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
};

REMARKS Remark_3_BValue[Num_Remark_3_BValue] = {
  { 3,  1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"B VALUES."}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "B VALUE TYPE :" },
    { "BVALUE", 10, 18, 3, 30, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
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
	"B22 (A**2) :" },
    { "BVALUE", 5, 17, 2, 10, 5, PDB_REMARK_LEFT_JUSTIFIED,
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
	"B12 (A**2) :" },
    { "BVALUE", 7, 17, 2, 10, 5, PDB_REMARK_LEFT_JUSTIFIED,
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
	"B23 (A**2) :" },
    { "BVALUE", 9, 17, 2, 10, 5, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_ESD[Num_Remark_3_ESD] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"ESTIMATED COORDINATE ERROR."}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"ESD FROM LUZZATI PLOT        (A) :" },
    { "REFERR", 2, 38, 2, 4, 2, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,

	"ESD FROM SIGMAA              (A) :" },
    { "REFERR", 3, 38, 2, 4, 2, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"LOW RESOLUTION CUTOFF        (A) :" },
    { "REFERR", 4, 38, 2, 4, 2, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_AtomCount[Num_Remark_3_AtomCount] = {
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT."}
  }
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",              0,  3, 3, 0, 0,PDB_REMARK_LEFT_JUSTIFIED,
        "PROTEIN ATOMS            :"},
    {"ATNUMS",    3,  30, 1, 6, 0,PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",              0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NUCLEIC ACID ATOMS       :"},
    {"ATNUMS",    4,  30, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",              0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "HETEROGEN ATOMS          :"},
    {"ATNUMS",    5,  30, 1, 5, 0, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",              0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SOLVENT ATOMS            :"},
    {"ATNUMS",    6,  30, 1, 5, 0, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  }
};

REMARKS Remark_3_Isotropic[Num_Remark_3_Isotropic] = {
  { 3, 1,PDB_REMARK_GENERAL, 0, {
    { "",               0,  2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"ISOTROPIC THERMAL FACTOR RESTRAINTS.    RMS    SIGMA"}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"MAIN-CHAIN BOND              (A**2) :"},
    { "RMSISO", 2, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 47, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSISO", 3, 49, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"MAIN-CHAIN ANGLE             (A**2) :"},
    { "RMSISO", 4, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 47, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSISO", 5, 49, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SIDE-CHAIN BOND              (A**2) :"},
    { "RMSISO", 6, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 47, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSISO", 7, 49, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SIDE-CHAIN ANGLE             (A**2) :"},
    { "RMSISO", 8, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 47, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSISO", 9, 49, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_Solvent[Num_Remark_3_Solvent] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BULK SOLVENT MODELING."}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "METHOD USED :"},
    { "SOLMOD", 2, 17, 3, 43, 3, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "KSOL        :"},
    { "SOLMOD", 3, 17, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BSOL        :"},
    { "SOLMOD", 4, 17, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  }
};

REMARKS Remark_3_Other = {
    3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"OTHER REFINEMENT REMARKS:"},
    { "REFREM",               3,  28, 3, 32, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"OTHER REFINEMENT REMARKS:"} }
};

REMARKS Remark_4 = {
    4,  2,  PDB_REMARK_GENERAL, 0, {
    { "PDBFIL",  2,  1,  3, 4, 0, PDB_REMARK_LEFT_JUSTIFIED, 
	  ""},
    { "",  0, 6,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED, 
	  "COMPLIES WITH FORMAT V. 3.30, 13-JUL-11"},
    }
};

REMARKS Remark_4_V315 = {
    4,  2,  PDB_REMARK_GENERAL, 0, {
      { "PDBFIL",  2,  1,  3, 4, 0, PDB_REMARK_LEFT_JUSTIFIED, 
	  ""},
      { "",  0, 6,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED, 
          "COMPLIES WITH FORMAT V. 3.15, 01-DEC-08"},
    }
};
