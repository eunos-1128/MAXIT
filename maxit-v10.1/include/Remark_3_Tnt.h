/*
FILE:     Remark_3_Tnt.h
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
REMARKS Remark_3_General_TNT[Num_Remark_3_General_TNT] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"USING DATA ABOVE SIGMA CUTOFF."}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"CROSS-VALIDATION METHOD          :"},
    { "FREERF",  7,   38, 3, 23, 0,  PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FREE R VALUE TEST SET SELECTION  :"},
    { "FREERF",  3,   38, 3, 23, 0,  PDB_REMARK_LEFT_JUSTIFIED,
	""}}
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
    { "RFACTR",  11,  38, 2, 6, 3,  PDB_REMARK_LEFT_JUSTIFIED,
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
  },
};

REMARKS Remark_3_General_Buster_TNT[Num_Remark_3_General_Buster_TNT] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FIT TO DATA USED IN REFINEMENT."}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"CROSS-VALIDATION METHOD           :"},
    { "FREERF",  7,   39, 3, 23, 0,  PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FREE R VALUE TEST SET SELECTION   :"},
    { "FREERF",  3,   39, 3, 23, 0,  PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_R_VALUE, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"R VALUE     (WORKING + TEST SET)  :"},
    { "RFACTR",  8,   39, 2, 6, 3,  PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_R_VALUE, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"R VALUE            (WORKING SET)  :"},
    { "RFACTR",  11,  39, 2, 6, 3,  PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_R_VALUE, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FREE R VALUE                      :"},
    { "FREERF",  2,   39, 2, 5, 3,  PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FREE R VALUE TEST SET SIZE   (%)  :"},
    { "FREERF",  4,   39, 2, 6, 3,  PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FREE R VALUE TEST SET COUNT       :"},
    { "FREERF",  5,   39, 1, 5, 0,  PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "ESTIMATED ERROR OF FREE R VALUE   :"},
    { "FREERF",  6,   39, 2, 6, 3,  PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  }
};

REMARKS Remark_3_TNT_Fit_Agreement[Num_Remark_3_Fit_Agreement] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "USING ALL DATA, NO SIGMA CUTOFF."}}
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
    { "RNOCUT",    5, 47, 2, 6, 2, PDB_REMARK_LEFT_JUSTIFIED,
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

REMARKS Remark_3_Tnt_Wilson = 
{   3,  2, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"WILSON B VALUE (FROM FCALC, A**2) :"},
    { "BVALUE", 3, 38, 2, 7, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
};

REMARKS Remark_3_Tnt_RMSD[Num_Remark_3_Tnt_RMSD] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"RMS DEVIATIONS FROM IDEAL VALUES.    RMS    WEIGHT  COUNT"}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BOND LENGTHS                 (A) :"},
    { "RMSTN1", 2, 38, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN1", 3, 46, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN1", 4, 54, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BOND ANGLES            (DEGREES) :"},
    { "RMSTN1", 5, 38, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN1", 6, 46, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN1", 7, 54, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"TORSION ANGLES         (DEGREES) :"},
    { "RMSTN1", 8, 38, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN1", 9, 46, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN1", 10, 54, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"PSEUDOROTATION ANGLES  (DEGREES) :"},
    { "RMSTN1", 11, 38, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN1", 12, 46, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN1", 13, 54, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"TRIGONAL CARBON PLANES       (A) :"},
    { "RMSTN2", 2, 38, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN2", 3, 46, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN2", 4, 54, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"GENERAL PLANES               (A) :"},
    { "RMSTN2", 5, 38, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN2", 6, 46, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN2", 7, 54, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"ISOTROPIC THERMAL FACTORS (A**2) :"},
    { "RMSTN2", 8, 38, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN2", 9, 46, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN2", 10, 54, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"NON-BONDED CONTACTS          (A) :"},
    { "RMSTN2", 11, 38, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN2", 12, 46, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN2", 13, 54, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_Tnt_Chiral = 
{   3,  2, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"INCORRECT CHIRAL-CENTERS (COUNT) :"},
    { "RMSVOL", 4, 37, 1, 5, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
};

REMARKS Remark_3_Tnt_Restraint[Num_Remark_3_Tnt_Restraint] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"RESTRAINT LIBRARIES."}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"STEREOCHEMISTRY :"},
    { "XFILE3", 2, 21, 3, 39, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"ISOTROPIC THERMAL FACTOR RESTRAINTS :"},
    { "RMSGEN", 7, 41, 3, 19, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_Buster_TNT_Fit_BIN[Num_Remark_3_Buster_TNT_Fit_BIN] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "FIT IN THE HIGHEST RESOLUTION BIN."}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "TOTAL NUMBER OF BINS USED               :" },
    { "SHELL",   12, 45, 1, 5, 0, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "BIN RESOLUTION RANGE HIGH   (ANGSTROMS) :" },
    { "SHELL",    2, 45, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "BIN RESOLUTION RANGE LOW    (ANGSTROMS) :" },
    { "SHELL",    3, 45, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "BIN COMPLETENESS (WORKING+TEST)     (%) :"},
    { "SHELL",    4, 45, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "REFLECTIONS IN BIN (WORKING + TEST SET) :"},
    { "SHELL",   13, 45, 1, 10, 0, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3,  2, PDB_REMARK_R_VALUE, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "BIN R VALUE        (WORKING + TEST SET) :"},
    { "SHELL",   14, 45, 2, 6, 4, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "REFLECTIONS IN BIN        (WORKING SET) :"},
    { "SHELL",    5, 45, 1, 10, 0, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3,  2, PDB_REMARK_R_VALUE, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "BIN R VALUE               (WORKING SET) :"},
    { "SHELL",    6, 45, 2, 6, 4, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3,  2, PDB_REMARK_R_VALUE, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "BIN FREE R VALUE                        :"},
    { "SHELL",   7, 45, 2, 6, 4, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "BIN FREE R VALUE TEST SET SIZE      (%) :"},
    { "SHELL",    8, 45, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "BIN FREE R VALUE TEST SET COUNT         :"},
    { "SHELL",    9, 45, 1, 5, 0, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "ESTIMATED ERROR OF BIN FREE R VALUE     :"},
    { "SHELL",    10, 45, 2, 5, 3, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  }
};

REMARKS Remark_3_Buster_TNT_Fit_BIN_X[Num_Remark_3_Buster_TNT_Fit_BIN_X] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "FIT IN THE HIGHEST RESOLUTION BIN."}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "TOTAL NUMBER OF BINS USED           :" },
    { "SHELL",   12, 41, 1, 5, 0, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "TOTAL NUMBER OF BINS USED               :" },
    { "SHELL",   12, 45, 1, 5, 0, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "BIN RESOLUTION RANGE HIGH   (ANGSTROMS) :" },
    { "SHELL",    2, 45, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "BIN RESOLUTION RANGE HIGH       (A) :" },
    { "SHELL",    2, 41, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "BIN RESOLUTION RANGE LOW    (ANGSTROMS) :" },
    { "SHELL",    3, 45, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "BIN RESOLUTION RANGE LOW        (A) :" },
    { "SHELL",    3, 41, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "BIN COMPLETENESS (WORKING+TEST)     (%) :"},
    { "SHELL",    4, 45, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "REFLECTIONS IN BIN (WORKING + TEST SET) :"},
    { "SHELL",   13, 45, 1, 10, 0, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "REFLECTIONS IN BIN   (WORKING+TEST) :"},
    { "SHELL",   13, 41, 1, 10, 0, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3,  2, PDB_REMARK_R_VALUE, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "BIN R VALUE        (WORKING + TEST SET) :"},
    { "SHELL",   14, 45, 2, 6, 4, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3,  2, PDB_REMARK_R_VALUE, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "BIN R VALUE          (WORKING+TEST) :" },
    { "SHELL",   14, 41, 2, 6, 4, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "REFLECTIONS IN BIN        (WORKING SET) :"},
    { "SHELL",    5, 45, 1, 10, 0, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3,  2, PDB_REMARK_R_VALUE, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "BIN R VALUE               (WORKING SET) :"},
    { "SHELL",    6, 45, 2, 6, 4, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3,  2, PDB_REMARK_R_VALUE, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "BIN FREE R VALUE                        :"},
    { "SHELL",   7, 45, 2, 6, 4, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "BIN FREE R VALUE TEST SET SIZE      (%) :"},
    { "SHELL",    8, 45, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "BIN FREE R VALUE TEST SET COUNT         :"},
    { "SHELL",    9, 45, 1, 5, 0, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
       "ESTIMATED ERROR OF BIN FREE R VALUE     :"},
    { "SHELL",    10, 45, 2, 5, 3, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  }
};

REMARKS Remark_3_Buster_ESD[Num_Remark_3_Buster_ESD] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"ESTIMATED COORDINATE ERROR."}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "ESD FROM LUZZATI PLOT                    (A) :" },
    { "REFERR", 2, 50, 2, 5, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "DPI (BLOW EQ-10) BASED ON R VALUE        (A) :" },
    { "REFDPI", 4, 50, 2, 5, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DPI (BLOW EQ-9) BASED ON FREE R VALUE    (A) :" },
    { "REFDPI", 5, 50, 2, 5, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DPI (CRUICKSHANK) BASED ON R VALUE       (A) :" },
    { "REFDPI", 2, 50, 2, 5, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DPI (CRUICKSHANK) BASED ON FREE R VALUE  (A) :" },
    { "REFDPI", 3, 50, 2, 5, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_Buster_Reference[Num_Remark_3_Buster_Reference] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"REFERENCES: BLOW, D. (2002) ACTA CRYST D58, 792-797"}}
  },
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"            CRUICKSHANK, D.W.J. (1999) ACTA CRYST D55, 583-601"}}
  }
};

REMARKS Remark_3_Buster_Tnt_RMSD[Num_Remark_3_Buster_Tnt_RMSD] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NUMBER OF GEOMETRIC FUNCTION TERMS DEFINED : 15" }}
  },
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TERM                          COUNT    WEIGHT   FUNCTION." }}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BOND LENGTHS              :" },
    { "RMSTN1", 4, 32, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 39, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN1", 3, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN1", 14, 50, 3, 16, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BOND ANGLES               :" },
    { "RMSTN1", 7, 32, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 39, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN1", 6, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN1", 15, 50, 3, 16, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TORSION ANGLES            :" },
    { "RMSTN1", 10, 32, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 39, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN1", 9, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN1", 16, 50, 3, 16, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TRIGONAL CARBON PLANES    :" },
    { "RMSTN2", 4, 32, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 39, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN2", 3, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN2", 14, 50, 3, 16, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "GENERAL PLANES            :" },
    { "RMSTN2", 7, 32, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 39, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN2", 6, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN2", 15, 50, 3, 16, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "ISOTROPIC THERMAL FACTORS :" },
    { "RMSTN2", 10, 32, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 39, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN2", 9, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN2", 16, 50, 3, 16, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BAD NON-BONDED CONTACTS   :" },
    { "RMSTN2", 13, 32, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 39, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN2", 12, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN2", 17, 50, 3, 16, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "IMPROPER TORSIONS         :" },
    { "RMSTN3",  4, 32, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 39, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN3",  3, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN3",  5, 50, 3, 16, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "PSEUDOROTATION ANGLES     :" },
    { "RMSTN1", 13, 32, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 39, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN1", 12, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN1", 17, 50, 3, 16, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "CHIRAL IMPROPER TORSION   :" },
    { "RMSTN3",  8, 32, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 39, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN3",  7, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN3",  9, 50, 3, 16, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SUM OF OCCUPANCIES        :" },
    { "RMSTN3", 12, 32, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 39, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN3", 11, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN3", 13, 50, 3, 16, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "UTILITY DISTANCES         :" },
    { "RMSTN4", 4, 32, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 39, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN4", 3, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN4", 14, 50, 3, 16, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "UTILITY ANGLES            :" },
    { "RMSTN4", 7, 32, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 39, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN4", 6, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN4", 15, 50, 3, 16, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "UTILITY TORSION           :" },
    { "RMSTN4", 10, 32, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 39, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN4", 9, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN4", 16, 50, 3, 16, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "IDEAL-DIST CONTACT TERM   :" },
    { "RMSTN4", 13, 32, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 39, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN4", 12, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSTN4", 17, 50, 3, 16, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_Buster_Tnt_RMSD1[Num_Remark_3_Buster_Tnt_RMSD1] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "RMS DEVIATIONS FROM IDEAL VALUES."}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BOND LENGTHS                       (A) :" },
    { "RMSTN1", 2, 45, 2, 5, 3, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BOND ANGLES                  (DEGREES) :" },
    { "RMSTN1", 5, 45, 2, 6, 2, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "PEPTIDE OMEGA TORSION ANGLES (DEGREES) :" },
    { "RMSTN3", 14, 45, 2, 6, 2, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "OTHER TORSION ANGLES         (DEGREES) :" },
    { "RMSTN3", 15, 45, 2, 6, 2, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  }
};

REMARKS Remark_3_BUSTER_TLS_GROUP[Num_Remark_3_BUSTER_TLS_GROUP] = {
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
    { "TLSGRO", 4, 29, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "TLSGRO", 5, 39, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "TLSGRO", 6, 49, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "T TENSOR"}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "T11:"},
    { "TTENSR", 3,  9, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "T22:"},
    { "TTENSR", 4, 24, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "T33:"},
    { "TTENSR", 5,  9, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "T12:"},
    { "TTENSR", 6, 24, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "T13:"},
    { "TTENSR", 7,  9, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "T23:"},
    { "TTENSR", 8, 24, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "L TENSOR"}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "L11:"},
    { "LTENSR", 3,  9, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "L22:"},
    { "LTENSR", 4, 24, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "L33:"},
    { "LTENSR", 5,  9, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "L12:"},
    { "LTENSR", 6, 24, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "L13:"},
    { "LTENSR", 7,  9, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "L23:"},
    { "LTENSR", 8, 24, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S TENSOR"}}
  },
  { 3, 6, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S11:"},
    { "STENSR", 3,  9, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S12:"},
    { "STENSR", 4, 24, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 35, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S13:"},
    { "STENSR", 5, 39, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 6, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S21:"},
    { "STENSR", 6,  9, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S22:"},
    { "STENSR", 7, 24, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 35, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S23:"},
    { "STENSR", 8, 39, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 6, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S31:"},
    { "STENSR", 9,  9, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S32:"},
    { "STENSR",10, 24, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 35, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S33:"},
    { "STENSR",11, 39, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  }
};

REMARKS Remark_3_BUSTER_TLS_GROUP_Read[Num_Remark_3_BUSTER_TLS_GROUP_Read] = {
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
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SELECTION :"},
    { "TLSGRO", 8, 16, 3, 44, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SET :"},
    { "TLSGRO", 8, 10, 3, 55, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "ORIGIN FOR THE GROUP (A):"},
    { "TLSGRO", 4, 29, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "TLSGRO", 5, 39, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "TLSGRO", 6, 49, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "T TENSOR"}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "T11:"},
    { "TTENSR", 3,  9, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "T22:"},
    { "TTENSR", 4, 24, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "T33:"},
    { "TTENSR", 5,  9, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "T12:"},
    { "TTENSR", 6, 24, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "T13:"},
    { "TTENSR", 7,  9, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "T23:"},
    { "TTENSR", 8, 24, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "L TENSOR"}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "L11:"},
    { "LTENSR", 3,  9, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "L22:"},
    { "LTENSR", 4, 24, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "L33:"},
    { "LTENSR", 5,  9, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "L12:"},
    { "LTENSR", 6, 24, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "L13:"},
    { "LTENSR", 7,  9, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "L23:"},
    { "LTENSR", 8, 24, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S TENSOR"}}
  },
  { 3, 6, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S11:"},
    { "STENSR", 3,  9, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S12:"},
    { "STENSR", 4, 24, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 35, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S13:"},
    { "STENSR", 5, 39, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 6, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S21:"},
    { "STENSR", 6,  9, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S22:"},
    { "STENSR", 7, 24, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 35, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S23:"},
    { "STENSR", 8, 39, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 6, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S31:"},
    { "STENSR", 9,  9, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 20, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S32:"},
    { "STENSR",10, 24, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 35, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "S33:"},
    { "STENSR",11, 39, 2, 10, 4, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  }
};

GROUP_REMARKS Remark_Tnt  =
{
  10,
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { "", "", "", "", "", "", "", "", "", "" },
  { Remark_3_General,  Remark_3_General_TNT,
    Remark_3_TNT_Fit_Agreement, Remark_3_AtomCount,
    &Remark_3_Tnt_Wilson, Remark_3_Tnt_RMSD,
    &Remark_3_Tnt_Chiral, Remark_3_Solvent,
    Remark_3_Tnt_Restraint, &Remark_3_Other
  },
  { Num_Remark_3_General, Num_Remark_3_General_TNT,
    Num_Remark_3_Fit_Agreement, Num_Remark_3_AtomCount,
    1, Num_Remark_3_Tnt_RMSD, 1, Num_Remark_3_Solvent,
    Num_Remark_3_Tnt_Restraint, 1
  }
};

GROUP_REMARKS Remark_Buster_Tnt = 
{
  13,
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0 },
  { "", "", "", "", "", "", "", "", "", "", "", "TLSGRO", "" },
  { Remark_3_General,  Remark_3_General_Buster_TNT,
    Remark_3_Buster_TNT_Fit_BIN, Remark_3_AtomCount,
    Remark_3_BValue, Remark_3_Buster_ESD, Remark_3_Buster_Reference,
    Remark_3_REFMAC_CORR, Remark_3_Buster_Tnt_RMSD, Remark_3_Buster_Tnt_RMSD1,
    Remark_3_REFMAC_TLS, Remark_3_BUSTER_TLS_GROUP, &Remark_3_Other
  },
  { Num_Remark_3_General, Num_Remark_3_General_Buster_TNT,
    Num_Remark_3_Buster_TNT_Fit_BIN, Num_Remark_3_AtomCount,
    Num_Remark_3_BValue, Num_Remark_3_Buster_ESD, Num_Remark_3_Buster_Reference,
    Num_Remark_3_REFMAC_CORR, Num_Remark_3_Buster_Tnt_RMSD,
    Num_Remark_3_Buster_Tnt_RMSD1, Num_Remark_3_REFMAC_TLS,
    Num_Remark_3_BUSTER_TLS_GROUP, 1
  }
};

GROUP_REMARKS Remark_Buster_Tnt_Read = 
{
  13,
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0 },
  { "", "", "", "", "", "", "", "", "", "", "", "TLSGRO", "" },
  { Remark_3_General,  Remark_3_General_Buster_TNT,
    Remark_3_Buster_TNT_Fit_BIN_X, Remark_3_AtomCount,
    Remark_3_BValue, Remark_3_Buster_ESD, Remark_3_Buster_Reference,
    Remark_3_REFMAC_CORR, Remark_3_Buster_Tnt_RMSD, Remark_3_Buster_Tnt_RMSD1,
    Remark_3_REFMAC_TLS, Remark_3_BUSTER_TLS_GROUP_Read, &Remark_3_Other
  },
  { Num_Remark_3_General, Num_Remark_3_General_Buster_TNT,
    Num_Remark_3_Buster_TNT_Fit_BIN_X, Num_Remark_3_AtomCount,
    Num_Remark_3_BValue, Num_Remark_3_Buster_ESD, Num_Remark_3_Buster_Reference,
    Num_Remark_3_REFMAC_CORR, Num_Remark_3_Buster_Tnt_RMSD,
    Num_Remark_3_Buster_Tnt_RMSD1, Num_Remark_3_REFMAC_TLS,
    Num_Remark_3_BUSTER_TLS_GROUP_Read, 1
  }
};
