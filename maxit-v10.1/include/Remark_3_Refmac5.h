/*
FILE:     Remark_3_Refmac5.h
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
REMARKS Remark_3_REFMAC_TARGET = {
    3,  2, PDB_REMARK_GENERAL, 0,{
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "REFINEMENT TARGET :"},
    { "XFILE3",    2, 24, 3, 36, 0, PDB_REMARK_LEFT_JUSTIFIED,
       ""} }
};

REMARKS Remark_3_REFMAC_BIN[Num_Remark_3_REFMAC_BIN] = {
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
	"BIN RESOLUTION RANGE HIGH       (A) :" },
    { "SHELL",    2, 41, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
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
	"REFLECTION IN BIN     (WORKING SET) :"},
    { "SHELL",    5, 41, 1, 10, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BIN COMPLETENESS (WORKING+TEST) (%) :"},
    { "SHELL",    4, 41, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3,  2, PDB_REMARK_R_VALUE, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BIN R VALUE           (WORKING SET) :"},
    { "SHELL",    6, 41, 2, 6, 4, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BIN FREE R VALUE SET COUNT          :"},
    { "SHELL",    9, 41, 1, 5, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_R_VALUE, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BIN FREE R VALUE                    :"},
    { "SHELL",   7, 41, 2, 6, 4, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_REFMAC_BIN_Read[Num_Remark_3_REFMAC_BIN_Read] = {
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
	"BIN RESOLUTION RANGE HIGH           :" },
    { "SHELL",    2, 41, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BIN RESOLUTION RANGE LOW            :" },
    { "SHELL",    3, 41, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
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
	"BIN RESOLUTION RANGE LOW        (A) :" },
    { "SHELL",    3, 41, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"REFLECTION IN BIN     (WORKING SET) :"},
    { "SHELL",    5, 41, 1, 10, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BIN COMPLETENESS (WORKING+TEST) (%) :"},
    { "SHELL",    4, 41, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3,  2, PDB_REMARK_R_VALUE, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BIN R VALUE           (WORKING SET) :"},
    { "SHELL",    6, 41, 2, 6, 4, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BIN FREE R VALUE SET COUNT          :"},
    { "SHELL",    9, 41, 1, 5, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_R_VALUE, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BIN FREE R VALUE                    :"},
    { "SHELL",   7, 41, 2, 6, 4, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_REFMAC_AtomCount[Num_Remark_3_REFMAC_AtomCount] = {
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT."}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",              0,  3, 3, 0, 0,PDB_REMARK_LEFT_JUSTIFIED,
        "ALL ATOMS                :"},
    {"ATNUMS",    7,  30, 1, 6, 0,PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  }
};

REMARKS Remark_3_REFMAC_CORR[Num_Remark_3_REFMAC_CORR] = {
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 1, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "CORRELATION COEFFICIENTS."}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",              0,  3, 3, 0, 0,PDB_REMARK_LEFT_JUSTIFIED,
        "CORRELATION COEFFICIENT FO-FC      :"},
    {"COFOFC",    2,  40, 2, 6, 3,PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",              0,  3, 3, 0, 0,PDB_REMARK_LEFT_JUSTIFIED,
        "CORRELATION COEFFICIENT FO-FC FREE :"},
    {"COFOFC",    3,  40, 2, 6, 3,PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
  }
};

REMARKS Remark_3_REFMAC_RMSD[Num_Remark_3_REFMAC_RMSD] = {
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "RMS DEVIATIONS FROM IDEAL VALUES        COUNT    RMS    WEIGHT"}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BOND LENGTHS REFINED ATOMS        (A):"},
    { "RMSRF1", 2, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 3, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 4, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BOND LENGTHS OTHERS               (A):"},
    { "RMSRF1", 5, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 6, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 7, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BOND ANGLES REFINED ATOMS   (DEGREES):"},
    { "RMSRF1", 8, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 9, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 10, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BOND ANGLES OTHERS          (DEGREES):"},
    { "RMSRF1", 11, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 12, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 13, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TORSION ANGLES, PERIOD 1    (DEGREES):"},
    { "RMSRF2", 2, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 3, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 4, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TORSION ANGLES, PERIOD 2    (DEGREES):"},
    { "RMSRF7", 5, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 6, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 7, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TORSION ANGLES, PERIOD 3    (DEGREES):"},
    { "RMSRF2", 5, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 6, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 7, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TORSION ANGLES, PERIOD 4    (DEGREES):"},
    { "RMSRF7", 8, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 9, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 10, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "CHIRAL-CENTER RESTRAINTS       (A**3):"},
    { "RMSRF2", 8, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 9, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 10, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "GENERAL PLANES REFINED ATOMS      (A):"},
    { "RMSRF2", 11, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 12, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 13, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "GENERAL PLANES OTHERS             (A):"},
    { "RMSRF3", 2, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 3, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 4, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NON-BONDED CONTACTS REFINED ATOMS (A):"},
    { "RMSRF3", 5, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 6, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 7, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NON-BONDED CONTACTS OTHERS        (A):"},
    { "RMSRF3", 8, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 9, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 10, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NON-BONDED TORSION REFINED ATOMS  (A):"},
    { "RMSRF8", 2, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 3, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 4, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NON-BONDED TORSION OTHERS         (A):"},
    { "RMSRF3", 11, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 12, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 13, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "H-BOND (X...Y) REFINED ATOMS      (A):"},
    { "RMSRF4", 2, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 3, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 4, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "H-BOND (X...Y) OTHERS             (A):"},
    { "RMSRF6", 11, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF6", 12, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF6", 13, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "POTENTIAL METAL-ION REFINED ATOMS (A):"},
    { "RMSRF8", 5, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 6, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 7, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "POTENTIAL METAL-ION OTHERS        (A):"},
    { "RMSRF8", 8, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 9, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 10, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY VDW REFINED ATOMS        (A):"},
    { "RMSRF4", 5, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 6, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 7, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY VDW OTHERS               (A):"},
    { "RMSRF4", 8, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 9, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 10, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY H-BOND REFINED ATOMS     (A):"},
    { "RMSRF4", 11, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 12, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 13, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY H-BOND OTHERS            (A):"},
    { "RMSRF7", 2, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 3, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 4, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY METAL-ION REFINED ATOMS  (A):"},
    { "RMSRF7", 11, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 12, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 13, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY METAL-ION OTHERS         (A):"},
    { "RMSRF9", 2, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 3, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 4, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  }
};

REMARKS Remark_3_REFMAC_RMSD_X[Num_Remark_3_REFMAC_RMSD_X] = {
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "RMS DEVIATIONS FROM IDEAL VALUES        COUNT    RMS    WEIGHT"}}
  },
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "RMS DEVIATIONS FROM IDEAL VALUES    COUNT    RMS    WEIGHT"}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BOND LENGTHS REFINED ATOMS        (A):"},
    { "RMSRF1", 2, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 3, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 4, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BOND LENGTHS REFINED           (A):"},
    { "RMSRF1", 2, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 3, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 4, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BOND LENGTHS REFINED          (A):"},
    { "RMSRF1", 2, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 3, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 4, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BOND LENGTHS OTHERS               (A):"},
    { "RMSRF1", 5, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 6, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 7, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BOND LENGTHS OTHERS            (A):"},
    { "RMSRF1", 5, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 6, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 7, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BOND LENGTHS OTHERS           (A):"},
    { "RMSRF1", 5, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 6, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 7, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BOND ANGLES REFINED ATOMS   (DEGREES):"},
    { "RMSRF1", 8, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 9, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 10, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BOND ANGLES REFINED      (DEGREES):"},
    { "RMSRF1", 8, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 9, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 10, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BOND ANGLES REFINED     (DEGREES):"},
    { "RMSRF1", 8, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 9, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 10, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BOND ANGLES OTHERS          (DEGREES):"},
    { "RMSRF1", 11, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 12, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 13, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BOND ANGLES OTHERS       (DEGREES):"},
    { "RMSRF1", 11, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 12, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 13, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BOND ANGLES OTHERS      (DEGREES):"},
    { "RMSRF1", 11, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 12, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF1", 13, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TORSION ANGLES, PERIOD 1    (DEGREES):"},
    { "RMSRF2", 2, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 3, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 4, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TORSION ANGLES, PERIOD 1 (DEGREES):"},
    { "RMSRF2", 2, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 3, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 4, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TORSION ANGLES, PERIOD 1(DEGREES):"},
    { "RMSRF2", 2, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 3, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 4, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TORSION ANGLES, PERIOD 2    (DEGREES):"},
    { "RMSRF7", 5, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 6, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 7, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TORSION ANGLES, PERIOD 2 (DEGREES):"},
    { "RMSRF7", 5, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 6, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 7, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TORSION ANGLES, PERIOD 2(DEGREES):"},
    { "RMSRF7", 5, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 6, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 7, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TORSION ANGLES, PERIOD 2   (DEGREES):"},
    { "RMSRF7", 5, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 6, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 7, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TORSION ANGLES, PERIOD 3    (DEGREES):"},
    { "RMSRF2", 5, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 6, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 7, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TORSION ANGLES, PERIOD 3 (DEGREES):"},
    { "RMSRF2", 5, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 6, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 7, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TORSION ANGLES, PERIOD 3(DEGREES):"},
    { "RMSRF2", 5, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 6, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 7, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TORSION ANGLES, PERIOD 4    (DEGREES):"},
    { "RMSRF7", 8, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 9, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 10, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TORSION ANGLES, PERIOD 4 (DEGREES):"},
    { "RMSRF7", 8, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 9, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 10, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TORSION ANGLES, PERIOD 4(DEGREES):"},
    { "RMSRF7", 8, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 9, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 10, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "CHIRAL-CENTER RESTRAINTS       (A**3):"},
    { "RMSRF2", 8, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 9, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 10, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "CHIRAL-CENTER RESTRAINTS    (A**3):"},
    { "RMSRF2", 8, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 9, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 10, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "CHIRAL-CENTER RESTRAINTS   (A**3):"},
    { "RMSRF2", 8, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 9, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 10, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "GENERAL PLANES REFINED ATOMS      (A):"},
    { "RMSRF2", 11, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 12, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 13, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "GENERAL PLANES REFINED         (A):"},
    { "RMSRF2", 11, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 12, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 13, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "GENERAL PLANES REFINED        (A):"},
    { "RMSRF2", 11, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 12, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF2", 13, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "GENERAL PLANES OTHERS             (A):"},
    { "RMSRF3", 2, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 3, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 4, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "GENERAL PLANES OTHERS          (A):"},
    { "RMSRF3", 2, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 3, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 4, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "GENERAL PLANES OTHERS         (A):"},
    { "RMSRF3", 2, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 3, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 4, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NON-BONDED CONTACTS REFINED ATOMS (A):"},
    { "RMSRF3", 5, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 6, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 7, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NON-BONDED CONTACTS REFINED    (A):"},
    { "RMSRF3", 5, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 6, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 7, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NON-BONDED CONTACTS REFINED   (A):"},
    { "RMSRF3", 5, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 6, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 7, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NON-BONDED CONTACTS OTHERS        (A):"},
    { "RMSRF3", 8, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 9, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 10, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NON-BONDED CONTACTS OTHERS     (A):"},
    { "RMSRF3", 8, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 9, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 10, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NON-BONDED CONTACTS OTHERS    (A):"},
    { "RMSRF3", 8, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 9, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 10, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NON-BONDED TORSION REFINED ATOMS  (A):"},
    { "RMSRF8", 2, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 3, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 4, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NON-BONDED TORSION REFINED     (A):"},
    { "RMSRF8", 2, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 3, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 4, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NON-BONDED TORSION REFINED    (A):"},
    { "RMSRF8", 2, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 3, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 4, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NON-BONDED TORSION OTHERS         (A):"},
    { "RMSRF3", 11, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 12, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 13, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NON-BONDED TORSION OTHERS      (A):"},
    { "RMSRF3", 11, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 12, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 13, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NON-BONDED TORSION OTHERS     (A):"},
    { "RMSRF3", 11, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 12, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF3", 13, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "H-BOND (X...Y) REFINED ATOMS      (A):"},
    { "RMSRF4", 2, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 3, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 4, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "H-BOND (X...Y) REFINED         (A):"},
    { "RMSRF4", 2, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 3, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 4, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "H-BOND (X...Y) REFINED        (A):"},
    { "RMSRF4", 2, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 3, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 4, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "H-BOND (X...Y) OTHERS             (A):"},
    { "RMSRF6", 11, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF6", 12, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF6", 13, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "H-BOND (X...Y) OTHERS          (A):"},
    { "RMSRF6", 11, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF6", 12, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF6", 13, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "H-BOND (X...Y) OTHERS         (A):"},
    { "RMSRF6", 11, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF6", 12, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF6", 13, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "POTENTIAL METAL-ION REFINED ATOMS (A):"},
    { "RMSRF8", 5, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 6, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 7, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "POTENTIAL METAL-ION REFINED    (A):"},
    { "RMSRF8", 5, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 6, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 7, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "POTENTIAL METAL-ION REFINED   (A):"},
    { "RMSRF8", 5, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 6, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 7, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "POTENTIAL METAL-ION OTHERS        (A):"},
    { "RMSRF8", 8, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 9, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 10, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "POTENTIAL METAL-ION OTHERS     (A):"},
    { "RMSRF8", 8, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 9, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 10, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "POTENTIAL METAL-ION OTHERS    (A):"},
    { "RMSRF8", 8, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 9, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 10, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY VDW REFINED ATOMS        (A):"},
    { "RMSRF4", 5, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 6, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 7, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY VDW REFINED           (A):"},
    { "RMSRF4", 5, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 6, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 7, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY VDW REFINED          (A):"},
    { "RMSRF4", 5, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 6, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 7, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY VDW OTHERS               (A):"},
    { "RMSRF4", 8, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 9, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 10, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY VDW OTHERS            (A):"},
    { "RMSRF4", 8, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 9, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 10, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY VDW OTHERS           (A):"},
    { "RMSRF4", 8, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 9, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 10, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY H-BOND REFINED ATOMS     (A):"},
    { "RMSRF4", 11, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 12, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 13, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY H-BOND REFINED        (A):"},
    { "RMSRF4", 11, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 12, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 13, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY H-BOND REFINED       (A):"},
    { "RMSRF4", 11, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 12, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF4", 13, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY H-BOND OTHERS            (A):"},
    { "RMSRF7", 2, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 3, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 4, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY H-BOND OTHERS         (A):"},
    { "RMSRF7", 2, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 3, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 4, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY H-BOND OTHERS        (A):"},
    { "RMSRF7", 2, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 3, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 4, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY METAL-ION REFINED ATOMS  (A):"},
    { "RMSRF7", 11, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 12, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 13, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY METAL-ION REFINED     (A):"},
    { "RMSRF7", 11, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 12, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 13, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY METAL-ION REFINED    (A):"},
    { "RMSRF7", 11, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 12, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF7", 13, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },

  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY METAL-ION OTHERS         (A):"},
    { "RMSRF9", 2, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 3, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 4, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY METAL-ION OTHERS      (A):"},
    { "RMSRF9", 2, 38, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 45, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 3, 46, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 53, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 4, 54, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SYMMETRY METAL-ION OTHERS     (A):"},
    { "RMSRF9", 2, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 3, 45, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 4, 53, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  }
};

REMARKS Remark_3_REFMAC_ISOTROPIC[Num_Remark_3_REFMAC_ISOTROPIC] = {
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "ISOTROPIC THERMAL FACTOR RESTRAINTS.     COUNT   RMS    WEIGHT"}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "MAIN-CHAIN BOND REFINED ATOMS  (A**2):"},
    { "RMSRF5", 2, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 3, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 4, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "MAIN-CHAIN BOND OTHER ATOMS    (A**2):"},
    { "RMSRF8", 11, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 12, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 13, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "MAIN-CHAIN ANGLE REFINED ATOMS (A**2):"},
    { "RMSRF5", 5, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 6, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 7, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "MAIN-CHAIN ANGLE OTHER ATOMS   (A**2):"},
    { "RMSRF9", 5, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 6, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 7, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SIDE-CHAIN BOND REFINED ATOMS  (A**2):"},
    { "RMSRF5", 8, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 9, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 10, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SIDE-CHAIN BOND OTHER ATOMS    (A**2):"},
    { "RMSRF9", 8, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 9, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 10, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SIDE-CHAIN ANGLE REFINED ATOMS (A**2):"},
    { "RMSRF5", 11, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 12, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 13, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SIDE-CHAIN ANGLE OTHER ATOMS   (A**2):"},
    { "RMSRF9", 11, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 12, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 13, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "LONG RANGE B REFINED ATOMS     (A**2):"},
    { "RMSRFA", 2, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRFA", 3, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRFA", 4, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "LONG RANGE B OTHER ATOMS       (A**2):"},
    { "RMSRFA", 5, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRFA", 6, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRFA", 7, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  }
};

REMARKS Remark_3_REFMAC_ISOTROPIC_X[Num_Remark_3_REFMAC_ISOTROPIC_X] = {
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "ISOTROPIC THERMAL FACTOR RESTRAINTS.     COUNT   RMS    WEIGHT"}}
  },
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "ISOTROPIC THERMAL FACTOR RESTRAINTS. COUNT   RMS    WEIGHT"}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "MAIN-CHAIN BOND REFINED ATOMS  (A**2):"},
    { "RMSRF5", 2, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 3, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 4, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "MAIN-CHAIN BOND REFINED   (A**2):"},
    { "RMSRF5", 2, 36, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 43, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 3, 44, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 51, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 4, 52, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "MAIN-CHAIN BOND OTHER ATOMS    (A**2):"},
    { "RMSRF8", 11, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 12, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 13, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "MAIN-CHAIN BOND OTHER     (A**2):"},
    { "RMSRF8", 11, 36, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 43, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 12, 44, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 51, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF8", 13, 52, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "MAIN-CHAIN ANGLE REFINED ATOMS (A**2):"},
    { "RMSRF5", 5, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 6, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 7, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "MAIN-CHAIN ANGLE REFINED  (A**2):"},
    { "RMSRF5", 5, 36, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 43, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 6, 44, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 51, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 7, 52, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "MAIN-CHAIN ANGLE OTHER ATOMS   (A**2):"},
    { "RMSRF9", 5, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 6, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 7, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "MAIN-CHAIN ANGLE OTHER ATOMS (A**2)  :"},
    { "RMSRF9", 5, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 6, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 7, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SIDE-CHAIN BOND REFINED ATOMS  (A**2):"},
    { "RMSRF5", 8, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 9, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 10, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SIDE-CHAIN BOND REFINED   (A**2):"},
    { "RMSRF5", 8, 36, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 43, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 9, 44, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 51, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 10, 52, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SIDE-CHAIN BOND OTHER ATOMS    (A**2):"},
    { "RMSRF9", 8, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 9, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 10, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SIDE-CHAIN BOND OTHER ATOMS  (A**2)  :"},
    { "RMSRF9", 8, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 9, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 10, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SIDE-CHAIN ANGLE REFINED ATOMS (A**2):"},
    { "RMSRF5", 11, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 12, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 13, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SIDE-CHAIN ANGLE REFINED  (A**2):"},
    { "RMSRF5", 11, 36, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 43, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 12, 44, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 51, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF5", 13, 52, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SIDE-CHAIN ANGLE OTHER ATOMS   (A**2):"},
    { "RMSRF9", 11, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 12, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 13, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SIDE-CHAIN ANGLE OTHER ATOMS (A**2)  :"},
    { "RMSRF9", 11, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 12, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF9", 13, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "LONG RANGE B REFINED ATOMS     (A**2):"},
    { "RMSRFA", 2, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRFA", 3, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRFA", 4, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "LONG RANGE B REFINED ATOMS (A**2)    :"},
    { "RMSRFA", 2, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRFA", 3, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRFA", 4, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "LONG RANGE B OTHER ATOMS       (A**2):"},
    { "RMSRFA", 5, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRFA", 6, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRFA", 7, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "LONG RANGE B OTHER ATOMS (A**2)      :"},
    { "RMSRFA", 5, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRFA", 6, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRFA", 7, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  }
};

REMARKS Remark_3_REFMAC_ANISOTROPIC[Num_Remark_3_REFMAC_ANISOTROPIC] = {
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 1, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "ANISOTROPIC THERMAL FACTOR RESTRAINTS.    COUNT   RMS   WEIGHT"}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "RIGID-BOND RESTRAINTS          (A**2):"},
    { "RMSRF6", 2, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF6", 3, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF6", 4, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SPHERICITY; FREE ATOMS         (A**2):"},
    { "RMSRF6", 5, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF6", 6, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF6", 7, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SPHERICITY; BONDED ATOMS       (A**2):"},
    { "RMSRF6", 8, 41, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF6", 9, 49, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 56, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RMSRF6", 10, 57, 2, 6, 3, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  }
};

REMARKS Remark_3_REFMAC_TWIN[Num_Remark_3_REFMAC_TWIN] = {
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TWIN DETAILS"}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NUMBER OF TWIN DOMAINS  :"},
    { "NCSTLS", 4, 29, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  }
};

REMARKS Remark_3_REFMAC_TWIN_OP[Num_Remark_3_REFMAC_TWIN_OP] = {
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 6, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TWIN DOMAIN   :"},
    { "TWIN", 5, 22, 3, 37, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 6, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TWIN OPERATOR :"},
    { "TWIN", 2, 22, 3, 37, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 6, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TWIN FRACTION :"},
    { "TWIN", 3, 22, 3, 37, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  }
};

REMARKS Remark_3_REFMAC_NCS[Num_Remark_3_REFMAC_NCS] = {
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NCS RESTRAINTS STATISTICS"}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NUMBER OF DIFFERENT NCS GROUPS :"},
    { "NCSTLS", 2, 36, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  }
};

REMARKS Remark_3_REFMAC_NCS_LOCAL[Num_Remark_3_REFMAC_NCS_LOCAL] = {
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NCS RESTRAINTS STATISTICS"}}
  },
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NCS TYPE: LOCAL"}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NUMBER OF DIFFERENT NCS PAIRS  :"},
    { "NCSTLS", 2, 36, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "GROUP  CHAIN1    RANGE     CHAIN2     RANGE    COUNT RMS  WEIGHT"}}
  }
};

REMARKS Remark_3_REFMAC_NCS_GROUP[Num_Remark_3_REFMAC_NCS_GROUP] = {
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NCS GROUP NUMBER               :"},
    { "NCSGRO", 5, 35, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "CHAIN NAMES                    :"},
    { "NCSGRO", 4, 38, 3, 30, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 5, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NUMBER OF COMPONENTS NCS GROUP :"},
    { "NCSGRO", 3, 38, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 7, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "COMPONENT C  SSSEQI  TO  C   SSSEQI   CODE"}}
  },
  { 3, 12, PDB_REMARK_GENERAL, -3, {
    { "", 0, 1, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "        "},
    { "CCPNCS", 7, 9, 1, 3, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 12, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "    "},
    { "CCPNCS", 3, 16, 3, 2, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 18, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "  "},
    { "CCPNCS", 4, 20, 1, 5, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 26, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "     "},
    { "CCPNCS", 5, 31, 3, 2, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 33, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "   "},
    { "CCPNCS", 6, 36, 1, 5, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 42, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "   "},
    { "CCPNCS", 8, 45, 1, 3, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 19, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "GROUP CHAIN        COUNT   RMS     WEIGHT"}}
  },
  { 3, 10, PDB_REMARK_GENERAL, -2, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TIGHT POSITIONAL"},
    { "RFTPOS", 2, 20, 1, 3, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 23, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "   "},
    { "RFTPOS", 3, 26, 3, 2, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 32, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "(A):"},
    { "RFTPOS", 4, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RFTPOS", 5, 45, 2, 6, 2, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RFTPOS", 6, 53, 2, 6, 2, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 10, PDB_REMARK_GENERAL, -2, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "MEDIUM POSITIONAL"},
    { "RFMPOS", 2, 20, 1, 3, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 23, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "   "},
    { "RFMPOS", 3, 26, 3, 2, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 32, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "(A):"},
    { "RFMPOS", 4, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RFMPOS", 5, 45, 2, 6, 2, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RFMPOS", 6, 53, 2, 6, 2, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 10, PDB_REMARK_GENERAL, -2, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "LOOSE POSITIONAL"},
    { "RFLPOS", 2, 20, 1, 3, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 23, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "   "},
    { "RFLPOS", 3, 26, 3, 2, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 32, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "(A):"},
    { "RFLPOS", 4, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RFLPOS", 5, 45, 2, 6, 2, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RFLPOS", 6, 53, 2, 6, 2, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 10, PDB_REMARK_GENERAL, -2, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TIGHT THERMAL"},
    { "RFTTHR", 2, 20, 1, 3, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 23, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "   "},
    { "RFTTHR", 3, 26, 3, 2, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 29, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "(A**2):"},
    { "RFTTHR", 4, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RFTTHR", 5, 45, 2, 6, 2, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RFTTHR", 6, 53, 2, 6, 2, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 10, PDB_REMARK_GENERAL, -2, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "MEDIUM THERMAL"},
    { "RFMTHR", 2, 20, 1, 3, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 23, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "   "},
    { "RFMTHR", 3, 26, 3, 2, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 29, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "(A**2):"},
    { "RFMTHR", 4, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RFMTHR", 5, 45, 2, 6, 2, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RFMTHR", 6, 53, 2, 6, 2, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  },
  { 3, 10, PDB_REMARK_GENERAL, -2, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "LOOSE THERMAL"},
    { "RFLTHR", 2, 20, 1, 3, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 23, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "   "},
    { "RFLTHR", 3, 26, 3, 2, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 29, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "(A**2):"},
    { "RFLTHR", 4, 37, 1, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 44, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RFLTHR", 5, 45, 2, 6, 2, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "", 0, 52, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ";"},
    { "RFLTHR", 6, 53, 2, 6, 2, PDB_REMARK_RIGHT_JUSTIFIED,
        ""}}
  }
};

REMARKS Remark_3_REFMAC_TLS[Num_Remark_3_REFMAC_TLS] = {
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TLS DETAILS"}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NUMBER OF TLS GROUPS  :"},
    { "NCSTLS", 3, 27, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  } /*,
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "ATOM RECORD CONTAINS RESIDUAL B FACTORS ONLY"}}
  }
*/
};

REMARKS Remark_3_REFMAC_TLS_GROUP[Num_Remark_3_REFMAC_TLS_GROUP] = {
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "TLS GROUP :"},
    { "TLSGRO", 2, 15, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NUMBER OF COMPONENTS GROUP :"},
    { "TLSGRO", 3, 33, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "COMPONENTS        C SSSEQI   TO  C SSSEQI"}}
  },
  { 3, 5, PDB_REMARK_GENERAL, -2, {
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "RESIDUE RANGE :"},
    { "TLSRNG", 3, 21, 3, 2, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "TLSRNG", 4, 24, 1, 5, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "TLSRNG", 5, 36, 3, 2, 0, PDB_REMARK_RIGHT_JUSTIFIED,
        ""},
    { "TLSRNG", 6, 39, 1, 5, 0, PDB_REMARK_RIGHT_JUSTIFIED,
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

REMARKS Remark_3_REFMAC_BULK[Num_Remark_3_REFMAC_BULK] = {
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "BULK SOLVENT MODELLING."}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "METHOD USED :"},
    { "SOLMOD", 2, 17, 3, 43, 3, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "PARAMETERS FOR MASK CALCULATION"}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "VDW PROBE RADIUS   :"},
    { "SOLMOD", 5, 24, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "ION PROBE RADIUS   :"},
    { "SOLMOD", 6, 24, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SHRINKAGE RADIUS   :"},
    { "SOLMOD", 7, 24, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  }
};

GROUP_REMARKS Remark_Refmac_5_0 = 
{
  18,
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, /* 4, */ 0, 2, 0, 0},
  { "", "", "", "", "", "", "", "", "", "", "", "", "NCSGRO", "", /* "", */ "", "TLSGRO", "", "" },
  { &Remark_3_REFMAC_TARGET, Remark_3_General,
    Remark_3_NUCLSQ_PROLSQ, Remark_3_REFMAC_BIN,
    /* Remark_3_REFMAC_AtomCount, */ Remark_3_AtomCount,
    Remark_3_BValue, Remark_3_ESU,
    Remark_3_REFMAC_CORR, Remark_3_REFMAC_RMSD,
    Remark_3_REFMAC_ISOTROPIC, Remark_3_REFMAC_ANISOTROPIC,
    Remark_3_REFMAC_NCS, Remark_3_REFMAC_NCS_GROUP,
    Remark_3_REFMAC_TWIN, // Remark_3_REFMAC_TWIN_OP,
    Remark_3_REFMAC_TLS, Remark_3_REFMAC_TLS_GROUP,
    Remark_3_REFMAC_BULK, &Remark_3_Other
  },
  { 1, Num_Remark_3_General, Num_Remark_3_NUCLSQ_PROLSQ,
    Num_Remark_3_REFMAC_BIN, /* Num_Remark_3_REFMAC_AtomCount, */
    Num_Remark_3_AtomCount,
    Num_Remark_3_BValue, Num_Remark_3_ESU, 
    Num_Remark_3_REFMAC_CORR, Num_Remark_3_REFMAC_RMSD,
    Num_Remark_3_REFMAC_ISOTROPIC, Num_Remark_3_REFMAC_ANISOTROPIC,
    Num_Remark_3_REFMAC_NCS, Num_Remark_3_REFMAC_NCS_GROUP,
    Num_Remark_3_REFMAC_TWIN, // Num_Remark_3_REFMAC_TWIN_OP,
    Num_Remark_3_REFMAC_TLS, Num_Remark_3_REFMAC_TLS_GROUP,
    Num_Remark_3_REFMAC_BULK, 1
  }
};

GROUP_REMARKS Remark_Refmac_5_0_local = 
{
  17,
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0},
  { "", "", "", "", "", "", "", "", "", "", "", "", "", "", "TLSGRO", "", "" },
  { &Remark_3_REFMAC_TARGET, Remark_3_General, Remark_3_NUCLSQ_PROLSQ,
    Remark_3_REFMAC_BIN, Remark_3_AtomCount, Remark_3_BValue, Remark_3_ESU,
    Remark_3_REFMAC_CORR, Remark_3_REFMAC_RMSD, Remark_3_REFMAC_ISOTROPIC,
    Remark_3_REFMAC_ANISOTROPIC, Remark_3_REFMAC_NCS_LOCAL,
    Remark_3_REFMAC_TWIN, Remark_3_REFMAC_TLS, Remark_3_REFMAC_TLS_GROUP,
    Remark_3_REFMAC_BULK, &Remark_3_Other
  },
  { 1, Num_Remark_3_General, Num_Remark_3_NUCLSQ_PROLSQ, Num_Remark_3_REFMAC_BIN,
    Num_Remark_3_AtomCount, Num_Remark_3_BValue, Num_Remark_3_ESU, 
    Num_Remark_3_REFMAC_CORR, Num_Remark_3_REFMAC_RMSD, Num_Remark_3_REFMAC_ISOTROPIC,
    Num_Remark_3_REFMAC_ANISOTROPIC, Num_Remark_3_REFMAC_NCS_LOCAL,
    Num_Remark_3_REFMAC_TWIN, Num_Remark_3_REFMAC_TLS, Num_Remark_3_REFMAC_TLS_GROUP,
    Num_Remark_3_REFMAC_BULK, 1
  }
};

GROUP_REMARKS Remark_Refmac_5_0_Mixed = 
{
  20,
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 4, 0, 2, 0, 0},
  { "", "", "", "", "", "", "", "", "", "", "", "", "", "NCSGRO", "", "", "", "TLSGRO", "", "" },
  { &Remark_3_REFMAC_TARGET, Remark_3_General,
    Remark_3_NUCLSQ_PROLSQ, Remark_3_REFMAC_BIN_Read,
    Remark_3_REFMAC_AtomCount, Remark_3_BValue, Remark_3_ESU,
    Remark_3_REFMAC_CORR, Remark_3_REFMAC_RMSD_X,
    Remark_3_REFMAC_ISOTROPIC_X, Remark_3_REFMAC_ANISOTROPIC,
    Remark_3_REFMAC_NCS_LOCAL, Remark_3_REFMAC_NCS, Remark_3_REFMAC_NCS_GROUP,
    Remark_3_REFMAC_TWIN, Remark_3_REFMAC_TWIN_OP,
    Remark_3_REFMAC_TLS, Remark_3_REFMAC_TLS_GROUP,
    Remark_3_REFMAC_BULK, &Remark_3_Other
  },
  { 1, Num_Remark_3_General, Num_Remark_3_NUCLSQ_PROLSQ,
    Num_Remark_3_REFMAC_BIN_Read, Num_Remark_3_REFMAC_AtomCount,
    Num_Remark_3_BValue, Num_Remark_3_ESU, 
    Num_Remark_3_REFMAC_CORR, Num_Remark_3_REFMAC_RMSD_X,
    Num_Remark_3_REFMAC_ISOTROPIC_X, Num_Remark_3_REFMAC_ANISOTROPIC,
    Num_Remark_3_REFMAC_NCS_LOCAL, Num_Remark_3_REFMAC_NCS, Num_Remark_3_REFMAC_NCS_GROUP,
    Num_Remark_3_REFMAC_TWIN, Num_Remark_3_REFMAC_TWIN_OP,
    Num_Remark_3_REFMAC_TLS, Num_Remark_3_REFMAC_TLS_GROUP,
    Num_Remark_3_REFMAC_BULK, 1
  }
};
