/*
FILE:     Remark_3_Xplors.h
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
REMARKS Remark_3_XPLOR[Num_Remark_3_XPLOR] = {
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
    { "RFACTR",    6, 39, 2, 6, 2, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,{
    { "",              0,  3, 3, 0, 0,PDB_REMARK_LEFT_JUSTIFIED,
        "DATA CUTOFF            (SIGMA(F)) :"},
    {"RFACTR",   5,  39, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
       ""},}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DATA CUTOFF HIGH         (ABS(F)) :" },
    { "RFACTR",   12, 39, 2, 15, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DATA CUTOFF LOW          (ABS(F)) :" },
    { "RFACTR",   10, 39, 2, 7, 4, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"COMPLETENESS (WORKING+TEST)   (%) :"},
    { "PEROBS",  2,  39, 2, 6, 1, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"NUMBER OF REFLECTIONS             :"},
    { "RFACTR",  3,   39, 1, 10, 0,  PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_CNS_Target =
{   3,  2, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "REFINEMENT TARGET :"},
    { "XFILE3",    2, 22, 3, 38, 0, PDB_REMARK_LEFT_JUSTIFIED,
       ""}}
};

REMARKS Remark_3_CNS[Num_Remark_3_CNS] = {
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
    { "RFACTR",    6, 39, 2, 6, 2, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,{
    { "",              0,  3, 3, 0, 0,PDB_REMARK_LEFT_JUSTIFIED,
        "DATA CUTOFF            (SIGMA(F)) :"},
    {"RFACTR",   5,  39, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
       ""},}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"OUTLIER CUTOFF HIGH (RMS(ABS(F))) :" },
    { "RFACTR",   13, 39, 2, 20, 5, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"COMPLETENESS (WORKING+TEST)   (%) :"},
    { "PEROBS",  2,  39, 2, 6, 1, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"NUMBER OF REFLECTIONS             :"},
    { "RFACTR",  3,   39, 1, 10, 0,  PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_CNS_Mixed[Num_Remark_3_CNS_Mixed] = {
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
    { "RFACTR",    6, 39, 2, 6, 2, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,{
    { "",              0,  3, 3, 0, 0,PDB_REMARK_LEFT_JUSTIFIED,
        "DATA CUTOFF            (SIGMA(F)) :"},
    {"RFACTR",   5,  39, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
       ""},}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DATA CUTOFF HIGH         (ABS(F)) :" },
    { "RFACTR",   12, 39, 2, 15, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DATA CUTOFF LOW          (ABS(F)) :" },
    { "RFACTR",   10, 39, 2, 7, 4, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"OUTLIER CUTOFF HIGH (RMS(ABS(F))) :" },
    { "RFACTR",   13, 39, 2, 20, 5, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"COMPLETENESS (WORKING+TEST)   (%) :"},
    { "PEROBS",  2,  39, 2, 6, 1, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
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
  }
};

REMARKS Remark_3_XPLOR_CNS_Fit[Num_Remark_3_XPLOR_CNS_Fit] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FIT TO DATA USED IN REFINEMENT."}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"CROSS-VALIDATION METHOD          :"},
    { "FREERF",  7,   38, 3, 22, 0,  PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FREE R VALUE TEST SET SELECTION  :"},
    { "FREERF",  3,   38, 3, 22, 0,  PDB_REMARK_LEFT_JUSTIFIED,
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
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"ESTIMATED ERROR OF FREE R VALUE  :"},
    { "FREERF",  6,   38, 2, 6, 3,  PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_CNX_Fit[Num_Remark_3_CNX_Fit] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FIT TO DATA USED IN REFINEMENT."}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"CROSS-VALIDATION METHOD          :"},
    { "FREERF",  7,   38, 3, 22, 0,  PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FREE R VALUE TEST SET SELECTION  :"},
    { "FREERF",  3,   38, 3, 22, 0,  PDB_REMARK_LEFT_JUSTIFIED,
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
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"ESTIMATED ERROR OF FREE R VALUE  :"},
    { "FREERF",  6,   38, 2, 6, 3,  PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_CNX_Fit_Agreement[Num_Remark_3_CNX_Fit_Agreement] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FIT/AGREEMENT OF MODEL WITH ALL DATA."}}
  },
  { 3, 2, PDB_REMARK_R_VALUE, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"R VALUE     (WORKING + TEST SET, NO CUTOFF) :"},
    { "RNOCUT",    2, 49, 2, 6, 4, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_R_VALUE, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"R VALUE            (WORKING SET, NO CUTOFF) :"},
    { "RNOCUT",    3, 49, 2, 6, 4, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_R_VALUE, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FREE R VALUE                    (NO CUTOFF) :"},
    { "RNOCUT",    4, 49, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FREE R VALUE TEST SET SIZE   (%, NO CUTOFF) :"},
    { "RNOCUT",    5, 49, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FREE R VALUE TEST SET COUNT     (NO CUTOFF) :"},
    { "RNOCUT",    6, 49, 1, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"ESTIMATED ERROR OF FREE R VALUE (NO CUTOFF) :"},
    { "RNOCUT",    8, 49, 2, 6, 4, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"TOTAL NUMBER OF REFLECTIONS     (NO CUTOFF) :"},
    { "RNOCUT",    7, 49, 2, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
};

REMARKS Remark_3_XPLOR_CNS_Bin[Num_Remark_3_XPLOR_CNS_Bin] = {
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
	"BIN COMPLETENESS (WORKING+TEST) (%) :"},
    { "SHELL",    4, 41, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"REFLECTIONS IN BIN    (WORKING SET) :"},
    { "SHELL",    5, 41, 1, 10, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_R_VALUE, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BIN R VALUE           (WORKING SET) :"},
    { "SHELL",    6, 41, 2, 6, 4, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_R_VALUE, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BIN FREE R VALUE                    :"},
    { "SHELL",   7, 41, 2, 6, 4, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BIN FREE R VALUE TEST SET SIZE  (%) :"},
    { "SHELL",    8, 41, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BIN FREE R VALUE TEST SET COUNT     :"},
    { "SHELL",    9, 41, 1, 5, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"ESTIMATED ERROR OF BIN FREE R VALUE :"},
    { "SHELL",    10, 41, 2, 5, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
};

REMARKS Remark_3_XPLOR_CNS_Bin_Mixed[Num_Remark_3_XPLOR_CNS_Bin_Mixed] = {
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
	"BIN COMPLETENESS (WORKING+TEST) (%) :"},
    { "SHELL",    4, 41, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"REFLECTIONS IN BIN    (WORKING SET) :"},
    { "SHELL",    5, 41, 1, 10, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"REFLECTION IN BIN     (WORKING SET) :"},
    { "SHELL",    5, 41, 1, 10, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_R_VALUE, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BIN R VALUE           (WORKING SET) :"},
    { "SHELL",    6, 41, 2, 6, 4, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_R_VALUE, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BIN FREE R VALUE                    :"},
    { "SHELL",   7, 41, 2, 6, 4, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BIN FREE R VALUE TEST SET SIZE  (%) :"},
    { "SHELL",    8, 41, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BIN FREE R VALUE TEST SET COUNT     :"},
    { "SHELL",    9, 41, 1, 5, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BIN FREE R VALUE SET COUNT          :"},
    { "SHELL",    9, 41, 1, 5, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"ESTIMATED ERROR OF BIN FREE R VALUE :"},
    { "SHELL",    10, 41, 2, 5, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
};

REMARKS Remark_3_XPLOR_CNS_Corss[Num_Remark_3_XPLOR_CNS_Corss] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"CROSS-VALIDATED ESTIMATED COORDINATE ERROR."}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"ESD FROM C-V LUZZATI PLOT    (A) :" },
    { "REFERR",    5, 38, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"ESD FROM C-V SIGMAA          (A) :" },
    { "REFERR",    6, 38, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_XPLOR_CNS_RMSD[Num_Remark_3_XPLOR_CNS_RMSD] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"RMS DEVIATIONS FROM IDEAL VALUES."}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BOND LENGTHS                 (A) :"},
    { "RMSGEN", 2, 38, 2, 5, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BOND ANGLES            (DEGREES) :"},
    { "RMSGEN", 3, 38, 2, 5, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DIHEDRAL ANGLES        (DEGREES) :"},
    { "RMSGEN", 5, 38, 2, 5, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"IMPROPER ANGLES        (DEGREES) :"},
    { "RMSGEN", 6, 38, 2, 5, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_XPLOR_CNS_ITM =
{   3, 2,PDB_REMARK_GENERAL, 0, {
    { "",  0,  2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"ISOTROPIC THERMAL MODEL :" },
    { "RMSGEN",  7, 28, 3, 25, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"" },}
};

REMARKS Remark_3_Model =
{   3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"NCS MODEL :"},
    { "XFILES", 5, 14, 3, 46, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
};

REMARKS Remark_3_XPlor1[Num_Remark_3_XPlor1] = {
  { 3, 1,PDB_REMARK_GENERAL, 0, {
    { "",               0,  2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"NCS RESTRAINTS.                         RMS   SIGMA/WEIGHT"}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 2, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"GROUP"},
    { "NCSMOD",       2,  9, 1, 2, 0, PDB_REMARK_RIGHT_JUSTIFIED,
	""},
    { "",             0,  13, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"POSITIONAL            (A) :"},
    { "NCSMOD", 3, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 47, 3, 0, 0, PDB_REMARK_RIGHT_JUSTIFIED,
	";"},
    { "NCSMOD", 4, 49, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  6, PDB_REMARK_GENERAL, 2, {
    { "",             0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"GROUP"},
    { "NCSMOD",       2,  9, 1, 2, 0, PDB_REMARK_RIGHT_JUSTIFIED,
	""},
    { "",             0,  13, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"B-FACTOR           (A**2) :"},
    { "NCSMOD", 5, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 47, 3, 0, 0, PDB_REMARK_RIGHT_JUSTIFIED,
	";"},
    { "NCSMOD", 6, 49, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_XPlor2[Num_Remark_3_XPlor2] = {
  { 3, 4, PDB_REMARK_GENERAL, 2, {
    { "",         0,  2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"PARAMETER FILE"},
    { "XFILES",   2,  17, 1, 2, 0, PDB_REMARK_RIGHT_JUSTIFIED,
	""},
    { "",         0, 21, 3, 0, 0, PDB_REMARK_RIGHT_JUSTIFIED,
	":"},
    { "XFILES",   3, 23, 3, 25, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 2, {
    { "",         0,  2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"TOPOLOGY FILE"},
    { "XFILES",   2,  16, 1, 2, 0, PDB_REMARK_RIGHT_JUSTIFIED,
	""},
    { "",         0, 21, 3, 0, 0, PDB_REMARK_RIGHT_JUSTIFIED,
	":"},
    { "XFILES",   4, 23, 3, 35, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_XPlor3 = {
    3, 4, PDB_REMARK_GENERAL, 2, {
    { "",         0,  2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"PARAMETER FILE"},
    { "XFILES",   2,  17, 1, 2, 0, PDB_REMARK_RIGHT_JUSTIFIED,
	""},
    { "",         0, 21, 3, 0, 0, PDB_REMARK_RIGHT_JUSTIFIED,
	":"},
    { "XFILES",   3, 23, 3, 25, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
};

REMARKS Remark_3_XPlor4 = {
    3, 4, PDB_REMARK_GENERAL, 2, {
    { "",         0,  2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"TOPOLOGY FILE"},
    { "XFILES",   2,  16, 1, 2, 0, PDB_REMARK_RIGHT_JUSTIFIED,
	""},
    { "",         0, 21, 3, 0, 0, PDB_REMARK_RIGHT_JUSTIFIED,
	":"},
    { "XFILES",   4, 23, 3, 35, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
};

GROUP_REMARKS Remark_XPlor = 
{
  14,
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0 },
  { "", "", "", "", "", "", "", "", "", "", "", "", "", "" },
  { Remark_3_XPLOR, Remark_3_XPLOR_CNS_Fit,
    Remark_3_XPLOR_CNS_Bin, Remark_3_AtomCount,
    Remark_3_BValue, Remark_3_ESD,
    Remark_3_XPLOR_CNS_Corss, Remark_3_XPLOR_CNS_RMSD,
    &Remark_3_XPLOR_CNS_ITM, Remark_3_Isotropic,
    &Remark_3_Model, Remark_3_XPlor1,
    Remark_3_XPlor2, &Remark_3_Other
  },
  { Num_Remark_3_XPLOR, Num_Remark_3_XPLOR_CNS_Fit,
    Num_Remark_3_XPLOR_CNS_Bin, Num_Remark_3_AtomCount,
    Num_Remark_3_BValue, Num_Remark_3_ESD,
    Num_Remark_3_XPLOR_CNS_Corss, Num_Remark_3_XPLOR_CNS_RMSD,
    1, Num_Remark_3_Isotropic, 1, Num_Remark_3_XPlor1,
    Num_Remark_3_XPlor2, 1
  }
};

GROUP_REMARKS Remark_CNS = 
{
  16,
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0},
  { "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" },
  { &Remark_3_CNS_Target, Remark_3_XPLOR, /* Remark_3_CNS, */ Remark_3_XPLOR_CNS_Fit,
    Remark_3_XPLOR_CNS_Bin, Remark_3_AtomCount, Remark_3_BValue,
    Remark_3_ESD, Remark_3_XPLOR_CNS_Corss, Remark_3_XPLOR_CNS_RMSD,
    &Remark_3_XPLOR_CNS_ITM, Remark_3_Isotropic, Remark_3_Solvent,
    &Remark_3_Model, Remark_3_XPlor1, Remark_3_XPlor2, &Remark_3_Other
  },
  { 1, Num_Remark_3_XPLOR, /* Num_Remark_3_CNS, */  Num_Remark_3_XPLOR_CNS_Fit,
    Num_Remark_3_XPLOR_CNS_Bin, Num_Remark_3_AtomCount,
    Num_Remark_3_BValue, Num_Remark_3_ESD,
    Num_Remark_3_XPLOR_CNS_Corss, Num_Remark_3_XPLOR_CNS_RMSD,
    1, Num_Remark_3_Isotropic, Num_Remark_3_Solvent, 1,
    Num_Remark_3_XPlor1, Num_Remark_3_XPlor2, 1
  }
};

GROUP_REMARKS Remark_CNX = 
{
  16,
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /**/ 0,/**/  -1, 0, 0},
  { "", "", "", "", "", "", "", "", "", "", "", "", /**/ "", /**/ "", "", "" },
  { Remark_3_XPLOR, Remark_3_CNX_Fit, Remark_3_CNX_Fit_Agreement,
    Remark_3_XPLOR_CNS_Bin, Remark_3_AtomCount, Remark_3_BValue,
    Remark_3_ESD, Remark_3_XPLOR_CNS_Corss, Remark_3_XPLOR_CNS_RMSD,
    &Remark_3_XPLOR_CNS_ITM, Remark_3_Isotropic, /**/ Remark_3_Solvent, /**/
    &Remark_3_Model, Remark_3_XPlor1, Remark_3_XPlor2, &Remark_3_Other
  },
  { Num_Remark_3_XPLOR, Num_Remark_3_CNX_Fit, Num_Remark_3_CNX_Fit_Agreement, 
    Num_Remark_3_XPLOR_CNS_Bin, Num_Remark_3_AtomCount,
    Num_Remark_3_BValue, Num_Remark_3_ESD,
    Num_Remark_3_XPLOR_CNS_Corss, Num_Remark_3_XPLOR_CNS_RMSD,
    1, Num_Remark_3_Isotropic, /**/ Num_Remark_3_Solvent, /**/ 1,
    Num_Remark_3_XPlor1, Num_Remark_3_XPlor2, 1
  }
};

GROUP_REMARKS Remark_XPlor_out = 
{
  15,
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, 0 },
  { "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" },
  { Remark_3_XPLOR, Remark_3_XPLOR_CNS_Fit,
    Remark_3_XPLOR_CNS_Bin, Remark_3_AtomCount,
    Remark_3_BValue, Remark_3_ESD,
    Remark_3_XPLOR_CNS_Corss, Remark_3_XPLOR_CNS_RMSD,
    &Remark_3_XPLOR_CNS_ITM, Remark_3_Isotropic,
    &Remark_3_Model, Remark_3_XPlor1,
    &Remark_3_XPlor3, &Remark_3_XPlor4, &Remark_3_Other
  },
  { Num_Remark_3_XPLOR, Num_Remark_3_XPLOR_CNS_Fit,
    Num_Remark_3_XPLOR_CNS_Bin, Num_Remark_3_AtomCount,
    Num_Remark_3_BValue, Num_Remark_3_ESD,
    Num_Remark_3_XPLOR_CNS_Corss, Num_Remark_3_XPLOR_CNS_RMSD,
    1, Num_Remark_3_Isotropic, 1, Num_Remark_3_XPlor1, 1, 1, 1
  }
};

GROUP_REMARKS Remark_CNS_out = 
{
  17,
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, 0},
  { "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" },
  { &Remark_3_CNS_Target, Remark_3_XPLOR, /* Remark_3_CNS, */ Remark_3_XPLOR_CNS_Fit,
    Remark_3_XPLOR_CNS_Bin, Remark_3_AtomCount, Remark_3_BValue,
    Remark_3_ESD, Remark_3_XPLOR_CNS_Corss, Remark_3_XPLOR_CNS_RMSD,
    &Remark_3_XPLOR_CNS_ITM, Remark_3_Isotropic, Remark_3_Solvent,
    &Remark_3_Model, Remark_3_XPlor1, &Remark_3_XPlor3, &Remark_3_XPlor4, &Remark_3_Other
  },
  { 1, Num_Remark_3_XPLOR, /* Num_Remark_3_CNS, */  Num_Remark_3_XPLOR_CNS_Fit,
    Num_Remark_3_XPLOR_CNS_Bin, Num_Remark_3_AtomCount,
    Num_Remark_3_BValue, Num_Remark_3_ESD,
    Num_Remark_3_XPLOR_CNS_Corss, Num_Remark_3_XPLOR_CNS_RMSD,
    1, Num_Remark_3_Isotropic, Num_Remark_3_Solvent, 1,
    Num_Remark_3_XPlor1, 1, 1, 1
  }
};

GROUP_REMARKS Remark_CNX_out = 
{
  17,
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /**/ 0,/**/  -1, -1, -1, 0},
  { "", "", "", "", "", "", "", "", "", "", "", "", /**/ "", /**/ "", "", "", "" },
  { Remark_3_XPLOR, Remark_3_CNX_Fit, Remark_3_CNX_Fit_Agreement,
    Remark_3_XPLOR_CNS_Bin, Remark_3_AtomCount, Remark_3_BValue,
    Remark_3_ESD, Remark_3_XPLOR_CNS_Corss, Remark_3_XPLOR_CNS_RMSD,
    &Remark_3_XPLOR_CNS_ITM, Remark_3_Isotropic, /**/ Remark_3_Solvent, /**/
    &Remark_3_Model, Remark_3_XPlor1, &Remark_3_XPlor3, &Remark_3_XPlor4, &Remark_3_Other
  },
  { Num_Remark_3_XPLOR, Num_Remark_3_CNX_Fit, Num_Remark_3_CNX_Fit_Agreement, 
    Num_Remark_3_XPLOR_CNS_Bin, Num_Remark_3_AtomCount,
    Num_Remark_3_BValue, Num_Remark_3_ESD,
    Num_Remark_3_XPLOR_CNS_Corss, Num_Remark_3_XPLOR_CNS_RMSD,
    1, Num_Remark_3_Isotropic, /**/ Num_Remark_3_Solvent, /**/ 1,
    Num_Remark_3_XPlor1, 1, 1, 1
  }
};

GROUP_REMARKS Remark_CNS_Mixed = 
{
  16,
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0},
  { "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" },
  { &Remark_3_CNS_Target, Remark_3_CNS_Mixed, Remark_3_XPLOR_CNS_Fit,
    Remark_3_XPLOR_CNS_Bin_Mixed, Remark_3_AtomCount, Remark_3_BValue,
    Remark_3_ESD, Remark_3_XPLOR_CNS_Corss, Remark_3_XPLOR_CNS_RMSD,
    &Remark_3_XPLOR_CNS_ITM, Remark_3_Isotropic, Remark_3_Solvent,
    &Remark_3_Model, Remark_3_XPlor1, Remark_3_XPlor2, &Remark_3_Other
  },
  { 1, Num_Remark_3_CNS_Mixed, Num_Remark_3_XPLOR_CNS_Fit,
    Num_Remark_3_XPLOR_CNS_Bin_Mixed, Num_Remark_3_AtomCount,
    Num_Remark_3_BValue, Num_Remark_3_ESD,
    Num_Remark_3_XPLOR_CNS_Corss, Num_Remark_3_XPLOR_CNS_RMSD,
    1, Num_Remark_3_Isotropic, Num_Remark_3_Solvent, 1,
    Num_Remark_3_XPlor1, Num_Remark_3_XPlor2, 1
  }
};

REMARKS Remark_3_PRIMEX[Num_Remark_3_PRIMEX] = {
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
    { "RFACTR",    6, 39, 2, 6, 2, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,{
    { "",              0,  3, 3, 0, 0,PDB_REMARK_LEFT_JUSTIFIED,
        "DATA CUTOFF            (SIGMA(F)) :"},
    {"RFACTR",   5,  39, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
       ""},}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "DATA CUTOFF HIGH  (ABS(F)/RMS(F)) :" },
    { "RFACTR",   12, 39, 2, 15, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DATA CUTOFF LOW          (ABS(F)) :" },
    { "RFACTR",   10, 39, 2, 7, 4, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"COMPLETENESS (WORKING+TEST)   (%) :"},
    { "PEROBS",  2,  39, 2, 6, 1, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,  {
    { "",             0,    3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"NUMBER OF REFLECTIONS             :"},
    { "RFACTR",  3,   39, 1, 10, 0,  PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_Solvent_PRIMEX[Num_Remark_3_Solvent_PRIMEX] = {
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
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SOLVENT VOLUME (%) :"},
    { "SOLMOD", 8, 24, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  }
};

REMARKS Remark_3_Parameter_PRIMEX[Num_Remark_3_Parameter_PRIMEX] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "PARAMETERS AND TOPOLOGY."}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "FORCE FIELD       :"},
    { "PAMTOP", 2, 23, 3, 37, 3, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "PLANARITY WEIGHT  :"},
    { "PAMTOP", 3, 23, 3, 37, 3, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SOLVENT MODEL     :"},
    { "PAMTOP", 4, 23, 3, 37, 3, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  }
};

REMARKS Remark_3_Clash_PRIMEX[Num_Remark_3_Clash_PRIMEX] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "CLASH COUNT PER 100 RESIDUES."}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SEVERE :"},
    { "PCLASH", 2, 12, 3, 48, 3, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 3,  2, PDB_REMARK_GENERAL, 0,{
    { "", 0, 3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "OTHER  :"},
    { "PCLASH", 3, 12, 3, 48, 3, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  }
};

GROUP_REMARKS Remark_PRIMEX = 
{
  17,
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0},
  { "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" },
  { &Remark_3_CNS_Target, Remark_3_PRIMEX, Remark_3_XPLOR_CNS_Fit,
    Remark_3_XPLOR_CNS_Bin, Remark_3_AtomCount, Remark_3_BValue,
    Remark_3_Solvent_PRIMEX, Remark_3_ESD, Remark_3_XPLOR_CNS_Corss,
    Remark_3_XPLOR_CNS_RMSD, &Remark_3_XPLOR_CNS_ITM, Remark_3_Isotropic,
    &Remark_3_Model, Remark_3_XPlor1, Remark_3_Parameter_PRIMEX,
    Remark_3_Clash_PRIMEX, &Remark_3_Other
  },
  { 1, Num_Remark_3_PRIMEX, Num_Remark_3_XPLOR_CNS_Fit,
    Num_Remark_3_XPLOR_CNS_Bin, Num_Remark_3_AtomCount,
    Num_Remark_3_BValue, Num_Remark_3_Solvent_PRIMEX, Num_Remark_3_ESD,
    Num_Remark_3_XPLOR_CNS_Corss, Num_Remark_3_XPLOR_CNS_RMSD,
    1, Num_Remark_3_Isotropic, 1, Num_Remark_3_XPlor1,
    Num_Remark_3_Parameter_PRIMEX,  Num_Remark_3_Clash_PRIMEX, 1
  }
};

GROUP_REMARKS *group_remarks[NUM_GROUP_REMARKS] = { &Remark_XPlor_out, &Remark_CNS_out, &Remark_Nuclsq,
                &Remark_Prolsq, &Remark_Refmac, &Remark_Refmac_5_0, &Remark_Tnt,
                &Remark_Buster_Tnt, &Remark_Shelxl, &Remark_CNX_out, &Remark_Phenix_Write, &Remark_PRIMEX
};

GROUP_REMARKS *group_remarks_in[NUM_GROUP_REMARKS] = { &Remark_XPlor, &Remark_CNS_Mixed, &Remark_Nuclsq,
                &Remark_Prolsq_Read, &Remark_Refmac, &Remark_Refmac_5_0_Mixed, &Remark_Tnt,
                &Remark_Buster_Tnt_Read, &Remark_Shelxl, &Remark_CNX, &Remark_Phenix_Read, &Remark_PRIMEX
};
