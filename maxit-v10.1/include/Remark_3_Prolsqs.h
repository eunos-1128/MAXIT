/*
FILE:     Remark_3_Prolsqs.h
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
REMARKS Remark_3_NUCLSQ_PROLSQ[Num_Remark_3_NUCLSQ_PROLSQ] = {
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
};

REMARKS Remark_3_ESU[Num_Remark_3_ESU] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"ESTIMATED OVERALL COORDINATE ERROR."}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"ESU BASED ON R VALUE                            (A):" },
    { "REFESU", 2, 56, 2, 7, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,

	"ESU BASED ON FREE R VALUE                       (A):" },
    { "REFESU", 3, 56, 2, 7, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"ESU BASED ON MAXIMUM LIKELIHOOD                 (A):" },
    { "REFESU", 4, 56, 2, 7, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 2, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"ESU FOR B VALUES BASED ON MAXIMUM LIKELIHOOD (A**2):" },
    { "REFESU", 5, 56, 2, 7, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_Nuclsq_RMSD[Num_Remark_3_Nuclsq_RMSD] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"RMS DEVIATIONS FROM IDEAL VALUES."}}
  },
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DISTANCE RESTRAINTS.                    RMS     SIGMA" },}
  },
  { 3,  4, PDB_REMARK_GENERAL, 0,{
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SUGAR-BASE BOND DISTANCE        (A) :"},
    { "RMSDIN", 2, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSDIN", 3, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SUGAR-BASE BOND ANGLE DISTANCE  (A) :"},
    { "RMSDIN", 4, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSDIN", 5, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"PHOSPHATE BONDS DISTANCE        (A) :"},
    { "RMSDIN", 6, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSDIN", 7, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"PHOSPHATE BOND ANGLE, H-BOND    (A) :"},
    { "RMSDIN", 8, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSDIN", 9, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_Nuclsq_NonBonded[Num_Remark_3_Nuclsq_NonBonded] = {
  { 3, 1,PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"NON-BONDED CONTACT RESTRAINTS."}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SINGLE TORSION CONTACT          (A) :"},
    { "RMSNBO", 2, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSNBO", 3, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"MULTIPLE TORSION CONTACT        (A) :"},
    { "RMSNBO", 4, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSNBO", 5, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_Nuclsq_Isotropic[Num_Remark_3_Nuclsq_Isotropic] = {
  { 3, 1,PDB_REMARK_GENERAL, 0, {
    { "",               0,  2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"ISOTROPIC THERMAL FACTOR RESTRAINTS.    RMS    SIGMA"}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SUGAR-BASE BONDS             (A**2) :"},
    { "RMSISN", 2, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 47, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSISN", 3, 49, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SUGAR-BASE ANGLES            (A**2) :"},
    { "RMSISN", 4, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 47, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSISN", 5, 49, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"PHOSPHATE BONDS              (A**2) :"},
    { "RMSISN", 6, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 47, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSISN", 7, 49, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"PHOSPHATE BOND ANGLE, H-BOND (A**2) :"},
    { "RMSISN", 8, 41, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 47, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSISN", 9, 49, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_Prolsq_RMSD[Num_Remark_3_Prolsq_RMSD] = {
  { 3,  1, PDB_REMARK_GENERAL, 0,{
    { "", 0, 2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"RMS DEVIATIONS FROM IDEAL VALUES."}}
  },
  { 3, 1, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DISTANCE RESTRAINTS.                    RMS    SIGMA" },}
  },
  { 3,  4, PDB_REMARK_GENERAL, 0,{
    { "", 0, 4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BOND LENGTH                     (A) :"},
    { "RMSDIS", 2, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSDIS", 3, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"ANGLE DISTANCE                  (A) :"},
    { "RMSDIS", 4, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSDIS", 5, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"INTRAPLANAR 1-4 DISTANCE        (A) :"},
    { "RMSDIS", 6, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSDIS", 7, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"H-BOND OR METAL COORDINATION    (A) :"},
    { "RMSDIS", 8, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSDIS", 9, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_Prolsq_Plane[Num_Remark_3_Prolsq_Plane] = {
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"PLANE RESTRAINT                  (A) :"},
    { "RMSPLA", 2, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSPLA", 3, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"CHIRAL-CENTER RESTRAINT       (A**3) :"},
    { "RMSVOL", 2, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSVOL", 3, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_Prolsq_NonBonded[Num_Remark_3_Prolsq_NonBonded] = {
  { 3, 1,PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"NON-BONDED CONTACT RESTRAINTS."}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SINGLE TORSION                  (A) :"},
    { "RMSNBO", 2, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSNBO", 3, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"MULTIPLE TORSION                (A) :"},
    { "RMSNBO", 4, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSNBO", 5, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"H-BOND (X...Y)                  (A) :"},
    { "RMSNBO", 8, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSNBO", 9, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"H-BOND (X-H...Y)                (A) :"},
    { "RMSNBO", 6, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSNBO", 7, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_Prolsq_Torsion[Num_Remark_3_Prolsq_Torsion] = {
  { 3, 1,PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"CONFORMATIONAL TORSION ANGLE RESTRAINTS."}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SPECIFIED                 (DEGREES) :"},
    { "RMSCON", 10, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSCON", 11, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"PLANAR                    (DEGREES) :"},
    { "RMSCON", 2, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSCON", 3, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"STAGGERED                 (DEGREES) :"},
    { "RMSCON", 4, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSCON", 5, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  4, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"TRANSVERSE                (DEGREES) :"},
    { "RMSCON", 8, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSCON", 9, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_Isotropic_Mixed[Num_Remark_3_Isotropic_Mixed] = {
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
	"MAIN-CHAIN BOND               (A**2) :"},
    { "RMSISO", 2, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSISO", 3, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
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
	"MAIN-CHAIN ANGLE              (A**2) :"},
    { "RMSISO", 4, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSISO", 5, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
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
	"SIDE-CHAIN BOND               (A**2) :"},
    { "RMSISO", 6, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSISO", 7, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
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
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SIDE-CHAIN ANGLE              (A**2) :"},
    { "RMSISO", 8, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSISO", 9, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

REMARKS Remark_3_Isotropic_Prolsq[Num_Remark_3_Isotropic_Prolsq] = {
  { 3, 1,PDB_REMARK_GENERAL, 0, {
    { "",               0,  2, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"ISOTROPIC THERMAL FACTOR RESTRAINTS.    RMS    SIGMA"}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"MAIN-CHAIN BOND               (A**2) :"},
    { "RMSISO", 2, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSISO", 3, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"MAIN-CHAIN ANGLE              (A**2) :"},
    { "RMSISO", 4, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSISO", 5, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SIDE-CHAIN BOND               (A**2) :"},
    { "RMSISO", 6, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSISO", 7, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 3, 4, PDB_REMARK_GENERAL, 0, {
    { "",               0,  3, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SIDE-CHAIN ANGLE              (A**2) :"},
    { "RMSISO", 8, 42, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""},
    { "", 0, 48, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	";"},
    { "RMSISO", 9, 50, 2, 6, 3, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }
};

GROUP_REMARKS Remark_Nuclsq = 
{
  11,
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { "", "", "", "", "", "", "", "", "", "", "" },
  { Remark_3_General,  Remark_3_NUCLSQ_PROLSQ, 
    Remark_3_Fit_Agreement, Remark_3_AtomCount,
    Remark_3_BValue, Remark_3_ESD, Remark_3_Nuclsq_RMSD,
    Remark_3_Prolsq_Plane, Remark_3_Nuclsq_NonBonded,
    Remark_3_Nuclsq_Isotropic, &Remark_3_Other
  },
  { Num_Remark_3_General, Num_Remark_3_NUCLSQ_PROLSQ, 
    Num_Remark_3_Fit_Agreement, Num_Remark_3_AtomCount,
    Num_Remark_3_BValue, Num_Remark_3_ESD, Num_Remark_3_Nuclsq_RMSD,
    Num_Remark_3_Prolsq_Plane, Num_Remark_3_Nuclsq_NonBonded,
    Num_Remark_3_Nuclsq_Isotropic, 1
  }
};

GROUP_REMARKS Remark_Prolsq_Read = 
{
  12,
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { "", "", "", "", "", "", "", "", "", "", "", "" },
  { Remark_3_General,  Remark_3_NUCLSQ_PROLSQ, 
    Remark_3_Fit_Agreement, Remark_3_AtomCount,
    Remark_3_BValue, Remark_3_ESD,
    Remark_3_Prolsq_RMSD, Remark_3_Prolsq_Plane,
    Remark_3_Prolsq_NonBonded, Remark_3_Prolsq_Torsion,
    Remark_3_Isotropic_Mixed, &Remark_3_Other
  },
  { Num_Remark_3_General, Num_Remark_3_NUCLSQ_PROLSQ, 
    Num_Remark_3_Fit_Agreement, Num_Remark_3_AtomCount,
    Num_Remark_3_BValue, Num_Remark_3_ESD,
    Num_Remark_3_Prolsq_RMSD, Num_Remark_3_Prolsq_Plane,
    Num_Remark_3_Prolsq_NonBonded, Num_Remark_3_Prolsq_Torsion,
    Num_Remark_3_Isotropic_Mixed, 1
  }
};

GROUP_REMARKS Remark_Prolsq = 
{
  12,
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { "", "", "", "", "", "", "", "", "", "", "", "" },
  { Remark_3_General,  Remark_3_NUCLSQ_PROLSQ, 
    Remark_3_Fit_Agreement, Remark_3_AtomCount,
    Remark_3_BValue, Remark_3_ESD,
    Remark_3_Prolsq_RMSD, Remark_3_Prolsq_Plane,
    Remark_3_Prolsq_NonBonded, Remark_3_Prolsq_Torsion,
    Remark_3_Isotropic_Prolsq, &Remark_3_Other
  },
  { Num_Remark_3_General, Num_Remark_3_NUCLSQ_PROLSQ, 
    Num_Remark_3_Fit_Agreement, Num_Remark_3_AtomCount,
    Num_Remark_3_BValue, Num_Remark_3_ESD,
    Num_Remark_3_Prolsq_RMSD, Num_Remark_3_Prolsq_Plane,
    Num_Remark_3_Prolsq_NonBonded, Num_Remark_3_Prolsq_Torsion,
    Num_Remark_3_Isotropic_Prolsq, 1
  }
};

GROUP_REMARKS Remark_Refmac = 
{
  11,
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { "", "", "", "", "", "", "", "", "", "", "" },
  { Remark_3_General,  Remark_3_NUCLSQ_PROLSQ, 
    Remark_3_AtomCount, Remark_3_BValue, Remark_3_ESU,
    Remark_3_Prolsq_RMSD, Remark_3_Prolsq_Plane,
    Remark_3_Prolsq_NonBonded, Remark_3_Prolsq_Torsion,
    Remark_3_Isotropic, &Remark_3_Other
  },
  { Num_Remark_3_General, Num_Remark_3_NUCLSQ_PROLSQ, 
    Num_Remark_3_AtomCount, Num_Remark_3_BValue, Num_Remark_3_ESU,
    Num_Remark_3_Prolsq_RMSD, Num_Remark_3_Prolsq_Plane,
    Num_Remark_3_Prolsq_NonBonded, Num_Remark_3_Prolsq_Torsion,
    Num_Remark_3_Isotropic, 1
  }
};
