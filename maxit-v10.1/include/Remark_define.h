/*
FILE:     Remark_define.h
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
#ifndef _H_REMARK_DEFINE_H_
#define _H_REMARK_DEFINE_H_

#define  NAT                               0
#define  GEN                               1
#define  SYN                               2

#define NUM_GROUP_REMARKS                 12

#define PDB_REMARK_GENERAL                 1
#define PDB_REMARK_R_VALUE                 2
#define PDB_REMARK_NNNMMM_SYMOP            3
#define PDB_REMARK_XYZT                    4

#define PDB_REMARK_LEFT_JUSTIFIED          1
#define PDB_REMARK_RIGHT_JUSTIFIED         2

#define Num_Remark_3_1                     3
#define Num_Remark_3_2                     5
#define Num_Remark_3_3                     2
#define Num_Remark_3_NMR                   5
#define Num_Remark_3_General               6
#define Num_Remark_3_Fit_Agreement         7
#define Num_Remark_3_BValue               11
#define Num_Remark_3_BValue_PHENIX        16
#define Num_Remark_3_ESD                   4
#define Num_Remark_3_ESU                   5
#define Num_Remark_3_AtomCount             5
#define Num_Remark_3_Isotropic             5
#define Num_Remark_3_Isotropic_Prolsq      5
#define Num_Remark_3_Isotropic_Mixed       9     
#define Num_Remark_3_Solvent               4
#define Num_Remark_3_XPLOR                 8
#define Num_Remark_3_CNS                   7
#define Num_Remark_3_CNS_Mixed            10
#define Num_Remark_3_CNX_Fit               9
#define Num_Remark_3_CNX_Fit_Agreement     8
#define Num_Remark_3_Cryo_Em              23
#define Num_Remark_3_Cryo_Em_In           29
#define Num_Remark_3_XPLOR_CNS_Fit         8
#define Num_Remark_3_XPLOR_CNS_Bin        11
#define Num_Remark_3_XPLOR_CNS_Bin_Mixed  13
#define Num_Remark_3_XPLOR_CNS_Corss       3
#define Num_Remark_3_XPLOR_CNS_RMSD        5
#define Num_Remark_3_XPlor1                3
#define Num_Remark_3_XPlor2                2
#define Num_Remark_3_PRIMEX                8
#define Num_Remark_3_Solvent_PRIMEX        5
#define Num_Remark_3_Parameter_PRIMEX      4
#define Num_Remark_3_Clash_PRIMEX          3
#define Num_Remark_3_NUCLSQ_PROLSQ         8
#define Num_Remark_3_Prolsq_RMSD           6
#define Num_Remark_3_Prolsq_Plane          2
#define Num_Remark_3_Prolsq_NonBonded      5
#define Num_Remark_3_Prolsq_Torsion        5
#define Num_Remark_3_Nuclsq_RMSD           6
#define Num_Remark_3_Nuclsq_NonBonded      3
#define Num_Remark_3_Nuclsq_Isotropic      5
#define Num_Remark_3_REFMAC_BIN            9
#define Num_Remark_3_REFMAC_BIN_Read      11
#define Num_Remark_3_REFMAC_AtomCount      2
#define Num_Remark_3_REFMAC_CORR           3
#define Num_Remark_3_REFMAC_RMSD          26
#define Num_Remark_3_REFMAC_RMSD_X        78
#define Num_Remark_3_REFMAC_ISOTROPIC     11
#define Num_Remark_3_REFMAC_ISOTROPIC_X   22
#define Num_Remark_3_REFMAC_ANISOTROPIC    4
#define Num_Remark_3_REFMAC_NCS            2
#define Num_Remark_3_REFMAC_NCS_LOCAL      4
#define Num_Remark_3_REFMAC_NCS_GROUP     12
#define Num_Remark_3_REFMAC_TWIN           2
#define Num_Remark_3_REFMAC_TWIN_OP        3
#define Num_Remark_3_REFMAC_TLS            2
#define Num_Remark_3_REFMAC_TLS_GROUP     17
#define Num_Remark_3_Phenix1               6
#define Num_Remark_3_Phenix2               6
#define Num_Remark_3_Phenix_BIN            2
#define Num_Remark_3_REFMAC_BULK           6
#define Num_Remark_3_PHENIX_BULK           6
#define Num_Remark_3_PHENIX_ESU            3
#define Num_Remark_3_PHENIX_RMSD           7
#define Num_Remark_3_PHENIX_TWINNING       3
#define Num_Remark_3_PHENIX_TLS_GROUP     15
#define Num_Remark_3_PHENIX_NCS            2
#define Num_Remark_3_PHENIX_NCS_GROUP_OP   5
#define Num_Remark_3_Shelxl_Refine         7
#define Num_Remark_3_Shelxl_NoCutOff       7
#define Num_Remark_3_Shelxl_4SIG           7
#define Num_Remark_3_Shelxl_Model          6
#define Num_Remark_3_Shelxl_RMSD          11
#define Num_Remark_3_Shelxl_BULK           2
#define Num_Remark_3_Shelxl_STEREO         2
#define Num_Remark_3_General_TNT           8
#define Num_Remark_3_General_Buster_TNT    9
#define Num_Remark_3_Buster_Reference      2
#define Num_Remark_3_Tnt_RMSD              9
#define Num_Remark_3_Tnt_Restraint         3
#define Num_Remark_3_Buster_TNT_Fit_BIN   13
#define Num_Remark_3_Buster_TNT_Fit_BIN_X 18
#define Num_Remark_3_Buster_ESD            6
#define Num_Remark_3_Buster_Tnt_RMSD      17
#define Num_Remark_3_Buster_Tnt_RMSD1      5
#define Num_Remark_3_BUSTER_TLS_GROUP     15
#define Num_Remark_3_BUSTER_TLS_GROUP_Read 17
#define Num_Remark_100                     2
#define Num_Remark_100_NoDate              2
#define Num_Remark_100_RCSB                2
#define Num_Remark_100_BMRB                2
#define Num_Remark_100_EBI                 2
#define Num_Remark_100_model               2
#define Num_Remark_200                    48
#define Num_Remark_205                     5
#define Num_Remark_210                    24
#define Num_Remark_215                     5
#define Num_Remark_217                     5
#define Num_Remark_220                     4
#define Num_Remark_225                     5
#define Num_Remark_230                    44
#define Num_Remark_240                    27
#define Num_Remark_245                    34
#define Num_Remark_247                     6
#define Num_Remark_250                     5
#define Num_Remark_265_2                  21
#define Num_Remark_265_3                   5
#define Num_Remark_265_4                   3
#define Num_Remark_280                     5
#define Num_Remark_290                    17
#define Num_Remark_300                     5
#define Num_Remark_300_Cryo_Em_A           6
#define Num_Remark_300_Cryo_Em_V           8
#define Num_Remark_350                     6
#define Num_Remark_375                     6
#define Num_Remark_465                     6
#define Num_Remark_465_NMR                 4
#define Num_Remark_470                     5
#define Num_Remark_470_NMR                 3
#define Num_Remark_475                     6
#define Num_Remark_480                     6
#define Num_Remark_500_A                   6
#define Num_Remark_500_S                  15
#define Num_Remark_500_NMR                 6
#define Num_Remark_500_BOND               15
#define Num_Remark_500_ANGLE              15
#define Num_Remark_500_TORSION            14
#define Num_Remark_500_NON_CIS             8
#define Num_Remark_500_SPLANE             12
#define Num_Remark_500_MPLANE             10
#define Num_Remark_500_CHIRAL             14
#define Num_Remark_525                    11
#define Num_Remark_610                     5
#define Num_Remark_615                     6
#define Num_Remark_620                     3
#define Num_Remark_800                     4
#define NUM_200_ALIGN                      6

#define NUM_METHOD  2

#endif
