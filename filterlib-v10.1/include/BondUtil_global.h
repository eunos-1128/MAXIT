/*
FILE:     BondUtil_global.h
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
#ifndef _H_ALLOWED_BOND_H_
#define _H_ALLOWED_BOND_H_

#include <map>
#include <string>

#define NUM_BOND  249

static _ALLOWED_BOND _allowed_bonds[NUM_BOND] = {
       { "AG_N",  "N",  "Ag", 1.85, 2.70, 2, false },
       { "AG_O",  "O",  "Ag", 1.85, 2.70, 2, false },
       { "AG_S",  "S",  "Ag", 2.00, 3.00, 2, false },
       { "AL_F",  "Al", "F",  1.70, 2.80, 3, false },
       { "AL_H",  "Al", "H",  1.35, 1.65, 1, false },
       { "AL_O",  "Al", "O",  1.70, 2.80, 3, false },
       { "AS_C",  "As", "C",  1.70, 2.80, 3, false },
       { "AS_H",  "As", "H",  1.20, 1.60, 1, false },
       { "AS_O",  "As", "O",  1.70, 2.80, 3, false },
       { "AS_PT", "As", "Pt", 1.70, 2.80, 3, false },
       { "AS_S",  "As", "S",  1.70, 2.80, 3, false },
       { "AU_C",  "Au", "C",  1.70, 2.80, 3, false },
       { "AU_CL", "Au", "Cl", 1.70, 2.80, 3, false },
       { "AU_N",  "N",  "Au", 1.80, 2.80, 2, false },
       { "AU_O",  "O",  "Au", 1.80, 2.80, 2, false },
       { "AU_P",  "Au", "P",  1.70, 2.80, 3, false },
       { "AU_S",  "S",  "Au", 1.80, 3.00, 2, false },
       { "BA_O",  "O",  "Ba", 2.50, 3.50, 2, false },
       { "BE_F",  "Be", "F",  1.70, 2.80, 3, false },
       { "BE_N",  "N",  "Be", 1.40, 2.20, 2, false },
       { "BE_O",  "O",  "Be", 1.40, 2.20, 2, false },
       { "BR_C",  "C",  "Br", 1.70, 2.15, 1, true  },
       { "BR_HG", "Br", "Hg", 1.70, 2.80, 3, false },
       { "BR_N",  "N",  "Br", 1.70, 2.15, 1, true  },
       { "BR_O",  "O",  "Br", 1.50, 2.00, 1, true  },
       { "BR_PB", "Br", "Pb", 1.70, 2.80, 3, false },
       { "BR_PD", "Br", "Pd", 1.70, 2.80, 3, false },
       { "BR_PT", "Br", "Pt", 1.70, 2.80, 3, false },
       { "BR_SN", "Br", "Sn", 1.70, 2.80, 2, false }, // added by DAOTHER-5962
       { "BR_TA", "Br", "Ta", 1.70, 2.80, 3, false },
       { "B_B",   "B",  "B",  1.45, 1.95, 1, true  },
       { "B_C",   "C",  "B",  1.20, 1.85, 1, false },
       { "B_CO",  "B",  "Co", 1.70, 2.80, 3, false },
       { "B_F",   "F",  "B",  1.20, 1.75, 1, false },
       { "B_H",   "B",  "H",  0.85, 1.20, 1, false }, // new added
       { "B_N",   "N",  "B",  1.20, 1.65, 1, false },
       { "B_O",   "O",  "B",  1.20, 1.65, 1, false },
       { "B_P",   "B",  "P",  1.50, 1.95, 1, false }, // new added
       { "B_S",   "B",  "S",  1.55, 2.00, 1, true  }, // new added
       { "CA_O",  "O",  "Ca", 1.85, 3.20, 2, false },
       { "CD_N",  "N",  "Cd", 1.75, 2.70, 2, false },
       { "CD_O",  "O",  "Cd", 1.75, 2.70, 2, false },
       { "CD_S",  "S",  "Cd", 1.90, 3.10, 2, false },
       { "CE_N",  "N",  "Ce", 2.00, 3.50, 2, false },
       { "CE_O",  "O",  "Ce", 2.00, 3.50, 2, false },
       { "CL_CU", "Cl", "Cu", 2.10, 2.80, 2, false },
       { "CL_FE", "Cl", "Fe", 1.70, 2.80, 3, false },
       { "CL_HG", "Cl", "Hg", 1.70, 2.80, 3, false },
       { "CL_IR", "Cl", "Ir", 1.70, 2.80, 3, false },
       { "CL_MN", "Mn", "Cl", 1.90, 2.80, 2, false },
       { "CL_N",  "N",  "Cl", 1.50, 1.95, 1, true  },
       { "CL_NI", "Cl", "Ni", 1.90, 2.70, 2, false },
       { "CL_O",  "O",  "Cl", 1.40, 1.95, 1, true  },
       { "CL_OS", "Cl", "Os", 1.70, 2.80, 3, false },
       { "CL_P",  "Cl", "P",  1.60, 2.10, 1, false }, // new added
       { "CL_PD", "Cl", "Pd", 1.70, 2.80, 3, false },
       { "CL_PT", "Cl", "Pt", 1.70, 2.80, 3, false },
       { "CL_RE", "Cl", "Re", 1.90, 2.70, 2, false },
       { "CL_RH", "Cl", "Rh", 1.70, 2.80, 3, false },
       { "CL_RU", "Cl", "Ru", 1.70, 2.80, 3, false },
       { "CL_S",  "Cl", "S",  1.55, 2.40, 1, true  }, // new added 
       { "CL_SN", "Cl", "Sn", 1.70, 2.80, 2, false }, // added by DAOTHER-5962
       { "CO_N",  "N",  "Co", 1.70, 2.80, 2, false },
       { "CO_O",  "O",  "Co", 1.70, 2.80, 2, false },
       { "CR_N",  "N",  "Cr", 1.70, 2.80, 2, false },
       { "CR_O",  "O",  "Cr", 1.70, 2.80, 2, false },
       { "CS_N",  "N",  "Cs", 2.30, 3.30, 2, false },
       { "CS_O",  "O",  "Cs", 2.30, 3.50, 2, false },
       { "CU_CU", "Cu", "Cu", 1.80, 2.80, 2, false },
       { "CU_N",  "N",  "Cu", 1.70, 2.70, 2, false },
       { "CU_O",  "O",  "Cu", 1.70, 2.70, 2, false },
       { "CU_S",  "S",  "Cu", 1.90, 2.70, 2, false },
       { "C_C",   "C",  "C",  1.05, 1.65, 1, false },
       { "C_CL",  "C",  "Cl", 1.50, 1.95, 1, true  },
       { "C_CO",  "C",  "Co", 1.70, 2.80, 3, false }, // removed by DAOTHER-5086
       { "C_F",   "C",  "F",  1.20, 1.65, 1, false },
       { "C_FE",  "C",  "Fe", 1.70, 2.80, 3, false },
       { "C_H",   "C",  "H",  0.85, 1.15, 1, false },
       { "C_HG",  "C",  "Hg", 1.70, 2.80, 3, false },
       { "C_I",   "C",  "I",  1.95, 2.30, 1, true  },
       { "C_IR",  "C",  "Ir", 1.70, 2.80, 3, false },
       { "C_MN",  "C",  "Mn", 1.70, 2.80, 3, false },
       { "C_MO",  "C",  "Mo", 1.70, 2.80, 3, false },
       { "C_N",   "C",  "N",  1.20, 1.60, 1, false },
       { "C_NI",  "C",  "Ni", 1.70, 2.80, 3, false },
       { "C_O",   "C",  "O",  1.10, 1.55, 1, false },
       { "C_OS",  "C",  "Os", 1.70, 2.80, 3, false },
       { "C_P",   "C",  "P",  1.65, 1.87, 1, true  },
       { "C_PB",  "C",  "Pb", 1.70, 2.80, 3, false },
       { "C_PD",  "C",  "Pd", 1.70, 2.80, 3, false },
       { "C_PT",  "C",  "Pt", 1.70, 2.80, 3, false },
       { "C_RE",  "C",  "Re", 1.80, 2.80, 3, false }, // removed by DAOTHER-5086
       { "C_RH",  "C",  "Rh", 1.70, 2.80, 3, false },
       { "C_RU",  "C",  "Ru", 1.70, 2.80, 3, false },
       { "C_S",   "C",  "S",  1.55, 2.00, 1, true  },
       { "C_SB",  "C",  "Sb", 1.70, 2.80, 3, false },
       { "C_SE",  "C",  "Se", 1.75, 2.15, 1, true  },
       { "C_SI",  "C",  "Si", 1.70, 2.80, 3, false },
       { "C_SN",  "C",  "Sn", 1.70, 2.80, 2, false }, // added by DAOTHER-5962
       { "C_TE",  "C",  "Te", 1.70, 2.80, 3, false },
       { "DY_N",  "N",  "Dy", 2.00, 3.50, 2, false },
       { "DY_O",  "O",  "Dy", 2.00, 3.50, 2, false },
       { "ER_N",  "N",  "Er", 2.00, 3.50, 2, false },
       { "ER_O",  "O",  "Er", 2.00, 3.50, 2, false },
       { "EU_N",  "N",  "Eu", 2.00, 3.50, 2, false },
       { "EU_O",  "O",  "Eu", 2.00, 3.50, 2, false },
       { "FE_FE", "Fe", "Fe", 1.70, 2.80, 3, false },
       { "FE_H",  "Fe", "H",  1.20, 1.65, 3, false },
       { "FE_N",  "N",  "Fe", 1.70, 2.80, 2, false },
       { "FE_NI", "Fe", "Ni", 1.70, 2.80, 3, false },
       { "FE_O",  "O",  "Fe", 1.70, 2.80, 2, false },
       { "FE_P",  "P",  "Fe", 1.90, 2.70, 2, false },
       { "FE_S",  "S",  "Fe", 1.90, 2.70, 2, false },
       { "FE_SE", "Se", "Fe", 2.10, 2.70, 2, false },
       { "F_MG",  "F",  "Mg", 1.70, 2.80, 3, false },
       { "F_MN",  "Mn", "F",  1.60, 2.50, 2, false },
       { "F_N",   "N",  "F",  1.10, 1.60, 1, false },
       { "F_P",   "P",  "F",  1.25, 1.75, 1, false },
       { "F_S",   "F",  "S",  1.55, 2.00, 1, true  }, // new added 
       { "F_SC",  "F",  "Sc", 1.70, 2.80, 3, false },
       { "GA_N",  "N",  "Ga", 1.70, 2.30, 2, false },
       { "GA_O",  "O",  "Ga", 1.70, 2.30, 2, false },
       { "GA_S",  "S",  "Ga", 2.00, 2.60, 2, false },
       { "GD_N",  "N",  "Gd", 2.00, 3.50, 2, false },
       { "GD_O",  "O",  "Gd", 2.00, 3.50, 2, false },
       { "HF_O",  "Hf", "O",  1.70, 2.80, 3, false },
       { "HG_HG", "Hg", "Hg", 1.70, 2.80, 3, false },
       { "HG_I",  "Hg", "I",  1.70, 2.80, 3, false },
       { "HG_N",  "N",  "Hg", 2.00, 2.80, 2, false },
       { "HG_O",  "O",  "Hg", 2.00, 3.20, 2, false },
       { "HG_S",  "S",  "Hg", 2.30, 3.30, 2, false },
       { "HO_N",  "N",  "Ho", 2.00, 3.50, 2, false },
       { "HO_O",  "O",  "Ho", 2.00, 3.50, 2, false },
       { "H_MO",  "H",  "Mo", 1.15, 1.65, 3, false },
       { "H_N",   "N",  "H",  0.85, 1.15, 1, false },
       { "H_NI",  "H",  "Ni", 1.15, 1.65, 3, false },
       { "H_O",   "O",  "H",  0.85, 1.15, 1, false },
       { "H_P",   "P",  "H",  1.10, 1.45, 1, false },
       { "H_S",   "S",  "H",  1.10, 1.45, 1, false },
       { "H_SE",  "H",  "Se", 1.20, 1.65, 1, false }, // new added
       { "H_SI",  "Si", "H",  1.20, 1.60, 1, false },
       { "IN_N",  "N",  "In", 2.00, 2.70, 2, false },
       { "IN_O",  "O",  "In", 2.00, 2.70, 2, false },
       { "IR_N",  "N",  "Ir", 1.70, 2.40, 2, false },
       { "IR_O",  "O",  "Ir", 1.70, 2.40, 2, false },
       { "I_I",   "I",  "I",  2.70, 3.10, 1, true  },
       { "I_N",   "N",  "I",  1.95, 2.30, 1, true  },
       { "I_O",   "O",  "I",  1.95, 2.30, 1, true  },
       { "I_PT",  "I",  "Pt", 1.70, 2.80, 3, false },
       { "I_S",   "S",  "I",  2.70, 3.10, 1, true  },
       { "K_O",   "O",  "K",  2.00, 3.50, 2, false },
       { "LA_N",  "N",  "La", 2.00, 3.50, 2, false },
       { "LA_O",  "O",  "La", 2.00, 3.50, 2, false },
       { "LI_O",  "O",  "Li", 1.70, 2.60, 2, false },
       { "LU_N",  "N",  "Lu", 2.00, 3.50, 2, false },
       { "LU_O",  "O",  "Lu", 2.00, 3.50, 2, false },
       { "MG_N",  "Mg", "N",  1.70, 3.00, 3, false },
       { "MG_O",  "O",  "Mg", 1.70, 3.00, 2, false },
       { "MN_MN", "Mn", "Mn", 2.80, 3.00, 2, false },
       { "MN_N",  "Mn", "N",  1.70, 2.80, 2, false },
       { "MN_O",  "Mn", "O",  1.70, 2.80, 2, false },
       { "MN_P",  "Mn", "P",  1.90, 2.60, 2, false },
       { "MN_S",  "Mn", "S",  2.00, 2.90, 2, false },
       { "MN_SE", "Mn", "Se", 2.20, 2.90, 2, false },
       { "MO_N",  "Mo", "N",  1.70, 2.80, 3, false },
       { "MO_O",  "Mo", "O",  1.70, 2.80, 3, false },
       { "MO_P",  "Mo", "P",  2.10, 2.80, 2, false },
       { "MO_S",  "Mo", "S",  2.00, 2.90, 2, false },
       { "MO_SE", "Mo", "Se", 2.10, 2.90, 2, false },
       { "NA_O",  "O",  "Na", 1.90, 3.20, 2, false },
       { "NI_O",  "O",  "Ni", 1.70, 2.80, 2, false },
       { "NI_S",  "S",  "Ni", 1.90, 2.70, 2, false },
       { "N_N",   "N",  "N",  1.15, 1.55, 1, false },
       { "N_NI",  "N",  "Ni", 1.70, 2.80, 2, false },
       { "N_O",   "N",  "O",  1.18, 1.50, 1, false },
       { "N_OS",  "N",  "Os", 1.70, 2.80, 2, false },
       { "N_P",   "N",  "P",  1.54, 1.77, 1, true  },
       { "N_PB",  "N",  "Pb", 1.90, 3.30, 2, false },
       { "N_PD",  "N",  "Pd", 1.80, 2.80, 2, false },
       { "N_PR",  "N",  "Pr", 2.00, 3.50, 2, false },
       { "N_PT",  "N",  "Pt", 1.80, 2.80, 2, false },
       { "N_RB",  "N",  "Rb", 2.40, 3.20, 2, false },
       { "N_RE",  "N",  "Re", 1.80, 2.80, 2, false },
       { "N_RH",  "N",  "Rh", 1.70, 2.80, 2, false },
       { "N_RU",  "N",  "Ru", 1.70, 2.80, 2, false },
       { "N_S",   "N",  "S",  1.45, 1.90, 1, true  },
       { "N_SM",  "N",  "Sm", 2.00, 3.50, 2, false },
       { "N_SN",  "N",  "Sn", 1.70, 2.80, 2, false }, // added by DAOTHER-5962
       { "N_SR",  "N",  "Sr", 2.45, 3.00, 2, false },
       { "N_TB",  "N",  "Tb", 2.00, 3.50, 2, false },
       { "N_TL",  "N",  "Tl", 2.20, 3.60, 2, false },
       { "N_U",   "N",  "U",  1.90, 3.90, 2, false },
       { "N_V",   "N",  "V",  1.70, 2.80, 2, false },
       { "N_W",   "N",  "W",  2.00, 3.30, 2, false },
       { "N_Y",   "N",  "Y",  2.00, 3.50, 2, false },
       { "N_YB",  "N",  "Yb", 2.00, 3.50, 2, false },
       { "N_ZN",  "N",  "Zn", 1.70, 2.70, 2, false },
       { "OS_S",  "Os", "S",  1.70, 2.80, 3, false },
       { "O_O",   "O",  "O",  1.30, 1.60, 1, false },
       { "O_P",   "O",  "P",  1.40, 1.75, 1, false },
       { "O_OS",  "O",  "Os", 1.70, 2.80, 2, false },
       { "O_PB",  "O",  "Pb", 1.90, 3.30, 2, false },
       { "O_PD",  "O",  "Pd", 1.80, 2.80, 2, false },
       { "O_PR",  "O",  "Pr", 2.00, 3.50, 2, false },
       { "O_PT",  "O",  "Pt", 1.80, 2.80, 2, false },
       { "O_RB",  "O",  "Rb", 2.40, 3.20, 2, false },
       { "O_RE",  "O",  "Re", 1.80, 2.80, 2, false },
       { "O_RH",  "O",  "Rh", 1.70, 2.80, 2, false },
       { "O_RU",  "O",  "Ru", 1.70, 2.80, 2, false },
       { "O_S",   "O",  "S",  1.35, 1.80, 1, false },
       { "O_SB",  "O",  "Sb", 1.70, 2.80, 3, false },
       { "O_SE",  "O",  "Se", 1.70, 2.80, 3, false },
       { "O_SI",  "O",  "Si", 1.70, 2.80, 3, false },
       { "O_SM",  "O",  "Sm", 2.00, 3.50, 2, false },
       { "O_SN",  "O",  "Sn", 1.70, 2.80, 2, false }, // added by DAOTHER-5962
       { "O_SR",  "O",  "Sr", 2.30, 3.00, 2, false },
       { "O_TB",  "O",  "Tb", 2.00, 3.50, 2, false },
       { "O_TE",  "O",  "Te", 1.80, 2.50, 2, false },
       { "O_TL",  "O",  "Tl", 2.20, 3.60, 2, false },
       { "O_U",   "O",  "U",  1.90, 3.90, 2, false },
       { "O_V",   "O",  "V",  1.70, 2.80, 2, false },
       { "O_W",   "O",  "W",  2.00, 3.30, 2, false },
       { "O_Y",   "O",  "Y",  2.00, 3.50, 2, false },
       { "O_YB",  "O",  "Yb", 2.00, 3.50, 2, false },
       { "O_ZN",  "O",  "Zn", 1.70, 2.70, 2, false },
       { "O_ZR",  "O",  "Zr", 1.70, 2.80, 3, false },
       { "PB_S",  "Pb", "S",  1.70, 2.80, 3, false },
       { "PD_S",  "S",  "Pd", 2.00, 2.80, 2, false },
       { "PT_S",  "S",  "Pt", 2.00, 2.90, 2, false },
       { "P_P",   "P",  "P",  1.95, 2.35, 1, true  },
       { "P_PD",  "P",  "Pd", 1.70, 2.80, 3, false },
       { "P_RE",  "P",  "Re", 2.00, 2.80, 2, false },
       { "P_PT",  "P",  "Pt", 1.70, 2.80, 3, false },
       { "P_RH",  "P",  "Rh", 1.70, 2.80, 3, false },
       { "P_RU",  "P",  "Ru", 1.70, 2.80, 3, false },
       { "P_S",   "P",  "S",  1.50, 2.00, 1, true  }, // new added
       { "P_SE",  "P",  "Se", 1.70, 2.80, 3, false },
       { "RE_S",  "S",  "Re", 2.10, 2.80, 2, false },
       { "RH_S",  "Rh", "S",  1.70, 2.80, 3, false },
       { "RU_RU", "Ru", "Ru", 1.70, 2.80, 3, false },
       { "RU_S",  "S",  "Ru", 2.20, 3.20, 2, false },
       { "SE_SE", "Se", "Se", 2.05, 2.55, 1, true  },
       { "S_SE",  "S",  "Se", 1.95, 2.40, 1, true  },
       { "S_S",   "S",  "S",  1.85, 2.20, 1, true  },
       { "S_V",   "S",  "V",  1.70, 2.80, 3, false },
       { "S_W",   "S",  "W",  2.00, 2.80, 2, false },
       { "S_ZN",  "S",  "Zn", 1.90, 3.00, 2, false },
       { "TA_TA", "Ta", "Ta", 1.70, 3.30, 3, false },
       { "V_V",   "V",  "V",  1.70, 3.30, 3, false }
};

static std::map<std::string, int> _bond_mapping;

#endif
