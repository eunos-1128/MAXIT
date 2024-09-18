/*
FILE:     Remark_200s.h
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

Remark_200_Format _remark_200_align[NUM_200_ALIGN] = {
       { "DIFPTL", 3, 1, { 2 }          },
       { "DTMEAS", 6, 2, { 1, 3 }       },
       { "DTMETH", 6, 4, { 1, 2, 3, 4 } },
       { "DTTEMP", 4, 1, { 1 }          },
       { "RADIAT", 1, 4, { 2, 3, 4, 5 } },
       { "WAVLEN", 1, 1, { 2 }          }
};

REMARKS  Remark_200[Num_Remark_200] =
{
 { 200, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"EXPERIMENTAL DETAILS"},}
 },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"EXPERIMENT TYPE                :"},
    { "EXPDTA",  3,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DATE OF DATA COLLECTION        :"},
    { "DTMEAS",  2,  35, 3, 11,0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"TEMPERATURE           (KELVIN) :"},
    { "DTTEMP",  2,  35, 3, 24, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"PH                             :"},
    { "PHVAL",  2,  35, 3, 24, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"NUMBER OF CRYSTALS USED        :"},
    { "CRMETH",  4,  35, 1, 10, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
 },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SYNCHROTRON              (Y/N) :"},
    { "DTMETH",  2,  35, 3, 2, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	"SYNCHROTRON"}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"RADIATION SOURCE               :"},
    { "DTMETH",  5,  35, 3, 16, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	"SYNCHROTRON NULL"}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BEAMLINE                       :"},
    { "DTMEAS",  4,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"X-RAY GENERATOR MODEL          :"},
    { "RADIAT",  6,  35, 3, 20, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"MONOCHROMATIC OR LAUE    (M/L) :"},
    { "RADIAT",  3,  35, 3, 4, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"WAVELENGTH OR RANGE        (A) :"},
    { "WAVLEN",  3,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"MONOCHROMATOR                  :"},
    { "RADIAT",  4,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"OPTICS                         :"},
    { "RADIAT",  5,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
 },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DETECTOR TYPE                  :"},
    { "DTMETH",  4,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DETECTOR MANUFACTURER          :"},
    { "DTMETH",  3,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"INTENSITY-INTEGRATION SOFTWARE :"},
    { "DTMEAS",  5,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DATA SCALING SOFTWARE          :"},
    { "DTMEAS",  6,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
 },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"NUMBER OF UNIQUE REFLECTIONS   :"},
    { "REFLEC",  5,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"RESOLUTION RANGE HIGH      (A) :"},
    { "RESTOT",  3,  35, 2, 7, 3,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"RESOLUTION RANGE LOW       (A) :"},
    { "RESTOT",  2,  35, 2, 7, 3,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"REJECTION CRITERIA  (SIGMA(I)) :"},
    { "REFLEC",  3,  35, 2, 7, 3,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
 },
 { 200, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"OVERALL."}}
 },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"COMPLETENESS FOR RANGE     (%) :"},
    { "PERCOM",  2,  35, 2, 7, 1,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DATA REDUNDANCY                :"},
    { "RESTOT",  7,  35, 2, 5, 3,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_R_VALUE, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"R MERGE                    (I) :"},
    { "RESTOT",  6,  35, 2, 10, 5,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_R_VALUE, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"R SYM                      (I) :"},
    { "RESTOT",  8,  35, 2, 7, 5,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"<I/SIGMA(I)> FOR THE DATA SET  :"},
    { "RESTOT",  9,  35, 2, 8, 4,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
 },
 { 200, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"IN THE HIGHEST RESOLUTION SHELL."}}
 },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"HIGHEST RESOLUTION SHELL, RANGE HIGH (A) :"},
    { "SHELC",  2,  45, 2, 5, 2,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"HIGHEST RESOLUTION SHELL, RANGE LOW  (A) :"},
    { "SHELC",  3,  45, 2, 5, 2,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"COMPLETENESS FOR SHELL     (%) :"},
    { "SHELC",  4,  35, 2, 5, 1,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DATA REDUNDANCY IN SHELL       :"},
    { "SHELC",  5,  35, 2, 5, 2,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_R_VALUE, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"R MERGE FOR SHELL          (I) :"},
    { "SHELC",  6,  35, 2, 10, 5,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_R_VALUE, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"R SYM FOR SHELL            (I) :"},
    { "SHELC",  8,  35, 2, 7, 5,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"<I/SIGMA(I)> FOR SHELL         :"},
    { "SHELC",  7,  35, 2, 5, 3,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
 },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "DIFFRACTION PROTOCOL:"},
    { "DIFPTL",  3, 23, 3, 39, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"METHOD USED TO DETERMINE THE STRUCTURE:"},
    { "XFILE1",  3,  41, 3, 21, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SOFTWARE USED:"},
    { "DTMEA1",  2,  16, 3, 44, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"STARTING MODEL:"},
    { "XFILE2",  2,  17, 3, 43, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
 },
 { 200, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
 },
 { 200, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"REMARK:"},
    { "CRMET1",  3,  9, 3, 51, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
};

REMARKS  Remark_205[Num_Remark_205] =
{
 { 205, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"FIBER DIFFRACTION"}}
 },
 { 205, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"THE COORDINATES IN THIS ENTRY WERE GENERATED FROM FIBER"}}
 },
 { 205, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DIFFRACTION DATA.  PROTEIN DATA BANK CONVENTIONS REQUIRE"}}
 },
 { 205, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"THAT CRYST1 AND SCALE RECORDS BE INCLUDED, BUT THE"}}
 },
 { 205, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"VALUES ON THESE RECORDS ARE MEANINGLESS."}}
 }
};

REMARKS  Remark_210[Num_Remark_210] =
{
 { 210, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"EXPERIMENTAL DETAILS"}}
 },
 { 210, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"EXPERIMENT TYPE                : NMR"}}
  },
 { 210, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"TEMPERATURE           (KELVIN) :"},
    { "NMRSMP",  4,  35, 3, 25, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
 { 210, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"PH                             :"},
    { "NMRSMP",  6,  35, 3, 25, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
 { 210, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"IONIC STRENGTH                 :"},
    { "NMRSMP",  7,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 210, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"PRESSURE                       :"},
    { "NMRSMP",  5,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 210, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SAMPLE CONTENTS                :"},
    { "NMRSDT",  3,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 210, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
 },
 { 210, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"NMR EXPERIMENTS CONDUCTED      :"},
    { "NMRECD",  3,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 210, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SPECTROMETER FIELD STRENGTH    :"},
    { "NMRSPM",  5,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 210, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SPECTROMETER MODEL             :"},
    { "NMRSPM",  3,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 210, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SPECTROMETER MANUFACTURER      :"},
    { "NMRSPM",  4,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 210, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
 },
 { 210, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"STRUCTURE DETERMINATION."},}
 },
 { 210, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SOFTWARE USED                 :"},
    { "NMRSFT",  3,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 210, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"METHOD USED                   :"},
    { "XFILE1",  3,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 210, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
 },
 { 210, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"CONFORMERS, NUMBER CALCULATED   :"},
    { "NMRSMB",  2,  35, 1, 4, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 210, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"CONFORMERS, NUMBER SUBMITTED    :"},
    { "NMRSMB",  3,  35, 1, 4, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 210, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"CONFORMERS, SELECTION CRITERIA  :"},
    { "NMRSMB",  4,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 210, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
 },
 { 210, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BEST REPRESENTATIVE CONFORMER IN THIS ENSEMBLE :"},
    { "NMRPST",  2,  50, 3, 10, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 210, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
 },
 { 210, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"REMARK:"},
    { "NMREMK",  3,  9, 3, 51, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
};  

REMARKS  Remark_215[Num_Remark_215] =
{
  { 215, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NMR STUDY"}}
  },
  { 215, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "THE COORDINATES IN THIS ENTRY WERE GENERATED FROM SOLUTION"}}
  },
  { 215, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "NMR DATA.  PROTEIN DATA BANK CONVENTIONS REQUIRE THAT"}}
  },
  { 215, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "CRYST1 AND SCALE RECORDS BE INCLUDED, BUT THE VALUES ON"}}
  },
  { 215, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "THESE RECORDS ARE MEANINGLESS."}}
  }
};

REMARKS  Remark_217[Num_Remark_217] =
{
  { 217, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SOLID STATE NMR STUDY"}}
  },
  { 217, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "THE COORDINATES IN THIS ENTRY WERE GENERATED FROM SOLID"}}
  },
  { 217, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "STATE NMR DATA. PROTEIN DATA BANK CONVENTIONS REQUIRE THAT"}}
  },
  { 217, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "CRYST1 AND SCALE RECORDS BE INCLUDED, BUT THE VALUES ON"}}
  },
  { 217, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "THESE RECORDS ARE MEANINGLESS."}}
  }
};

REMARKS  Remark_220[Num_Remark_220] =
{
  { 220, 1, PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "EXPERIMENTAL DETAILS"}}
  },
  { 220, 1, PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "EXPERIMENT TYPE                : THEORETICAL MODELLING"}}
  },
  { 220, 1, PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 220, 2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "REMARK:"},
    { "MODTLS",  3,  9, 3, 51, 0,  PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  }
};

REMARKS  Remark_225[Num_Remark_225] =
{
  { 225, 1, PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "THEORETICAL MODEL"}}
  },
  { 225, 1, PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "THE COORDINATES IN THIS ENTRY REPRESENT A MODEL STRUCTURE."}}
  },
  { 225, 1, PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "PROTEIN DATA BANK CONVENTIONS REQUIRE THAT CRYST1 AND"}}
  },
  { 225, 1, PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SCALE RECORDS BE INCLUDED, BUT THE VALUES ON THESE"}}
  },
  { 225, 1, PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "RECORDS ARE MEANINGLESS."}}
  }
};

REMARKS  Remark_230[Num_Remark_230] =
{
 { 230, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"EXPERIMENTAL DETAILS"}}
 },
 { 230, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"EXPERIMENT TYPE                : NEUTRON DIFFRACTION"}}
 },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DATE OF DATA COLLECTION        :"},
    { "DTMEAS",  2,  35, 3, 11,0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
 },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"TEMPERATURE           (KELVIN) :"},
    { "DTTEMP",  2,  35, 2, 7, 1,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
 },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"PH                             :"},
    { "PHVAL",  2,  35, 2, 5, 2,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"NUMBER OF CRYSTALS USED        :"},
    { "CRMETH",  4,  35, 1, 4, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
 },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"NEUTRON SOURCE                 :"},
    { "DTMETH",  5,  35, 3, 16, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"BEAMLINE                       :"},
    { "DTMEAS",  4,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"WAVELENGTH OR RANGE        (A) :"},
    { "WAVLEN",  3,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"MONOCHROMATOR                  :"},
    { "RADIAT",  4,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"OPTICS                         :"},
    { "RADIAT",  5,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
 },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DETECTOR TYPE                  :"},
    { "DTMETH",  4,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DETECTOR MANUFACTURER          :"},
    { "DTMETH",  3,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"INTENSITY-INTEGRATION SOFTWARE :"},
    { "DTMEAS",  5,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DATA SCALING SOFTWARE          :"},
    { "DTMEAS",  6,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
 },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"NUMBER OF UNIQUE REFLECTIONS   :"},
    { "REFLEC",  5,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"RESOLUTION RANGE HIGH      (A) :"},
    { "RESTOT",  3,  35, 2, 7, 3,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"RESOLUTION RANGE LOW       (A) :"},
    { "RESTOT",  2,  35, 2, 7, 3,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"REJECTION CRITERIA  (SIGMA(I)) :"},
    { "REFLEC",  3,  35, 2, 7, 3,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
 },
 { 230, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"OVERALL."}}
 },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"COMPLETENESS FOR RANGE     (%) :"},
    { "PERCOM",  2,  35, 2, 7, 1,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DATA REDUNDANCY                :"},
    { "RESTOT",  7,  35, 2, 5, 3,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 2,PDB_REMARK_R_VALUE, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"R MERGE                    (I) :"},
    { "RESTOT",  6,  35, 2, 10, 5,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 2,PDB_REMARK_R_VALUE, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"R SYM                      (I) :"},
    { "RESTOT",  8,  35, 2, 7, 5,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"<I/SIGMA(I)> FOR THE DATA SET  :"},
    { "RESTOT",  9,  35, 2, 8, 4,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
 },
 { 230, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"IN THE HIGHEST RESOLUTION SHELL."}}
 },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"HIGHEST RESOLUTION SHELL, RANGE HIGH (A) :"},
    { "SHELC",  2,  45, 2, 5, 2,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"HIGHEST RESOLUTION SHELL, RANGE LOW  (A) :"},
    { "SHELC",  3,  45, 2, 5, 2,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"COMPLETENESS FOR SHELL     (%) :"},
    { "SHELC",  4,  35, 2, 5, 1,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DATA REDUNDANCY IN SHELL       :"},
    { "SHELC",  5,  35, 2, 5, 2,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 2,PDB_REMARK_R_VALUE, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"R MERGE FOR SHELL          (I) :"},
    { "SHELC",  6,  35, 2, 10, 5,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 2,PDB_REMARK_R_VALUE, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"R SYM FOR SHELL            (I) :"},
    { "SHELC",  8,  35, 2, 7, 5,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"<I/SIGMA(I)> FOR SHELL         :"},
    { "SHELC",  7,  35, 2, 5, 3,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
 },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"METHOD USED TO DETERMINE THE STRUCTURE:"},
    { "XFILE1",  3,  41, 3, 21, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SOFTWARE USED :"},
    { "DTMEA1",  2,  17, 3, 44, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"STARTING MODEL:"},
    { "XFILE2",  2,  17, 3, 43, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
 },
 { 230, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
 },
 { 230, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"REMARK:"},
    { "CRMET1",  3,  9, 3, 51, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  }
};

REMARKS  Remark_240[Num_Remark_240] =
{
 { 240, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"EXPERIMENTAL DETAILS"},}
 },
 { 240, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "RECONSTRUCTION METHOD          :"},
    { "EMEXPT",  2,  36, 3, 23, 0,  PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
 },
 { 240, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SAMPLE TYPE                    :"},
    { "EMSAMP",  2,  36, 3, 23, 0,  PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
 },
 { 240, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "SPECIMEN TYPE                  :"},
    { "EMEXPT",  3,  36, 3, 23, 0,  PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
 },
 { 240, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "DATA ACQUISITION"},}
 },
 { 240, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DATE OF DATA COLLECTION        :"},
    { "DTMEAS",  2,  36, 3, 11,0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 240, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"TEMPERATURE           (KELVIN) :"},
    { "DTTEMP",  2,  36, 2, 7, 1,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 240, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"PH                             :"},
    { "PHVAL",  2,  36, 2, 5, 2,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 240, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"NUMBER OF CRYSTALS USED        :"},
    { "CRMETH",  4,  36, 1, 4, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 240, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "MICROSCOPE MODEL               :"},
    { "EMDATA",  5,  36, 3, 20, 0,  PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
 { 240, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DETECTOR TYPE                  :"},
    { "EMDATA",  6,  36, 3, 23, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 240, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "ACCELERATION VOLTAGE (KV)      :"},
    { "EMDATA",  18, 36, 1, 6, 0,  PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
 { 240, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"NUMBER OF UNIQUE REFLECTIONS   :"},
    { "REFLEC",  5,  36, 3, 23, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 240, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"RESOLUTION RANGE HIGH      (A) :"},
    { "RESTOT",  3,  36, 2, 7, 3,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 240, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"RESOLUTION RANGE LOW       (A) :"},
    { "RESTOT",  2,  36, 2, 7, 3,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 240, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DATA SCALING SOFTWARE          :"},
    { "DTMEAS",  6,  36, 3, 23, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 240, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"COMPLETENESS FOR RANGE     (%) :"},
    { "PERCOM",  2,  36, 2, 7, 1,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 240, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DATA REDUNDANCY                :"},
    { "RESTOT",  7,  36, 2, 5, 3,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 240, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"IN THE HIGHEST RESOLUTION SHELL"}}
 },
 { 240, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"HIGHEST RESOLUTION SHELL, RANGE HIGH (A) :"},
    { "SHELC",  2,  45, 2, 5, 2,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 240, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"HIGHEST RESOLUTION SHELL, RANGE LOW  (A) :"},
    { "SHELC",  3,  45, 2, 5, 2,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 240, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"COMPLETENESS FOR SHELL     (%) :"},
    { "SHELC",  4,  36, 2, 5, 1,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 240, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DATA REDUNDANCY IN SHELL       :"},
    { "SHELC",  5,  36, 2, 5, 2,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 240, 2,PDB_REMARK_R_VALUE, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"R MERGE FOR SHELL          (I) :"},
    { "SHELC",  6,  36, 2, 10, 5,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 240, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"METHOD USED TO DETERMINE THE STRUCTURE:"},
    { "XFILE1",  3,  43, 3, 16, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 240, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SOFTWARE USED                  :"},
    { "DTMEA1",  2,  36, 3, 23, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
 { 240, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  3,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"STARTING MODEL                 :"},
    { "XFILE2",  2,  36, 3, 23, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
 }
};

REMARKS  Remark_250[Num_Remark_250] =
{
  { 250, 1, PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "EXPERIMENTAL DETAILS"}}
  },
  { 250, 2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"EXPERIMENT TYPE                :"},
    { "EXPDTA",  3,  35, 3, 25, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
  { 250, 2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"DATE OF DATA COLLECTION        :"},
    { "DTMEAS",  2,  35, 3, 11,0,  PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
  { 250, 1, PDB_REMARK_GENERAL, 0, {
    { "",             0,  2,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        ""}}
  },
  { 250, 2, PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "REMARK:"},
    { "CRMET1",  3,  9, 3, 51, 0,  PDB_REMARK_LEFT_JUSTIFIED, 
        ""}}
  }
};

REMARKS  Remark_280[Num_Remark_280] =
{
  { 280, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"CRYSTAL"},}
  },
  { 280, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SOLVENT CONTENT, VS   (%):"},
    { "PERWAT",  2,  28, 2, 5, 2, PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
  { 280, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"MATTHEWS COEFFICIENT, VM (ANGSTROMS**3/DA):"},
    { "DENSTY",  3,  45, 2, 4, 2, PDB_REMARK_LEFT_JUSTIFIED, 
	""}}
  },
  { 280, 1, PDB_REMARK_GENERAL, 0, {
    { "", 0, 1, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED, ""}}
  },
  { 280, 2, PDB_REMARK_GENERAL, 0, {
    { "", 0, 1, 3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
        "CRYSTALLIZATION CONDITIONS:"},
    { "CRDTLS", 3, 29, 3, 31, 0, PDB_REMARK_LEFT_JUSTIFIED, ""}}
  }
};

REMARKS Remark_290[Num_Remark_290] = {
  { 290, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"CRYSTALLOGRAPHIC SYMMETRY"}}
  },
  { 290, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SYMMETRY OPERATORS FOR SPACE GROUP: "},
    { "CRYST1",            11,  37,  3, 15, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 290, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 290, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  6,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SYMOP   SYMMETRY"}}
  },
  { 290, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  5,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"NNNMMM   OPERATOR"}}
  },
  { 290, 2,PDB_REMARK_NNNMMM_SYMOP, 0, {
    { "XXX",             0,  5,  3, 6, 0, PDB_REMARK_RIGHT_JUSTIFIED,
	"NNNMMM"},
    { "XXX",             0,  14,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SYMOP"}  }
  },
  { 290, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 290, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  5,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"WHERE NNN -> OPERATOR NUMBER"}}
  },
  { 290, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  11,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"MMM -> TRANSLATION VECTOR"}}
  },
  { 290, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 290, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"CRYSTALLOGRAPHIC SYMMETRY TRANSFORMATIONS"}}
  },
  { 290, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"THE FOLLOWING TRANSFORMATIONS OPERATE ON THE ATOM/HETATM"}}
  },
  { 290, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"RECORDS IN THIS ENTRY TO PRODUCE CRYSTALLOGRAPHICALLY"}}
  },
  { 290, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"RELATED MOLECULES."}}
  },
  { 290, 6,PDB_REMARK_XYZT, 0, {
    { "XXX",             0,  3,  3, 6, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"SMTRY"},
    { "XXX",             0,  11,  1, 2, 0, PDB_REMARK_RIGHT_JUSTIFIED,
	"Serial_No"},
    { "XXX",             0,  13,  2, 10, 6, PDB_REMARK_RIGHT_JUSTIFIED,
	"X"},
    { "XXX",             0,  23,  2, 10, 6, PDB_REMARK_RIGHT_JUSTIFIED,
	"Y"},
    { "XXX",             0,  33,  2, 10, 6, PDB_REMARK_RIGHT_JUSTIFIED,
	"Z"},
    { "XXX",             0,  49,  2,  9, 5, PDB_REMARK_RIGHT_JUSTIFIED,
	"T"} }
  },
  { 290, 1,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  },
  { 290, 2,PDB_REMARK_GENERAL, 0, {
    { "",             0,  1,  3, 0, 0, PDB_REMARK_LEFT_JUSTIFIED,
	"REMARK:" },
    { "XXX",  0, 9, 3, 0, 0,  PDB_REMARK_LEFT_JUSTIFIED,
	""}}
  }, 
};
