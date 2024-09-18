/*
FILE:     BaseParameter_global.h
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
#ifndef _H_BASE_PARAMETER_GLOBAL_H_
#define _H_BASE_PARAMETER_GLOBAL_H_

typedef struct {
        int num_atom;
        char atomname[12][4];
        double ref[12][3];
} REF_NA;

static REF_NA STAND_R = { 9,
        { "N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4", "", "", "" },
        { { -2.250,  5.000,  0.250 },
          { -2.250,  0.500,  0.250 },
          { -2.250,  0.500, -0.250 },
          { -2.250,  5.000, -0.250 },
          {  2.250,  5.000,  0.250 },
          {  2.250,  0.500,  0.250 },
          {  2.250,  0.500, -0.250 },
          {  2.250,  5.000, -0.250 },
          {  0.000,  0.000,  0.000 },
          {  0.000,  0.000,  0.000 },
          {  0.000,  0.000,  0.000 },
          {  0.000,  0.000,  0.000 }
        }
};

static REF_NA STAND_Y = { 6,
        { "N1", "C2", "N3", "C4", "C5", "C6", "", "", "", "", "", "" }, 
        { { -2.250,  5.000,  0.250 },
          { -2.250,  2.000,  0.250 },
          { -2.250,  2.000, -0.250 },
          { -2.250,  5.000, -0.250 },
          {  2.250,  5.000,  0.250 },
          {  2.250,  2.000,  0.250 },
          {  2.250,  2.000, -0.250 },
          {  2.250,  5.000, -0.250 },
          {  0.000,  0.000,  0.000 },
          {  0.000,  0.000,  0.000 },
          {  0.000,  0.000,  0.000 },
          {  0.000,  0.000,  0.000 }
        }
};

static REF_NA STAND_A = { 11, 
        { "C1'", "C2", "C4", "C5", "C6", "C8", "N1", "N3", "N6", "N7", "N9", "" },
        { { -2.479,  5.346,  0.000 },
          { -1.912,  1.023,  0.000 },
          { -1.267,  3.124,  0.000 },
          {  0.071,  2.771,  0.000 },
          {  0.369,  1.398,  0.000 },
          {  0.024,  4.897,  0.000 },
          { -0.668,  0.532,  0.000 },
          { -2.320,  2.290,  0.000 },
          {  1.611,  0.909,  0.000 },
          {  0.877,  3.902,  0.000 },
          { -1.291,  4.498,  0.000 },
          {  0.000,  0.000,  0.000 }
        }
};

static REF_NA STAND_C = { 9,
        { "C1'", "C2", "C4", "C5", "C6", "N1", "N3", "N4", "O2", "", "", "" },
        { {  -2.477,  5.402,  0.000 },
          {  -1.472,  3.158,  0.000 },
          {   0.837,  2.868,  0.000 },
          {   1.056,  4.275,  0.000 },
          {  -0.023,  5.068,  0.000 },
          {  -1.285,  4.542,  0.000 },
          {  -0.391,  2.344,  0.000 },
          {   1.875,  2.027,  0.001 },
          {  -2.628,  2.709,  0.001 },
          {   0.000,  0.000,  0.000 },
          {   0.000,  0.000,  0.000 },
          {   0.000,  0.000,  0.000 }
        }
};

static REF_NA STAND_G = { 12,
        { "C1'", "C2", "C4", "C5", "C6", "C8", "N1", "N2", "N3", "N7", "N9", "O6" },
        { {  -2.477,  5.399,  0.000 },
          {  -1.999,  1.087,  0.000 },
          {  -1.265,  3.177,  0.000 },
          {   0.071,  2.833,  0.000 },
          {   0.424,  1.460,  0.000 },
          {   0.023,  4.962,  0.000 },
          {  -0.700,  0.641,  0.000 },
          {  -2.949,  0.139, -0.001 },
          {  -2.342,  2.364,  0.001 },
          {   0.870,  3.969,  0.000 },
          {  -1.289,  4.551,  0.000 },
          {   1.554,  0.955,  0.000 }
        }
};

static REF_NA STAND_P = { 9,
        { "C1'", "C2", "C4", "C5", "C6", "N1", "N3", "O2", "O4", "", "", "" },
        {
          {  -2.506,  5.371,  0.000 },
          {   1.037,  2.915,  0.000 },
          {  -1.422,  3.076,  0.000 },
          {  -1.284,  4.500,  0.000 },
          {  -0.064,  5.048,  0.000 },
          {   1.087,  4.295,  0.000 },
          {  -0.229,  2.383,  0.000 },
          {   2.036,  2.217,  0.000 },
          {  -2.485,  2.453,  0.000 },
          {   0.000,  0.000,  0.000 },
          {   0.000,  0.000,  0.000 },
          {   0.000,  0.000,  0.000 }
        }
};

static REF_NA STAND_T = { 10,
        { "C1'", "C2", "C4", "C5", "C7", "C6", "N1",  "N3", "O2", "O4", "", "" },
        {
          {  -2.481,  5.354,  0.000 },
          {  -1.462,  3.135,  0.000 },
          {   0.994,  2.897,  0.000 },
          {   1.106,  4.338,  0.000 },
          {   2.466,  4.961,  0.001 },
          {  -0.024,  5.057,  0.000 },
          {  -1.284,  4.500,  0.000 },
          {  -0.298,  2.407,  0.000 },
          {  -2.562,  2.608,  0.000 },
          {   1.944,  2.119,  0.000 },
          {   0.000,  0.000,  0.000 },
          {   0.000,  0.000,  0.000 }
        }
};

static REF_NA STAND_U = { 9,
        { "C1'",  "C2", "C4", "C5", "C6", "N1", "N3" ,  "O2", "O4", "", "", "" },
        {
          {  -2.481,  5.354,  0.000 },
          {  -1.462,  3.131,  0.000 },
          {   0.989,  2.884,  0.000 },
          {   1.089,  4.311,  0.000 },
          {  -0.024,  5.053,  0.000 },
          {  -1.284,  4.500,  0.000 },
          {  -0.302,  2.397,  0.000 },
          {  -2.563,  2.608,  0.000 },
          {   1.935,  2.094, -0.001 },
          {   0.000,  0.000,  0.000 },
          {   0.000,  0.000,  0.000 },
          {   0.000,  0.000,  0.000 }
        }
};

typedef struct {
       char na_type[4];
       REF_NA *ref_type;
} REF_NA_TYPE;

#define NUM_REF_NA  9

static REF_NA_TYPE ref_na_type[NUM_REF_NA] = {
       { "A", &STAND_A },
       { "C", &STAND_C },
       { "G", &STAND_G },
       { "I", &STAND_G },
       { "P", &STAND_P },
       { "R", &STAND_R },
       { "T", &STAND_T },
       { "U", &STAND_U },
       { "Y", &STAND_Y }
};

#define NUM_RES_MAPPING   32

static const char *_res_code_mapping[NUM_RES_MAPPING][2] = {
       { "1AP", "G" },
       { "6MP", "A" },
       { "ADE", "A" },
       { "AMD", "A" },
       { "AMP", "A" },
       { "APN", "A" },
       { "ATP", "A" },
       { "AZT", "T" },
       { "CMP", "A" },
       { "CPN", "C" },
       { "CYT", "C" },
       { "DAD", "A" },
       { "DG3", "G" },
       { "GNP", "G" },
       { "GPN", "G" },
       { "GUA", "G" },
       { "GUN", "G" },
       { "HPA", "A" },
       { "INO", "I" },
       { "IPN", "T" },
       { "MGT", "G" },
       { "N",   "U" },
       { "P2U", "P" },
       { "PSU", "P" },
       { "QSI", "A" },
       { "SAM", "A" },
       { "TPN", "T" },
       { "TSP", "T" },
       { "TTP", "T" },
       { "UCP", "U" },
       { "URA", "U" },
       { "VAA", "A" }
};

#define NUM_CODE_MAPPING  7

static const char *_code_code_mapping[NUM_CODE_MAPPING][2] = {
       { "A", "R" },
       { "C", "Y" },
       { "G", "R" },
       { "I", "R" },
       { "P", "Y" },
       { "T", "Y" },
       { "U", "Y" }
};

#endif
