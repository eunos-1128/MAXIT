/*
FILE:     SeqCodeUtil_global.h
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
#ifndef _H_SEQCODEUTIL_GLOBAL_H_
#define _H_SEQCODEUTIL_GLOBAL_H_


#include <map>
#include <set>
#include <string>

#define AANUMBER 31

typedef struct {        
        const char* a1let;   /* one letter code   */
        const char* a3let;   /* three letter code */
        const char* fname;   /* full name         */
        const char* type;    /* residue type      */
} __CODES;

static __CODES aacodes[AANUMBER]  = {
       { "X", "UNK", "Unknown",        "POLYPEPTIDE(L)"          }, /*  0 */
       { "A", "ALA", "Alanine",        "POLYPEPTIDE(L)"          }, /*  1 */
       { "R", "ARG", "Arginine",       "POLYPEPTIDE(L)"          }, /*  2 */
       { "N", "ASN", "Asparagine",     "POLYPEPTIDE(L)"          }, /*  3 */
       { "D", "ASP", "Aspartic Acid",  "POLYPEPTIDE(L)"          }, /*  4 */
       { "C", "CYS", "Cysteine",       "POLYPEPTIDE(L)"          }, /*  5 */
       { "Q", "GLN", "Glutamine",      "POLYPEPTIDE(L)"          }, /*  6 */
       { "E", "GLU", "Glutamic Acid",  "POLYPEPTIDE(L)"          }, /*  7 */
       { "G", "GLY", "Glycine",        "POLYPEPTIDE(L)"          }, /*  8 */
       { "H", "HIS", "Histidine",      "POLYPEPTIDE(L)"          }, /*  9 */
       { "I", "ILE", "Isoleucine",     "POLYPEPTIDE(L)"          }, /* 10 */
       { "L", "LEU", "Leucine",        "POLYPEPTIDE(L)"          }, /* 11 */
       { "K", "LYS", "Lysine",         "POLYPEPTIDE(L)"          }, /* 12 */
       { "M", "MET", "Methionine",     "POLYPEPTIDE(L)"          }, /* 13 */
       { "F", "PHE", "Phenylalanine",  "POLYPEPTIDE(L)"          }, /* 14 */
       { "P", "PRO", "Proline",        "POLYPEPTIDE(L)"          }, /* 15 */
       { "S", "SER", "Serine",         "POLYPEPTIDE(L)"          }, /* 16 */
       { "T", "THR", "Threonine",      "POLYPEPTIDE(L)"          }, /* 17 */
       { "W", "TRP", "Tryptophane",    "POLYPEPTIDE(L)"          }, /* 18 */
       { "Y", "TYR", "Tyrosine",       "POLYPEPTIDE(L)"          }, /* 19 */
       { "V", "VAL", "Valine",         "POLYPEPTIDE(L)"          }, /* 20 */
       { "O", "PYL", "Pyrrolysine",    "POLYPEPTIDE(L)"          }, /* 21 */
       { "U", "SEC", "Selenocysteine", "POLYPEPTIDE(L)"          }, /* 22 */
       { "-", "---", "Gap",            "POLYPEPTIDE(L)"          }, /* 23 */
       { "*", "END", "End",            "POLYPEPTIDE(L)"          }, /* 24 */
       { "A", "DA",  "DEOXYADENOSINE", "POLYDEOXYRIBONUCLEOTIDE" }, /* 25 */
       { "C", "DC",  "DEOXYCYTIDINE",  "POLYDEOXYRIBONUCLEOTIDE" }, /* 26 */
       { "G", "DG",  "DEOXYGUANOSINE", "POLYDEOXYRIBONUCLEOTIDE" }, /* 27 */
       { "I", "DI",  "DEOXYINOSINE",   "POLYDEOXYRIBONUCLEOTIDE" }, /* 28 */
       { "T", "DT",  "THYMIDINE",      "POLYDEOXYRIBONUCLEOTIDE" }, /* 29 */
       { "U", "DU",  "DEOXYURIDINE",   "POLYDEOXYRIBONUCLEOTIDE" }  /* 30 */
};

#define NUM_DNA_RESIDUE      6

static const char *_dna_residues[NUM_DNA_RESIDUE] = {
       "DA", "DC", "DG", "DI", "DT", "DU"
};

static const char *_s_dna_residues[NUM_DNA_RESIDUE] = {
       "A", "C", "G", "I", "T", "U"
};

#define NUM_RNA_RESIDUE      5

static const char *_rna_residues[NUM_RNA_RESIDUE] = {
       "A", "C", "G", "I", "U"
};

#define NUM_AA_RESIDUE      22

static const char *_aa_residues[NUM_AA_RESIDUE] = {
       "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
       "LYS", "MET", "PHE", "PRO", "PYL", "SEC", "SER", "THR", "TRP", "TYR", "VAL"
};

#define NUM_DAA_RESIDUE     19

static const char *_daa_residues[NUM_DAA_RESIDUE] = {
       "DAL", "DAR", "DSG", "DAS", "DCY", "DGN", "DGL", "DHI", "DIL", "DLE",
       "DLY", "MED", "DPN", "DPR", "DSN", "DTH", "DTR", "DTY", "DVA"
};

#define NUM_MIXED_RESIDUE    2

static const char *_mixed_residues[NUM_MIXED_RESIDUE] = { "ASX", "GLX" };

#define NUM_UNKNOWN_RESIDUE  2

static const char *_unknown_residues[NUM_UNKNOWN_RESIDUE] = { "N", "UNK" };

#define NUM_WATER_RESIDUE    7 

static const char *_water_residues[NUM_WATER_RESIDUE] = {
       "H2O", "WAT", "TIP", "WTR", "OH2", "HOH", "DOD"
};

#define NUM_NN_RESIDUE       5

static const char *_nn_residues[NUM_NN_RESIDUE] = {
       "URI", "THY", "ADE", "CYT", "GUA"
};

// #define NUM_SPECIAL_NN_RESIDUE  14
#define NUM_SPECIAL_NN_RESIDUE   6

static const char *_special_nn_residues[NUM_SPECIAL_NN_RESIDUE] = {
       /* "AD", "CD", "GD", "TD", "AR", "CR", "GR", "UR", */ "URA", "URI", "THY", "ADE", "CYT", "GUA"
};

// #define NUM_MAPPING_NN_RESIDUE  11
#define NUM_MAPPING_NN_RESIDUE   3

static const char *_mapping_nn_residues[NUM_MAPPING_NN_RESIDUE] = {
       /* "DA", "DC", "DG", "DT", "A",  "C",  "G",  "U", */ "U",  "U",  "DT"
};

#define NUM_SIMILAR_RESIDUE  2

static const char *_similar_residues[NUM_SIMILAR_RESIDUE][2] = {
       { "ASN", "ASP" },
       { "GLN", "GLU" }
};

static std::map<std::string, std::string> _a1_a3_mapping;
static std::map<std::string, std::string> _a3_a1_mapping;
static std::set<std::string> _standard_residue_list_set1;
static std::set<std::string> _standard_residue_list_set2;
static std::set<std::string> _standard_aa_residue_list_set;
static std::set<std::string> _standard_Laa_residue_list_set;
static std::set<std::string> _standard_Daa_residue_list_set;
static std::set<std::string> _standard_na_residue_list_set;
static std::set<std::string> _standard_dna_residue_list_set;
static std::set<std::string> _standard_rna_residue_list_set;
static std::set<std::string> _water_residue_set;
static std::map<std::string, std::string> _d1_d2_mapping;
static std::set<std::string> _similar_residue_set;
static std::set<std::string> _special_na_residue_list_set;
static std::map<std::string, std::string> _special_na_residue_mapping;

#endif
