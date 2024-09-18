/*
FILE:     ReservedWord_global.h
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
#ifndef _H_RESERVED_WORD_GLOBAL_H_
#define _H_RESERVED_WORD_GLOBAL_H_

#include <map>
#include <string>

static std::map<std::string, std::string> _reserved_word_mapping;

#define NUM_RESERVED_WORD    141

static const char *_reserved_word[NUM_RESERVED_WORD][2] = {
       { "A",         "a"         },
       { "Al",        "al"        },
       { "Amp",       "AMP"       },
       { "An",        "an"        },
       { "And",       "and"       },
       { "Arac",      "AraC"      },
       { "Are",       "are"       },
       { "As",        "as"        },
       { "Asv",       "ASV"       },
       { "At",        "at"        },
       { "Atf",       "ATF"       },
       { "Atpase",    "ATPase"    },
       { "Atp",       "ATP"       },
       { "at-Rich",   "AT-Rich"   },
       { "At-Rich",   "AT-Rich"   },
       { "Bamhi",     "BamHI"     },
       { "Be",        "be"        },
       { "Between",   "between"   },
       { "Bhlh",      "bHLH"      },
       { "Bp",        "BP"        },
       { "By",        "by"        },
       { "Bzip",      "bZIP"      },
       { "Cap",       "CAP"       },
       { "Creb",      "CREB"      },
       { "Cspa",      "CspA"      },
       { "Dcoh",      "DCoH"      },
       { "Ddctp",     "ddCTP"     },
       { "dna",       "DNA"       },
       { "Dna",       "DNA"       },
       { "DNa",       "DNA"       },
       { "Dnase",     "DNase"     },
       { "dsdna",     "DsDNA"     },
       { "Dsdna",     "DsDNA"     },
       { "Dtmp",      "dTMP"      },
       { "Due",       "due"       },
       { "Ebna",      "EBNA"      },
       { "Ecori",     "EcoRI"     },
       { "Ecorv",     "EcoRV"     },
       { "Ef",        "EF"        },
       { "Et",        "et"        },
       { "Febs",      "FEBS"      },
       { "Fis",       "FIS"       },
       { "For",       "for"       },
       { "From",      "from"      },
       { "Gal4",      "GAL4"      },
       { "Gcn4",      "GCN4"      },
       { "Gli",       "GLI"       },
       { "Gtp",       "GTP"       },
       { "Gvp-ssdna", "GVP-ssDNA" },
       { "Gvp-Ssdna", "GVP-ssDNA" },
       { "Hck",       "HCK"       },
       { "Hept",      "HEPT"      },
       { "Hha1",      "HhaI"      },
       { "Hhai",      "HhaI"      },
       { "Hhh",       "HHH"       },
       { "Hiv-1",     "HIV-1"     },
       { "hiv",       "HIV"       },
       { "Hiv",       "HIV"       },
       { "Hlh",       "HLH"       },
       { "Hnf1",      "HNF1"      },
       { "Hnf",       "HNF"       },
       { "Ihf",       "IHF"       },
       { "Ii",        "II"        },
       { "Iii",       "III"       },
       { "In",        "in"        },
       { "Is",        "is"        },
       { "Its",       "its"       },
       { "Kappab",    "kappaB"    },
       { "Kda",       "kDa"       },
       { "Laci",      "LacI"      },
       { "Lfb1",      "LFB1"      },
       { "Lviii",     "LVIII"     },
       { "Mad",       "MAD"       },
       { "Mata1",     "MATa1"     },
       { "Matalpha2", "MATalpha2" },
       { "Mat",       "MAT"       },
       { "Mckenna",   "McKenna"   },
       { "Mckercher", "McKercher" },
       { "Mcpherson", "McPherson" },
       { "Mer",       "mer"       },
       { "Met",       "met"       },
       { "Mmlv",      "MMLV"      },
       { "Ms2",       "MS2"       },
       { "Myod",      "MyoD"      },
       { "Nf",        "NF"        },
       { "Nh2",       "NH2"       },
       { "nmr",       "NMR"       },
       { "Nmr",       "NMR"       },
       { "Not",       "not"       },
       { "Ny",        "NY"        },
       { "Of",        "of"        },
       { "Ompr",      "OmpR"      },
       { "On",        "on"        },
       { "or1",       "OR1"       },
       { "Or1",       "OR1"       },
       { "or2",       "OR2"       },
       { "Or2",       "OR2"       },
       { "Or",        "or"        },
       { "P50",       "p50"       },
       { "Pcna",      "PCNA"      },
       { "Pdtp",      "pdTp"      },
       { "Phix174",   "PhiX174"   },
       { "Pl8",       "PL8"       },
       { "Pou",       "POU"       },
       { "Ppr1",      "PPR1"      },
       { "Purr",      "PurR"      },
       { "Pvuii",     "PvuII"     },
       { "Rap",       "RAP"       },
       { "rna",       "RNA"       },
       { "Rna",       "RNA"       },
       { "Rsref",     "RSRef"     },
       { "Rt",        "RT"        },
       { "ruvc",      "RuvC"      },
       { "Ruvc",      "RuvC"      },
       { "Sh3",       "SH3"       },
       { "Sos",       "SOS"       },
       { "Tafs",      "TAFs"      },
       { "Taqi",      "TaqI"      },
       { "Tata",      "TATA"      },
       { "Tbp",       "TBP"       },
       { "Tfiia",     "TFIIA"     },
       { "Tfiib",     "TFIIB"     },
       { "Tfiid",     "TFIID"     },
       { "That",      "that"      },
       { "The",       "the"       },
       { "Tibo",      "TIBO"      },
       { "To",        "to"        },
       { "Trnaphe",   "tRNAphe"   },
       { "trna",      "tRNA"      },
       { "Trna",      "tRNA"      },
       { "Tu",        "TU"        },
       { "Umud'",     "UmuD'"     },
       { "Usa",       "USA"       },
       { "Usf",       "USF"       },
       { "Which",     "which"     },
       { "Within",    "within"    },
       { "Without",   "without"   },
       { "With",      "with"      },
       { "x-ray",     "X-Ray"     },
       { "X-ray",     "X-Ray"     },
       { "Yy1",       "YY1"       }
};

#endif
