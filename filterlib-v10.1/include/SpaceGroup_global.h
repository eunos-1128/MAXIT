/*
FILE:     SpaceGroup_global.h
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
#ifndef _H_SPACE_GROUP_GLOBAL_H_
#define _H_SPACE_GROUP_GLOBAL_H_

#include <map>
// #include <set>
#include <string>
#include <vector>

#define spgbinfile            "space_group.odb"

#define NUM_ITEMS        15

static const char *_items[NUM_ITEMS] = {
       "matrix[1][1]",
       "matrix[1][2]",
       "matrix[1][3]",
       "matrix[2][1]",
       "matrix[2][2]",
       "matrix[2][3]",
       "matrix[3][1]",
       "matrix[3][2]",
       "matrix[3][3]",
       "trans[1][1]",
       "trans[1][2]",
       "trans[2][1]",
       "trans[2][2]",
       "trans[3][1]",
       "trans[3][2]"
};
/*
#define NUM_SPACE_GROUP  95

static const char *space_group_name[NUM_SPACE_GROUP][2] = {
       { "A1",         "A 1"         },
       { "A121",       "A 1 2 1"     },
       { "A2",         "A 2"         },
       { "B112",       "B 1 1 2"     },
       { "B2",         "B 2"         },
       { "B2212",      "B 2 21 2"    },
       { "C121",       "C 1 2 1"     },
       { "C1211",      "C 1 21 1"    },
       { "C2",         "C 2"         },
       { "C2(A112)",   "C 2(A 112)"  },
       { "C21",        "C 21"        },
       { "C222",       "C 2 2 2"     },
       { "C2221",      "C 2 2 21"    },
       { "C4212",      "C 4 21 2"    },
       { "F222",       "F 2 2 2"     },
       { "F23",        "F 2 3"       },
       { "F4132",      "F 41 3 2"    },
       { "F422",       "F 4 2 2"     },
       { "F432",       "F 4 3 2"     },
       { "H3",         "H 3"         },
       { "H32",        "H 3 2"       },
       { "I121",       "I 1 2 1"     },
       { "I1211",      "I 1 21 1"    },
       { "I2",         "I 2"         },
       { "I21",        "I 21"        },
       { "I212121",    "I 21 21 21"  },
       { "I213",       "I 21 3"      },
       { "I222",       "I 2 2 2"     },
       { "I23",        "I 2 3"       },
       { "I4",         "I 4"         },
       { "I41",        "I 41"        },
       { "I4122",      "I 41 2 2"    },
       { "I4132",      "I 41 3 2"    },
       { "I422",       "I 4 2 2"     },
       { "I432",       "I 4 3 2"     },
       { "P-1",        "P -1"        },
       { "P1",         "P 1"         },
       { "P1-",        "P -1"        },
       { "P112",       "P 1 1 2"     },
       { "P1121",      "P 1 1 21"    },
       { "P1211",      "P 1 21 1"    },
       { "P121",       "P 1 2 1"     },
       { "P2",         "P 2"         },
       { "P21",        "P 21"        },
       { "P21(C)",     "P 21(C)"     },
       { "P21212",     "P 21 21 2"   },
       { "P212121",    "P 21 21 21"  },
       { "P21212A",    "P 21 21 2 A" },
       { "P21221",     "P 21 2 21"   },
       { "P213",       "P 21 3"      },
       { "P22121",     "P 2 21 21"   },
       { "P222",       "P 2 2 2"     },
       { "P2221",      "P 2 2 21"    },
       { "P23",        "P 2 3"       },
       { "P3",         "P 3"         },
       { "P31",        "P 31"        },
       { "P3112",      "P 31 1 2"    },
       { "P312",       "P 3 1 2"     },
       { "P3121",      "P 31 2 1"    },
       { "P32",        "P 32"        },
       { "P321",       "P 3 2 1"     },
       { "P3212",      "P 32 1 2"    },
       { "P3221",      "P 32 2 1"    },
       { "P4",         "P 4"         },
       { "P41",        "P 41"        },
       { "P41212",     "P 41 21 2"   },
       { "P4122",      "P 41 2 2"    },
       { "P4132",      "P 41 3 2"    },
       { "P42",        "P 42"        },
       { "P4212",      "P 4 21 2"    },
       { "P422",       "P 4 2 2"     },
       { "P42212",     "P 42 21 2"   },
       { "P4222",      "P 42 2 2"    },
       { "P4232",      "P 42 3 2"    },
       { "P43",        "P 43"        },
       { "P432",       "P 4 3 2"     },
       { "P43212",     "P 43 21 2"   },
       { "P4322",      "P 43 2 2"    },
       { "P4332",      "P 43 3 2"    },
       { "P6",         "P 6"         },
       { "P61",        "P 61"        },
       { "P6122",      "P 61 2 2"    },
       { "P62",        "P 62"        },
       { "P622",       "P 6 2 2"     },
       { "P6222",      "P 62 2 2"    },
       { "P63",        "P 63"        },
       { "P6322",      "P 63 2 2"    },
       { "P64",        "P 64"        },
       { "P6422",      "P 64 2 2"    },
       { "P65",        "P 65"        },
       { "P6522",      "P 65 2 2"    },
       { "R3",         "R 3"         },
       { "R3(RHOM)",   "R 3(RHOM)"   },
       { "R32",        "R 3 2"       },
       { "R32(RHOM)",  "R 3 2(RHOM)" }
};
*/
std::map<std::string, std::vector<std::vector<int> > > _Data;
// std::set<std::string> _Allowed_names;
std::map<std::string, std::string> _Number_mapping;
std::map<std::string, std::string> _Name_mapping;
std::vector<std::string> _Name_list;

#endif
