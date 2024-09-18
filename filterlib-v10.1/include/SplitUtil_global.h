/*
FILE:     SplitUtil_global.h
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
#ifndef _H_SPLIT_UTIL_GLOBAL_H_
#define _H_SPLIT_UTIL_GLOBAL_H_

static _SPLIT_RESIDUE_LIST _split_residue_list[NUM_SPLIT_RESIDUE_LIST] = {
       { "1CU", 2,
          {
             { "CU",  1,
                {
                   { "CU",  "CU"  }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             }
          }
       },
       { "2OF", 3,
          {
             { "FE",  1,
                {
                   { "FE",  "FE",  }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             }
          }
       },
       { "3OF", 3,
          {
             { "FE",  1,
                {
                   { "FE",  "FE",  }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             }
          }
       },
       { "543", 8,
          {
             { "CA",  1,
                {
                   { "CA",   "CA"   }
                }
             },
             { "EOH", 9,
                {
                   { "C1",   "C1"   },
                   { "C2",   "C2"   },
                   { "O",    "O"    },
                   { "H11",  "H11"  },
                   { "H12",  "H12"  },
                   { "H21",  "H21"  },
                   { "H22",  "H22"  },
                   { "H23",  "H23"  },
                   { "HO",   "HO"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O4",   "O"   },
                   { "HO41", "H1"  },
                   { "HO42", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O5",   "O"   },
                   { "HO51", "H1"  },
                   { "HO52", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O6",   "O"   },
                   { "HO61", "H1"  },
                   { "HO62", "H2"  }
                }
             }
          }
       },
       { "CD1", 2,
          {
             { "CD",  1,
                {
                   { "CD",   "CD"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             }
          }
       },
       { "CD3", 4,
          {
             { "CD",  1,
                {
                   { "CD",   "CD"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             }
          }
       },
       { "CD5", 6,
          {
             { "CD",  1,
                {
                   { "CD",   "CD"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O4",   "O"   },
                   { "HO41", "H1"  },
                   { "HO42", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O5",   "O"   },
                   { "HO51", "H1"  },
                   { "HO52", "H2"  }
                }
             }
          }
       },
       { "CO5", 6,
          {
             { "CO",  1,
                {
                   { "CO",   "CO"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O4",   "O"   },
                   { "HO41", "H1"  },
                   { "HO42", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O5",   "O"   },
                   { "HO51", "H1"  },
                   { "HO52", "H2"  }
                }
             }
          }
       },
       { "KO4", 5,
          {
             { "K",   1,
                {
                   { "K",    "K"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "H11",  "H1"  },
                   { "H12",  "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "H21",  "H1"  },
                   { "H22",  "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "H31",  "H1"  },
                   { "H32",  "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O4",   "O"   },
                   { "H41",  "H1"  },
                   { "H42",  "H2"  }
                }
             }
          }
       },
       { "MH3", 2,
          {
             { "MN",  1,
                {
                   { "MN",   "MN"   }
                }
             },
             { "HOH", 2,
                {
                   { "O1",   "O"   },
                   { "H2",   "H1"  }
                }
             }
          }
       },
       { "MN5", 6,
          {
             { "MN",  1,
                {
                   { "MN",   "MN"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O4",   "O"   },
                   { "HO41", "H1"  },
                   { "HO42", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O5",   "O"   },
                   { "HO51", "H1"  },
                   { "HO52", "H2"  }
                }
             }
          }
       },
       { "MN6", 7,
          {
             { "MN",  1,
                {
                   { "MN",   "MN"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O4",   "O"   },
                   { "HO41", "H1"  },
                   { "HO42", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O5",   "O"   },
                   { "HO51", "H1"  },
                   { "HO52", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O6",   "O"   },
                   { "HO61", "H1"  },
                   { "HO62", "H2"  }
                }
             }
          }
       },
       { "MO1", 2,
          {
             { "MG",  1,
                {
                   { "MG",   "MG"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             }
          }
       },
       { "MO2", 3,
          {
             { "MG",  1,
                {
                   { "MG",   "MG"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             }
          }
       },
       { "MO3", 4,
          {
             { "MG",  1,
                {
                   { "MG",   "MG"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             }
          }
       },
       { "MO4", 5,
          {
             { "MG",  1,
                {
                   { "MG",   "MG"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O4",   "O"   },
                   { "HO41", "H1"  },
                   { "HO42", "H2"  }
                }
             }
          }
       },
       { "MO5", 6,
          {
             { "MG",  1,
                {
                   { "MG",   "MG"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O4",   "O"   },
                   { "HO41", "H1"  },
                   { "HO42", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O5",   "O"   },
                   { "HO51", "H1"  },
                   { "HO52", "H2"  }
                }
             }
          }
       },
       { "MO6", 7,
          {
             { "MG",  1,
                {
                   { "MG",   "MG"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O4",   "O"   },
                   { "HO41", "H1"  },
                   { "HO42", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O5",   "O"   },
                   { "HO51", "H1"  },
                   { "HO52", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O6",   "O"   },
                   { "HO61", "H1"  },
                   { "HO62", "H2"  }
                }
             }
          }
       },
       { "MW1", 2,
          {
             { "MN",  1,
                {
                   { "MN",   "MN"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             }
          }
       },
       { "MW2", 3, 
          {
             { "MN",  1,
                {
                   { "MN",   "MN"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             }
          }
       },
       { "MW3", 4,
          {
             { "MN",  1,
                {
                   { "MN",   "MN"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             }
          }
       },
       { "NA2", 3,
          {
             { "NA",  1,
                {
                   { "NA",   "NA"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             }
          }
       },
       { "NA5", 6,
          {
             { "NA",  1,
                {
                   { "NA",   "NA"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O4",   "O"   },
                   { "HO41", "H1"  },
                   { "HO42", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O5",   "O"   },
                   { "HO51", "H1"  },
                   { "HO52", "H2"  }
                }
             }
          }
       },
       { "NA6", 7,
          {
             { "NA",  1,
                {
                   { "NA",   "NA"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O4",   "O"   },
                   { "HO41", "H1"  },
                   { "HO42", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O5",   "O"   },
                   { "HO51", "H1"  },
                   { "HO52", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O6",   "O"   },
                   { "HO61", "H1"  },
                   { "HO62", "H2"  }
                }
             }
          }
       },
       { "NAO", 2,
          {
             { "NA",  1,
                {
                   { "NA",   "NA"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             }
          }
       },
       { "NAW", 4,
          {
             { "NA",  1,
                {
                   { "NA",   "NA"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             }
          }
       },
       { "NI1", 2,
          {
             { "NI",  1,
                {
                   { "NI",   "NI"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             }
          }
       },
       { "NI2", 3,
          {
             { "NI",  1,
                {
                   { "NI",   "NI"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "3HO1", "H1"  },
                   { "4HO1", "H2"  }
                }
             }
          }
       },
       { "NI3", 4,
          {
             { "NI",  1,
                {
                   { "NI",   "NI"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             }
          }
       },
       { "NIK", 6,
          {
             { "NI",   1,
                {
                   { "NI",   "NI"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "H11",  "H1"  },
                   { "H12",  "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "H21",  "H1"  },
                   { "H22",  "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "H31",  "H1"  },
                   { "H32",  "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O4",   "O"   },
                   { "H41",  "H1"  },
                   { "H42",  "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O5",   "O"   },
                   { "H51",  "H1"  },
                   { "H52",  "H2"  }
                }
             }
          }
       },
       { "O4M", 5,
          {
             { "MN",  1,
                {
                   { "MN",   "MN"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O4",   "O"   },
                   { "HO41", "H1"  },
                   { "HO42", "H2"  }
                }
             }
          }
       },
       { "OC1", 2,
          {
             { "CA",  1,
                {
                   { "CA",   "CA"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             }
          }
       },
       { "OC2", 3,
          {
             { "CA",  1,
                {
                   { "CA",   "CA"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             }
          }
       },
       { "OC3", 4,
          {
             { "CA",  1,
                {
                   { "CA",   "CA"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             }
          }
       },
       { "OC4", 5,
          {
             { "CA",  1,
                {
                   { "CA",   "CA"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O4",   "O"   },
                   { "HO41", "H1"  },
                   { "HO42", "H2"  }
                }
             }
          }
       },
       { "OC5", 6,
          {
             { "CA",  1,
                {
                   { "CA",   "CA"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O4",   "O"   },
                   { "HO41", "H1"  },
                   { "HO42", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O5",   "O"   },
                   { "HO51", "H1"  },
                   { "HO52", "H2"  }
                }
             }
          }
       },
       { "OC6", 7,
          {
             { "CA",  1,
                {
                   { "CA",   "CA"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O4",   "O"   },
                   { "HO41", "H1"  },
                   { "HO42", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O5",   "O"   },
                   { "HO51", "H1"  },
                   { "HO52", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O6",   "O"   },
                   { "HO61", "H1"  },
                   { "HO62", "H2"  }
                }
             }
          }
       },
       { "OC7", 8,
          {
             { "CA",  1,
                {
                   { "CA",   "CA"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O4",   "O"   },
                   { "HO41", "H1"  },
                   { "HO42", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O5",   "O"   },
                   { "HO51", "H1"  },
                   { "HO52", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O6",   "O"   },
                   { "HO61", "H1"  },
                   { "HO62", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O7",   "O"   },
                   { "HO71", "H1"  },
                   { "HO72", "H2"  }
                }
             }
          }
       },
       { "OC8", 9,
          {
             { "CA",  1,
                {
                   { "CA",   "CA"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O4",   "O"   },
                   { "HO41", "H1"  },
                   { "HO42", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O5",   "O"   },
                   { "HO51", "H1"  },
                   { "HO52", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O6",   "O"   },
                   { "HO61", "H1"  },
                   { "HO62", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O7",   "O"   },
                   { "HO71", "H1"  },
                   { "HO72", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O8",   "O"   },
                   { "HO81", "H1"  },
                   { "HO82", "H2"  }
                }
             }
          }
       },
       { "OCL", 2,
          {
             { "CO",  1,
                {
                   { "CO",   "CO"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             }
          }
       },
       { "OCM", 4,
          {
             { "CO",  1,
                {
                   { "CO",   "CO"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             }
          }
       },
       { "OCN", 3,
          {
             { "CO",  1,
                {
                   { "CO",   "CO"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             }
          }
       },
       { "OCO", 7,
          {
             { "CO",  1,
                {
                   { "CO",   "CO"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O4",   "O"   },
                   { "HO41", "H1"  },
                   { "HO42", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O5",   "O"   },
                   { "HO51", "H1"  },
                   { "HO52", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O6",   "O"   },
                   { "HO61", "H1"  },
                   { "HO62", "H2"  }
                }
             }
          }
       },
       { "OF1", 2,
          {
             { "FE",  1,
                {
                   { "FE",  "FE",  }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             }
          }
       },
       { "OF3", 2,
          {
             { "FE",  1,
                {
                   { "FE",  "FE",  }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             }
          }
       },
       { "YH", 2,
          {
             { "Y",   1,
                {
                   { "Y",    "Y"   }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "H1",   "H1"  },
                   { "H2",   "H2"  }
                }
             }
          }
       },
       { "ZH3", 4,
          {
             { "ZN",  1,
                {
                   { "ZN",  "ZN"   }
                }
             },
             { "HOH", 2,
                {
                   { "OH1",  "O"   },
                   { "HH1",  "H1"  }
                }
             },
             { "HOH", 2,
                {
                   { "OH2",  "O"   },
                   { "HH2",  "H1"  }
                }
             },
             { "HOH", 2,
                {
                   { "OH3",  "O"   },
                   { "HH3",  "H1"  }
                }
             }
          }
       },
       { "ZN3", 2,
          {
             { "ZN",  1,
                {
                   { "ZN",  "ZN",  }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             }
          }
       },
       { "ZNO", 3,
          {
             { "ZN",  1,
                {
                   { "ZN",  "ZN",  }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "3HO1", "H1"  },
                   { "4HO1", "H2"  }
                }
             }
          }
       },
       { "ZO3", 4,
          {
             { "ZN",  1,
                {
                   { "ZN",  "ZN",  }
                }
             },
             { "HOH", 3,
                {
                   { "O1",   "O"   },
                   { "HO11", "H1"  },
                   { "HO12", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O2",   "O"   },
                   { "HO21", "H1"  },
                   { "HO22", "H2"  }
                }
             },
             { "HOH", 3,
                {
                   { "O3",   "O"   },
                   { "HO31", "H1"  },
                   { "HO32", "H2"  }
                }
             }
          }
       }
};

static std::map<std::string, int> _splitIndex;

#endif
