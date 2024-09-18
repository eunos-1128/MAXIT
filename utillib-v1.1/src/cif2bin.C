/*
FILE:     cif2bin.C
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
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "CifFile.h"
#include "DICParserBase.h"

#define NUM_FILES    20

#include "utillib.h"

typedef struct {
        const char *textfile;
        const char *binfile;
} FILE_PAIR;

static FILE_PAIR file_pair[NUM_FILES] =
{
       { "batch_mandatory_items.cif",       "batch_mandatory_items.odb"       },
       { "corres_info_items.cif",           "corres_info_items.odb"           },
       { "dict_config.cif",                 "dict_config.odb"                 },
       { "em_v4_v5_mapping.cif",            "em_v4_v5_mapping.odb"            },
       { "master_internal_cif_mapping.cif", "master_internal_cif_mapping.odb" },
       { "merge_meta_data_tasks.cif",       "merge_meta_data_tasks.odb"       },
       { "nmr_star_update_items.cif",       "nmr_star_update_items.odb"       },
       { "nmr_v4_v5_mapping.cif",           "nmr_v4_v5_mapping.odb"           },
       { "ndb_cif.cif",                     "ndb_cif.odb"                     },
       { "ndb_cif_rms.cif",                 "ndb_cif_rms.odb"                 },
       { "ndb_refn.cif",                    "ndb_refn.odb"                    },
       { "ndb_tokens.cif",                  "ndb_tokens.odb"                  },
       { "old_new_version_mapping.cif",     "old_new_version_mapping.odb"     },
       { "pdb_rev_record_types.cif",        "pdb_rev_record_types.odb"        },
       { "refinement.cif",                  "refinement.odb"                  },
       { "sg_center.cif",                   "sg_center.odb"                   },
       { "space_group.cif",                 "space_group.odb"                 },
       { "v4_v5_clean_up.cif",              "v4_v5_clean_up.odb"              },
       { "version_content_types.cif",       "version_content_types.odb"       },
       { "xray_v4_v5_mapping.cif",          "xray_v4_v5_mapping.odb"          }
};

static char defalut_dir[] = "/ndbapps";

int main(int argc, char **argv)
{
       char *root_dir = getenv("RCSBROOT");
       if (root_dir == NULL) root_dir = defalut_dir;

       for (int i = 1; i < argc; i++) {
            if (argv[i][0] == '-' ) {
                 if (!strcasecmp(argv[i], "-path")) {
                      i++; root_dir = argv[i];
                 }
            }
       }

       struct stat stattext, statbin;
       std::string textfile, binaryfile, diags;
       for (int i = 0; i < NUM_FILES; i++) {
            bool create_flag = false;
            textfile = root_dir; textfile += "/data/ascii/";
            textfile += file_pair[i].textfile;
            binaryfile = root_dir; binaryfile += "/data/binary/";
            binaryfile += file_pair[i].binfile;
            if (stat(textfile.c_str(), &stattext) == 0) {
                 create_flag = true;
                 if (stat(binaryfile.c_str(), &statbin) == 0) {
                      if (statbin.st_mtime > stattext.st_mtime) create_flag = false;
                 }
            }
            if (!create_flag) continue;

            CifFile *fobj = read_cif_file(binaryfile, textfile, diags);
            if (!fobj) {
                 if (!diags.empty()) {
                      fprintf(stdout, "Read cif file %s failed\n%s\n", textfile.c_str(),
                              diags.c_str());
                      exit (-1);
                 }
            }
            fobj->Close();
            delete fobj;
            chmod(binaryfile.c_str(), 00644);
       }

       return 0;
}
