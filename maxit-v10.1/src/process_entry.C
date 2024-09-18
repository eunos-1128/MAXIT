/*
FILE:     process_entry.C
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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "BondUtil.h"
#include "CategoryMapping.h"
#include "ConnectDic.h"
#include "Element.h"
#include "Maxit.h"
#include "NdbToken.h"
#include "SeqCodeUtil.h"
#include "SgCenter.h"
#include "SpaceGroup.h"
#include "SplitUtil.h"
#include "TypeDef.h"

static void usage(const std::string& program, LogUtil& logutil);

int main(int argc, char **argv)
{
       std::string logfile = argv[0];
       logfile += ".log";
       std::string::size_type p = logfile.find_last_of('/');
       if (p != std::string::npos) logfile = logfile.substr(p + 1);
       LogUtil logutil(logfile);

       char* rcsbroot = getenv("RCSBROOT");

       if (rcsbroot == NULL) {
            logutil.messageWarning("Environment variable 'RCSBROOT' not defined");
            fprintf(stdout, "Environment variable 'RCSBROOT' not defined\n");
            return 0;
       }

       std::string inputfile, outputfile, user_logfile;
       inputfile.clear();
       outputfile.clear();
       user_logfile.clear();

       int input_format = 0;
       int output_format = 0;
       bool skip_flag = false;

       for (int i = 1; i < argc; i++) {
            if (argv[i][0] == '-' ) {
                 if (!strcasecmp(argv[i], "-i") || !strcasecmp(argv[i], "-input")) {
                      i++;
                      if (i < argc) inputfile = argv[i];
                      else usage(argv[0], logutil);
                 } else if (!strcasecmp(argv[i], "-input_format")) {
                      i++;
                      if (i < argc) {
                           if (!strcasecmp(argv[i], "pdb"))
                                input_format = NDB_FILE_FORMAT_PDB;
                           else if (!strcasecmp(argv[i], "cif"))
                                input_format = NDB_FILE_FORMAT_MMCIF;
                      } else usage(argv[0], logutil);
                 } else if (!strcasecmp(argv[i], "-o") || !strcasecmp(argv[i], "-output")) {
                      i++;
                      if (i < argc) outputfile = argv[i];
                      else usage(argv[0], logutil);
                 } else if (!strcasecmp(argv[i], "-output_format")) {
                      i++;
                      if (i < argc) {
                           if (!strcasecmp(argv[i], "pdb"))
                                output_format = NDB_FILE_FORMAT_PDB;
                           else if (!strcasecmp(argv[i], "cif"))
                                output_format = NDB_FILE_FORMAT_MMCIF;
                      } else usage(argv[0], logutil);
                 } else if (!strcasecmp(argv[i], "-log")) {
                      i++;
                      if (i < argc) user_logfile = argv[i];
                      else usage(argv[0], logutil);
                 } else if (!strcasecmp(argv[i], "-keep_original_numbering")) {
                      skip_flag = true;
                 } else {
                      usage(argv[0], logutil);
                      logutil.close();
                      return (1);
                 }

            } else {
                 usage(argv[0], logutil);
                 logutil.close();
                 return (1);
            }
       }

       if (inputfile.empty() || outputfile.empty() || (input_format == 0) || (output_format == 0)) {
            usage(argv[0], logutil);
            logutil.close();
            return (1);
       }

       if (!user_logfile.empty()) {
            logutil.close();
            logutil.open(user_logfile);
       }

       struct stat statbuf;
       if (stat(inputfile.c_str(), &statbuf) < 0) {
            std::string message = "File " + inputfile + " does not exist.\n";
            logutil.messageWarning(message);
            fprintf(stdout, "%s\n", message.c_str());
            logutil.close();
            return (1);
       }

       try {
            BondUtil::initialize();
            Element::initialize();
            NdbToken::initialize();
            CategoryMapping::initialize();
            SeqCodeUtil::initialize();
            SgCenter::initialize();
            SpaceGroup::initialize();
            SplitUtil::initialize();

            NdbToken::Read(logutil, rcsbroot);
            SgCenter::Read(rcsbroot);
            SpaceGroup::Read(logutil, rcsbroot);
            CategoryMapping::Read(logutil, rcsbroot);
            
            ConnectDic ccDic;
            ccDic.setRCSBROOT(rcsbroot);
            if (!ccDic.OpenFile()) {
                 logutil.messageError("Read chemical component dictionary failed.");
                 logutil.close();
                 return (-1);
            }

            Maxit obj;
            obj.setLog(&logutil);
            obj.setCCDic(&ccDic);
            obj.setRCSBROOT(rcsbroot);
            obj.set_input_filename(inputfile);
            obj.set_input_format(input_format);
            obj.set_output_filename(outputfile);
            obj.set_output_format(output_format);
            obj.set_include_date_original_flag();

            if (input_format == NDB_FILE_FORMAT_PDB)
                 obj.read_pdb_file();
            else if (input_format == NDB_FILE_FORMAT_MMCIF)
                 obj.read_mmcif_file();

            obj.build_molecule();

            if (input_format == NDB_FILE_FORMAT_PDB)
                 obj.Read_PDB_StructuralFeatures();
            else if (input_format == NDB_FILE_FORMAT_MMCIF)
                 obj.ReadStructuralFeatures();

            obj.CorrectAtomName(true);
            obj.SwitchLabeling();
            // moving water & calculate distance water
            obj.Moving_Waters(skip_flag);
            obj.CheckDistantWaters();
            // calculate link, ssbond & base pair information
            obj.Calculate_Link_and_SSBond(FIND_LINK | FIND_SSBOND, true);
            obj.CalculateBasePairInfo(true);
            obj.Update_Base_Pair_and_Step();
            // calculate secondary structures
            obj.CalculateSecondStruct("", true);
            obj.CheckGeometry();
            obj.CheckChiralityAndPlanarity();
            obj.calculate_all_contact();
            obj.AssignEntityId();
            obj.AssignAsymId();

            if (input_format == NDB_FILE_FORMAT_PDB) {
                 obj.pdb_to_ndb();
                 obj.ndb_to_cif();
            }

            // calculate cis peptide bond
            obj.UpdateStructMonProtCis();
            obj.WriteStructuralFeatures(true, true, false);

            if (output_format == NDB_FILE_FORMAT_PDB) {
                 obj.clear_pdb_records();
                 obj.cif_to_ndb();
                 obj.ndb_to_pdb();
                 obj.write_pdb_file(false, true);
            } else if (output_format == NDB_FILE_FORMAT_MMCIF) {
                 obj.UpdateCoordinateCategory();
                 obj.Output_public_cif_only();
                 obj.write_mmcif_file(false, true, false);
            }

            logutil.finishMessage();
       } catch (const std::exception& exc) {
            std::cerr << exc.what();
            logutil.messageError(exc.what());
       }

       logutil.close();

       return 0;
}

static void usage(const std::string& program, LogUtil& logutil)
{
       std::string log = "Usage: " + program + " -input inputfile -input_format pdb/cif -output outputfile -output_format pdb/cif ";
       log += "[ -log logfile -keep_original_numbering ]\n";
       logutil.messageWarning(log);
       fprintf(stdout, "%s\n", log.c_str());
}
