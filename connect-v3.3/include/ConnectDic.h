/*
FILE:     ConnectDic.h
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

/*!
** \file ConnectDic.h
**
** \brief Header file for ConnectDic class.
*/

#ifndef _H_CONNECTDIC_H_
#define _H_CONNECTDIC_H_

#include <map>
#include <set>
#include <string>
#include <vector>

#include "CifFile.h"
#include "ConnectFormat.h"

#define MONOMER_TYPE_UNK            0
#define MONOMER_TYPE_DNA            1
#define MONOMER_TYPE_RNA            2
#define MONOMER_TYPE_A_A            3

#define componenttextfile      "component.cif"
#define componentbinaryfile    "component.odb"
#define varianttextfile        "variant.cif"
#define variantbinaryfile      "variant.odb"

/**
**  \class ConnectDic
** 
**  \brief Public class that respresents a chemical component dictionary.
** 
**  The class provides methods for construction, destruction and retrieval.
*/

class ConnectDic {
   private:
       std::string _empty_string;
       std::set<std::string> _empty_set;
       std::vector<std::string> _empty_vector;

       // chemical components array
       std::vector<ConnectFormat> _drugs;
       // ID/index mapping for chemical components
       std::map<std::string, int> _ConnMap;
       // ID/index mapping for variant chemical components
       std::multimap<std::string, int> _VariantConnMap;
       // component.odb CifFile object
       CifFile *_fObj;
       // COMP_PATH path
       std::string _rcsbrootpath;
       // COMP_PATH path
       std::string _cvscomponentpath, _additional_comp_path;
       // peptide type set
       std::set<std::string> _peptide_type;
       // nucleic acid type set
       std::set<std::string> _nucleic_type;
       // DNA type set
       std::set<std::string> _dna_type;
       // RNA type set
       std::set<std::string> _rna_type;
       // ATOM type set
       std::set<std::string> _atom_type;

       std::map<std::string, std::set<std::string> > _tailTerminalAtoms, _headTerminalAtoms;

       // key: CCD ID
       // value: The user provided CCD file name with full path.
       std::map<std::string, std::string> _userProvidedCCDFiles;

       /**
       **  Clear data storage
       **
       **  \param: Not applicable
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void clear();

       /**
       **  Initialize type set
       **
       **  \param: Not applicable
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void init();

       /**
       **  Read chemical component from data block and add it to dictionary
       **
       **  \param[in]: block - reference to chemical component data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void add_drug(Block &block, const bool& with_preferred_name = false);

       /**
       **  Read chemical component file and add it to dictionary
       **
       **  \param[in]: drugID   - chemical component ID
       **  \param[in]: compfile - chemical component file name
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void add_drug(const std::string& drugID, const std::string& compfile, const bool& with_preferred_name = false);

       /**
       **  Find terminal atom set
       **
       **  \param[in]: drugID   - chemical component ID
       **  \param[in]: searchType - search atom type
       **  \param[in]: excludeType - exclude atom type
       **
       **  \return terminal atom set
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       std::set<std::string> _get_terminal_atoms(const std::string& drugID, const std::string& searchType, const std::string& excludeType);

   public:
       /**
       **  Constructs ConnectDic
       **
       **  \param: Not applicable
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       ConnectDic();

       /**
       **  Destructs ConnectDic
       **
       **  \param: Not applicable
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       ~ConnectDic();

       /**
       **  Set RCSBROOT path
       **
       **  \param[in]: path - RCSBROOT path 
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void setRCSBROOT(const std::string& path);

       /**
       **  Set COMP_PATH path
       **
       **  \param[in]: path - COMP_PATH path 
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void setCOMP_PATH(const std::string& path);

       /**
       **  Set additional COMP_PATH path
       **
       **  \param[in]: path - COMP_PATH path 
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void setAdditionalCompPath(const std::string& path);

       /**
       **  Get COMP_PATH path
       **
       **  \param: Not applicable
       **
       **  \return _cvscomponentpath
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string& getCOMP_PATH() const;

       /**
       **  Get additional COMP_PATH path
       **
       **  \param: Not applicable
       **
       **  \return _additional_comp_path
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string& getAdditionalCompPath() const;

       /**
       **  Read user defined chemical component dictionaries - Only read file name list and read the content later as needed.
       **                                                      The file name must be ${CCDID}.cif
       **
       **  \param[in]: cc_directory_path - The directory path which contains user defined chemical component dictionaries
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void readUserProvidedCCDictionary(const std::string& cc_directory_path);

       /**
       **  Read user defined chemical component dictionaries - Read all at once in the beginning
       **
       **  \param[in]: cc_directory_path - The directory path which contains user defined chemical component dictionaries
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void readUserProvidedCCDictionaryAll(const std::string& cc_directory_path);

       /**
       **  Read binary chemical component dictionaries
       **
       **  \param: Not applicable
       **
       **  \return true  - if successful
       **          false - if failed
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool OpenFile();

       /**
       **  Find chemical component based on ID
       **
       **  \param[in]: drugID   - chemical component ID
       **  \param[in]: compfile - chemical component file name
       **
       **  \return found chemical component 
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: NotFoundException - if chemical component with ID drugID does not exist
       */
       const ConnectFormat& find_drug(const std::string& drugID, const std::string& compfile = "");

       /**
       **  Find chemical component based on ID
       **
       **  \param[in]: drugID   - chemical component ID
       **
       **  \return found chemical component 
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: NotFoundException - if chemical component with ID drugID does not exist
       */
       const ConnectFormat& find_author_defined_drug(const std::string& drugID);

       /**
       **  Add chemical component based on ID
       **
       **  \param[in]: drugID   - Input chemical component ID
       **  \param[in]: drug     - Input chemical component definition
       **
       **  \return None
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void add_author_defined_drug(const std::string& drugID, const ConnectFormat& drug);

       /**
       **  Find residue's chain type based on ID
       **
       **  \param[in]: drugID - chemical component ID
       **
       **  \return found chain type
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string find_residue_type(const std::string& drugID); 

       /**
       **  Find residue's monomer type based on ID
       **
       **  \param[in]: drugID - chemical component ID
       **
       **  \return found monomer type
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const int find_monomer_type(const std::string& drugID);

       /**
       **  Find atom type for given chemical component ID & atom name
       **
       **  \param[in]: drugID - chemical component ID
       **  \param[in]: atomName - atom name
       **
       **  \return found atomn type
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string getAtomType(const std::string& drugID, const std::string& atomName);

       /**
       **  Find molecule weight for given chemical component ID
       **
       **  \param[in]: drugID - chemical component ID
       **
       **  \return found molecule weight
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const double get_molecule_weight(const std::string& drugID);

       /**
       **  Find if a chemical component is L-type amino acid
       **
       **  \param[in]: drugID - chemical component ID
       **
       **  \return true  - if it's L-type amino acid
       **          false - if it's not
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool is_L_aa_residue(const std::string& drugID);

       /**
       **  Find if a chemical component is D-type amino acid
       **
       **  \param[in]: drugID - chemical component ID
       **
       **  \return true  - if it's D-type amino acid
       **          false - if it's not
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool is_D_aa_residue(const std::string& drugID);

       /**
       **  Find Chemical Component's Meta data based on ID & cif item name
       **
       **  \param[in]: residueName - chemical component ID
       **  \param[in]: cifitem     - cif item name
       **
       **  \return found Meta data
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string find_cc_metadata(const std::string& residueName, const std::string& cifitem);

       /**
       **  Find first terminal atom based on ID & atom type
       **
       **  \param[in]: residueName - chemical component ID
       **  \param[in]: type        - atom type
       **
       **  \return found first terminal atom name
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string& find_terminal_atom(const std::string& residueName, const std::string& type);

       /**
       **  Find terminal atom(s) based on ID & atom type
       **
       **  \param[in]: residueName - chemical component ID
       **  \param[in]: type        - atom type
       **
       **  \return found terminal atom name(s)
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::vector<std::string>& find_terminal_atoms(const std::string& residueName, const std::string& type);

       /**
       **  Find tail terminal atom set ( "C" for peptide polymer, "O3'" nucleic acid polymer )
       **
       **  \param[in]: type        - polymer type (ATOMP/ATOMN)
       **  \param[in]: residueName - chemical component ID
       **  \param[out]: terminalAtoms - terminal atom set
       **
       **  \return standard terminal atom ( "C"/"O3'" )
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       std::string getTailTerminalAtoms(const std::string& type, const std::string& residueName, std::set<std::string>& terminalAtoms);

       /**
       **  Find head terminal atom set ( "N" for peptide polymer, "P" nucleic acid polymer )
       **
       **  \param[in]: type        - polymer type (ATOMP/ATOMN)
       **  \param[in]: residueName - chemical component ID
       **  \param[out]: terminalAtoms - terminal atom set
       **
       **  \return standard terminal atom ( "N"/"P" )
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       std::string getHeadTerminalAtoms(const std::string& type, const std::string& residueName, std::set<std::string>& terminalAtoms);

       /**
       **  Check if the given atom is terminal atom for given component ID
       **
       **  \param[in]: residueName - chemical component ID
       **  \param[in]: atomName    - atom name
       **
       **  \return true/false
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const bool is_terminal_atom(const std::string& residueName, const std::string& atomName);

       /**
       **  Find all chemical components related to given ID, includes variants
       **
       **  \param[in]:  drugID    - chemical component ID
       **  \param[out]: drug_list - list of found chemical components
       **
       **  \return None
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception None
       */
       void find_drugs(const std::string& drugID, std::vector<ConnectFormat>& drug_list);

       /**
       **  Find all chemical component variants related to given ID
       **
       **  \param[in]:  drugID    - chemical component ID
       **  \param[out]: drug_list - list of found chemical component variants
       **
       **  \return None
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception None
       */
       void find_variants(const std::string& drugID, std::vector<ConnectFormat>& drug_list);

       /**
       **  Find all atoms with leaving atoms attached
       **
       **  \param[in]:  drugID          - chemical component ID
       **
       **  \return found atom set
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception None
       */
       const std::set<std::string>& get_connected_atoms(const std::string& drugID);

       /**
       **  Convert textual chemical component dictionary into binary
       **
       **  \param[in]: textfile   - textual dictionary filename
       **  \param[in]: binaryfile - binary dictionary filename
       **
       **  \return None
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception None
       */
       bool BinaryConvertor(const std::string& textfile, const std::string& binaryfile);

       /**
       **  Find if a nucleic acid component has valid phosphoryl group 
       **
       **  \param[in]: drugID - chemical component ID
       **
       **  \return true  - if it has valid phosphoryl group
       **          false - if it has not
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool hasValidPhosphorylGroup(const std::string& drugID);
};

#endif
