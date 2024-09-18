/*
FILE:     CifCorrector.C
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


#include <string>
#include <iostream>

#include "GenCont.h"
#include "RcsbFile.h"
#include "CifFile.h"
#include "DataInfo.h"
#include "CifCorrector.h"


using std::string;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;


CifCorrector::CifCorrector(CifFile& cifFile, DataInfo& dataInfo,
  DataInfo& pdbxDataInfo, CifFile& configFile, const bool verbose) :
  _cifFile(cifFile), _dataInfo(dataInfo), _pdbxDataInfo(pdbxDataInfo),
  _configFile(configFile), _verbose(verbose)
{
    Block& block = _configFile.GetBlock("cif_corrector");

    _configTableP = block.GetTablePtr("config");

    if (_configTableP == NULL)
    {
        throw EmptyValueException("No configuration table defined.",
          "CifCorrector::CifCorrector");
    }

    ValidateConfigTable();
}


CifCorrector::~CifCorrector()
{

}


void CifCorrector::MakeOutputCifFileName(string& outCifFileName,
  const string& inCifFileName)
{
    string relInFileName;
    RcsbFile::RelativeFileName(relInFileName, inCifFileName);

    outCifFileName = relInFileName + ".corrected";
}


CifFile* CifCorrector::CreateConfigFile()
{
    CifFile* configFileP = new CifFile();

    configFileP->AddBlock("cif_corrector");

    Block& block = configFileP->GetBlock(configFileP->GetFirstBlockName());

    ISTable* configTableP = new ISTable("config");
    configTableP->AddColumn("oper");
    configTableP->AddColumn("item");
    configTableP->AddColumn("item_value");
    configTableP->AddColumn("ref_item");
    configTableP->AddColumn("ref_item_value");


    // Note: Order of operations matter. Item names related operations
    // must come first, followed by value related operations in the order
    // of complexity (e.g., uppercasing before complex value changes)

    vector<string> configRow(configTableP->GetNumColumns());

    configRow.assign(configRow.size(), CifString::UnknownValue);

    // Set the items that are to be removed
    configRow[0] = "remove";

    configRow[1] = "_diffrn_source.pdbx_wavelength_list";
    configTableP->AddRow(configRow);


    configRow.assign(configRow.size(), CifString::UnknownValue);

    // Set the items that are to be renamed
    configRow[0] = "rename";

    configRow[1] = "_struct_asym.ndb_pdb_id";
    configRow[3] = "_struct_asym.ndb_PDB_id";
    configTableP->AddRow(configRow);


    configRow.assign(configRow.size(), CifString::UnknownValue);

    // Set the items that are to be uppercased. This must be the first
    // operation.

    configRow[0] = "upper_case";

    configRow[1] = "_pdbx_entry_details.compound_details";
    configTableP->AddRow(configRow);

    configRow[1] = "_pdbx_entry_details.source_details";
    configTableP->AddRow(configRow);

    configRow[1] = "_pdbx_entry_details.nonpolymer_details";
    configTableP->AddRow(configRow);

    configRow[1] = "_pdbx_entry_details.sequence_details";
    configTableP->AddRow(configRow);

    configRow[1] = "_pdbx_database_status.dep_release_code_coordinates";
    configTableP->AddRow(configRow);

    configRow[1] = "_pdbx_database_status.dep_release_code_struct_fact";
    configTableP->AddRow(configRow);

    configRow[1] = "_pdbx_database_status.dep_release_code_sequence";
    configTableP->AddRow(configRow);

    configRow[1] = "_ndb_database_status.dep_release_code_struct_fact";
    configTableP->AddRow(configRow);

    configRow[1] = "_ndb_database_status.dep_release_code_coordinates";
    configTableP->AddRow(configRow);


    // Set the items that are to be value changed
    configRow.assign(configRow.size(), CifString::UnknownValue);

    configRow[0] = "value_change";

    configRow[1] = "_pdbx_entity_name.name_type";
    configRow[2] = CifString::InapplicableValue;
    configRow[4] = "RCSB_NAME";
    configTableP->AddRow(configRow);

    configRow[1] = "_em_assembly.aggregation_state";
    configRow[2] = "single particle";
    configRow[4] = "PARTICLE";
    configTableP->AddRow(configRow);

    configRow[1] = "_diffrn_radiation.rcsb_diffrn_protocol";
    configRow[2] = "SAD";
    configRow[4] = "SINGLE WAVELENGTH";
    configTableP->AddRow(configRow);

    configRow[1] = "_ndb_database_status.dep_release_code_sequence";
    configRow[2] = "HOLD FOR PUBLICATION";
    configRow[4] = "HOLD FOR RELEASE";
    configTableP->AddRow(configRow);

    configRow[1] = "_pdbx_database_status.dep_release_code_coordinates";
    configRow[2] = "HOLD FOR RELEASE";
    configRow[4] = "HOLD FOR PUBLICATION";
    configTableP->AddRow(configRow);

    configRow[1] = "_pdbx_database_status.dep_release_code_struct_fact";
    configRow[2] = "HOLD FOR RELEASE";
    configRow[4] = "HOLD FOR PUBLICATION";
    configTableP->AddRow(configRow);

    configRow[1] = "_ndb_database_status.dep_release_code_struct_fact";
    configRow[2] = "HOLD FOR RELEASE";
    configRow[4] = "HOLD FOR PUBLICATION";
    configTableP->AddRow(configRow);

    configRow[1] = "_ndb_database_status.dep_release_code_coordinates";
    configRow[2] = "HOLD FOR RELEASE";
    configRow[4] = "HOLD FOR PUBLICATION";
    configTableP->AddRow(configRow);


    // Set the items which numeric list values are to be fixed
    configRow.assign(configRow.size(), CifString::UnknownValue);

    configRow[0] = "numeric_list";
    configRow[1] = "_diffrn_source.rcsb_wavelength_list";
    configTableP->AddRow(configRow);


    // Set the items which missing values are to be fixed
    configRow.assign(configRow.size(), CifString::UnknownValue);

    configRow[0] = "missing_values";

    configRow[1] = "_refine_ls_shell.d_res_high";
    configRow[3] = "_reflns_shell.d_res_high";
    configTableP->AddRow(configRow);


    // Set the items which labeling are to be fixed
    configRow.assign(configRow.size(), CifString::UnknownValue);

    configRow[0] = "labeling";

    configRow[1] = "_refine_ls_restr_ncs.dom_id";
    configRow[3] = "_struct_ncs_dom.id";
    configTableP->AddRow(configRow);


    // Set the items which bad sequence are to be fixed
    configRow.assign(configRow.size(), CifString::UnknownValue);

    configRow[0] = "value_change_complex";

    configRow[1] = "_refine_ls_restr_ncs.pdbx_asym_id";
    configRow[3] =
      "_refine_ls_restr_ncs.pdbx_ens_id,_refine_ls_restr_ncs.dom_id;"\
      "_refine_ls_restr_ncs.pdbx_ens_id,_refine_ls_restr_ncs.dom_id|"\
      "_struct_ncs_dom.pdbx_ens_id,_struct_ncs_dom.id;"\
      "_struct_ncs_dom.details|"\
      "_struct_asym.ndb_PDB_id,_struct_asym.ndb_type;_struct_asym.id";

    configRow[4] =
      ".,.;"\
      ".,.;"\
      ".,ATOMP";

    configTableP->AddRow(configRow);

    block.WriteTable(configTableP);

    return (configFileP);
} // End of CifCorrector::CreateConfigFile()


void CifCorrector::Correct()
{
    for (unsigned int confRowI = 0; confRowI < _configTableP->GetNumRows();
      ++confRowI)
    {
        const string& oper = (*_configTableP)(confRowI, "oper");
        const string& item = (*_configTableP)(confRowI, "item");
        const string& itemValue = (*_configTableP)(confRowI, "item_value");
        const string& refItem = (*_configTableP)(confRowI, "ref_item");
        const string& refItemValue = (*_configTableP)(confRowI,
          "ref_item_value");

        if (oper == "upper_case")
            CorrectUpperCase(item);
        else if (oper == "rename")
            RenameItem(item, refItem);
        else if (oper == "remove")
            RemoveItem(item);
        else if (oper == "value_change")
            CorrectValues(item, itemValue, refItemValue);
        else if (oper == "numeric_list")
            CorrectNumericList(item);
        else if (oper == "missing_values")
            CorrectMissingValues(item, refItem);
        else if (oper == "labeling")
            CorrectLabeling(item, refItem);
        else if (oper == "value_change_complex")
            CorrectBadSequence(item, refItem, refItemValue);
        else
            cerr << "Warning: Bad operation \"" << oper << "\" in row# " <<
              confRowI + 1 << " of table \"" << _configTableP->GetName() <<
              "\"" << endl;
    }

    CorrectNotApplicableValues();
    CorrectEnums();

} // End of CifCorrector::Correct()


void CifCorrector::CheckAliases()
{
    // Get list of all categories in the CIF file
    Block& block = _cifFile.GetBlock(_cifFile.GetFirstBlockName());

    vector<string> catNames;
    block.GetTableNames(catNames);

    for (unsigned int catI = 0; catI < catNames.size(); ++catI)
    {
        // Get category table pointer.
        ISTable* catTableP = block.GetTablePtr(catNames[catI]);

        // Get category attributes.
        const vector<string>& attrNames = catTableP->GetColumnNames();

        for (unsigned int attrI = 0; attrI < attrNames.size(); ++attrI)
        {
            // Make a CIF item.
            string item;
            CifString::MakeCifItem(item, catNames[catI], attrNames[attrI]);

            // Is it in internal? 
            if (!_dataInfo.IsItemDefined(item))
            {
                cerr << "Warning: Item \"" << item << "\" is not defined in "\
                  "the internal dictionary." << endl;

                // Is it in PDBx dictionary?
                if (_pdbxDataInfo.IsItemDefined(item))
                {
                    const vector<string>& aliases =
                      _dataInfo.GetItemAttribute(item,
                      CifString::CIF_DDL_CATEGORY_ITEM_ENUMERATION,
                      CifString::CIF_DDL_ITEM_VALUE);

                    if (aliases.empty() || aliases[0].empty())
                    {
                        // No alias.
                        cerr << "Warning: Item \"" << item << "\" is "\
                          "defined in the PDBx dictionary but does not "\
                          "have an alias." << endl;
                    }
                }
                else
                {
                    // Unknown item.
                    cerr << "Warning: Item \"" << item << "\" is not "\
                      "defined either in the internal or in the PDBx "\
                      "dictionary." << endl;
                }
            } // if (item not in internal dictionary) 
        } // for (all attributes)
    } // for (all categories)
}


void CifCorrector::CorrectUpperCase(const string& item)
{
    if (!_dataInfo.IsItemDefined(item))
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Item \"" << item << "\" is not defined in the "\
          "internal dictionary." << endl;
#endif
    }

    // From item name, get category name
    string catName;
    CifString::GetCategoryFromCifItem(catName, item);

    Block& block = _cifFile.GetBlock(_cifFile.GetFirstBlockName());

    // Get a pointer to a category table
    ISTable* catTableP = block.GetTablePtr(catName);
    if (catTableP == NULL)
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Table \"" << catName << "\" not found." << endl;
#endif
        return;
    }

    // From item name, get attribute name
    string attribName;
    CifString::GetItemFromCifItem(attribName, item);

    // Check to see if column is present
    if (!catTableP->IsColumnPresent(attribName))
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Table \"" << catName << "\" does not "\
          "have column \"" << attribName << "\"." << endl;
#endif
        return;
    }

    for (unsigned int rowI = 0; rowI < catTableP->GetNumRows(); ++rowI)
    {
        string cellValue = (*catTableP)(rowI, attribName);

        if (CifString::IsEmptyValue(cellValue))
        {
            continue;
        }

        String::UpperCase(cellValue);

        catTableP->UpdateCell(rowI, attribName, cellValue);

        if (_verbose)
        {
            cout << "Info: Uppercasing item \"" << item <<
              "\" to make it \"" << cellValue << "\"" << endl;
        }
    } // for (all rows in category table)
} // End of CifCorrector::CorrectUpperCase()


void CifCorrector::RenameItem(const string& item, const string& refItem)
{
    if (!_dataInfo.IsItemDefined(item))
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Item \"" << item << "\" is not defined in the "\
          "internal dictionary." << endl;
#endif
    }

    if (!_dataInfo.IsItemDefined(refItem))
    {
        cerr << "Warning: Item \"" << refItem << "\" is not defined in the "\
          "internal dictionary." << endl;
        return;
    }

    // From item name, get category name
    string catName;
    CifString::GetCategoryFromCifItem(catName, item);

    // From reference item name, get category name
    string refCatName;
    CifString::GetCategoryFromCifItem(refCatName, refItem);

    if (catName != refCatName)
    {
        cerr << "Warning: Cannot rename item \"" << item << "\" to item \"" <<
          refItem << "\", since they belong to different categories." << endl;
        return;
    }

    Block& block = _cifFile.GetBlock(_cifFile.GetFirstBlockName());

    // Get a pointer to a category table
    ISTable* catTableP = block.GetTablePtr(catName);
    if (catTableP == NULL)
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Table \"" << catName << "\" not found." << endl;
#endif
        return;
    }

    // From item name, get attribute name
    string attribName;
    CifString::GetItemFromCifItem(attribName, item);

    // Check to see if column is present
    if (!catTableP->IsColumnPresent(attribName))
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Table \"" << catName << "\" does not "\
          "have column \"" << attribName << "\"." << endl;
#endif
        return;
    }

    string refAttribName;
    CifString::GetItemFromCifItem(refAttribName, refItem);

    // Check to see if column is present
    if (catTableP->IsColumnPresent(refAttribName))
    {
        cerr << "Warning: Table \"" << catName << "\" already has "\
          "column \"" << refAttribName << "\"." << endl;
        return;
    }

    if (_verbose)
    {
        cout << "Info: Renaming item \"" << item << "\"";
    }

    catTableP->RenameColumn(attribName, refAttribName);

    if (_verbose)
    {
        cout << " to \"" << refItem << "\"" << endl;
    }
} // End of CifCorrector::RenameItem()


void CifCorrector::RemoveItem(const string& item)
{
    if (!_dataInfo.IsItemDefined(item))
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Item \"" << item << "\" is not defined in the "\
          "internal dictionary." << endl;
#endif
    }

    // From item name, get category name
    string catName;
    CifString::GetCategoryFromCifItem(catName, item);

    Block& block = _cifFile.GetBlock(_cifFile.GetFirstBlockName());

    // Get a pointer to a category table
    ISTable* catTableP = block.GetTablePtr(catName);
    if (catTableP == NULL)
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Table \"" << catName << "\" not found." << endl;
#endif
        return;
    }

    // From item name, get attribute name
    string attribName;
    CifString::GetItemFromCifItem(attribName, item);

    // Check to see if column is present
    if (!catTableP->IsColumnPresent(attribName))
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Table \"" << catName << "\" does not "\
          "have column \"" << attribName << "\"." << endl;
#endif
        return;
    }

    catTableP->DeleteColumn(attribName);

    if (_verbose)
    {
        cout << "Info: Deleting item \"" << item << "\"" << endl;
    }
} // End of CifCorrector::RemoveItem()


void CifCorrector::CorrectValues(const string& item, const string& itemValue,
  const string& refItemValue)
{
    if (!_dataInfo.IsItemDefined(item))
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Item \"" << item << "\" is not defined in the "\
          "internal dictionary." << endl;
#endif
    }

    // From item name, get category name
    string catName;
    CifString::GetCategoryFromCifItem(catName, item);

    Block& block = _cifFile.GetBlock(_cifFile.GetFirstBlockName());

    // Get a pointer to a category table
    ISTable* catTableP = block.GetTablePtr(catName);
    if (catTableP == NULL)
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Table \"" << catName << "\" not found." << endl;
#endif
        return;
    }

    // From item name, get attribute name
    string attribName;
    CifString::GetItemFromCifItem(attribName, item);

    // Check to see if column is present
    if (!catTableP->IsColumnPresent(attribName))
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Table \"" << catName << "\" does not "\
          "have column \"" << attribName << "\"." << endl;
#endif
        return;
    }

    for (unsigned int rowI = 0; rowI < catTableP->GetNumRows(); ++rowI)
    {
        const string& currCellValue = (*catTableP)(rowI, attribName);

        if (currCellValue != itemValue)
        {
            continue;
        }

        catTableP->UpdateCell(rowI, attribName, refItemValue);

        if (_verbose)
        {
            cout << "Info: Changing item \"" << item <<
              "\" value from \"" << itemValue << "\" to \"" <<
              refItemValue << "\"" << endl;
        }
    } // for (all rows in category table)
} // End of CifCorrector::CorrectValues()


void CifCorrector::CorrectEnums()
{
    // Get list of all categories in the CIF file
    Block& block = _cifFile.GetBlock(_cifFile.GetFirstBlockName());

    vector<string> catNames;
    block.GetTableNames(catNames);

    for (unsigned int catI = 0; catI < catNames.size(); ++catI)
    {
        const string& catName = catNames[catI];

        // Get category table pointer.
        ISTable* catTableP = block.GetTablePtr(catName);

        // Get category attributes.
        const vector<string>& attrNames = catTableP->GetColumnNames();

        for (unsigned int attrI = 0; attrI < attrNames.size(); ++attrI)
        {
            const string& attribName = attrNames[attrI];

            // Make a CIF item.
            string item;
            CifString::MakeCifItem(item, catName, attribName);

            if (!_dataInfo.IsItemDefined(item))
            {
#ifdef VLAD_DELETE
                cerr << "Warning: Item \"" << item << "\" is not defined in "\
                  "the internal dictionary." << endl;
#endif
            }

            // Check to see if column is present
            if (!catTableP->IsColumnPresent(attribName))
            {
#ifdef VLAD_DELETED
                cerr << "Warning: Table \"" << catName << "\" does not "\
                  "have column \"" << attribName << "\"." << endl;
#endif
                continue;
            }

            vector<string> enums = _dataInfo.GetItemAttribute(item,
              CifString::CIF_DDL_CATEGORY_ITEM_ENUMERATION,
              CifString::CIF_DDL_ITEM_VALUE);

            if (enums.empty())
            {
                enums = _pdbxDataInfo.GetItemAttribute(item,
                  CifString::CIF_DDL_CATEGORY_ITEM_ENUMERATION,
                  CifString::CIF_DDL_ITEM_VALUE);
            }

            // Check if the attribute is defined as an enumeration
            if (enums.empty())
            {
#ifdef VLAD_DELETE
                cerr << "Warning: Item \"" << item << "\" does not have "\
                  "enumerations defined." << endl;
#endif
                continue;
            }

            for (unsigned int rowI = 0; rowI < catTableP->GetNumRows(); ++rowI)
            {
                const string& currCellValue = (*catTableP)(rowI, attribName);

                // Find what enumeration it is to be replaced with by finding
                // enumInd within enumerations.

                unsigned int enumInd =
                  CifCorrector::FindEnumIndex(currCellValue, enums);

                if (enumInd == enums.size())
                {
                    // No enumerations fix is done
                    continue;
                }

                if (currCellValue == enums[enumInd])
                {
                    continue;
                }

                if (_verbose)
                {
                    cout << "Info: Changing item \"" << item <<
                      "\" value from \"" << currCellValue << "\"";
                }

                catTableP->UpdateCell(rowI, attribName, enums[enumInd]);

                if (_verbose)
                {
                    cout << " to \"" << enums[enumInd] << "\"" << endl;
                }
            } // for (all rows in category table)
        } // for (all attributes of a category)
    } // for (all categories in CIF file)
} // End of CifCorrector::CorrectEnums()


void CifCorrector::CorrectNumericList(const string& item)
{
    if (!_dataInfo.IsItemDefined(item))
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Item \"" << item << "\" is not defined in the "\
          "internal dictionary." << endl;
#endif
    }

    // From item name, get category name
    string catName;
    CifString::GetCategoryFromCifItem(catName, item);

    Block& block = _cifFile.GetBlock(_cifFile.GetFirstBlockName());

    // Get a pointer to a category table
    ISTable* catTableP = block.GetTablePtr(catName);
    if (catTableP == NULL)
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Table \"" << catName << "\" not found." << endl;
#endif
        return;
    }

    // From item name, get attribute name
    string attribName;
    CifString::GetItemFromCifItem(attribName, item);

    // Check to see if column is present
    if (!catTableP->IsColumnPresent(attribName))
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Table \"" << catName << "\" does not have "\
          "column \"" << attribName << "\"." << endl;
#endif
        return;
    }

    for (unsigned int rowI = 0; rowI < catTableP->GetNumRows(); ++rowI)
    {
        const string& currCellValue = (*catTableP)(rowI, attribName);

        string fixedCellValue;
        FixNumericList(fixedCellValue, currCellValue);

        if (fixedCellValue == currCellValue)
        {
            // No fix is done
            continue;
        }

        if (_verbose)
        {
            cout << "Info: Changing item \"" << item <<
              "\" value from \"" << currCellValue << "\"";
        }

        catTableP->UpdateCell(rowI, attribName, fixedCellValue);

        if (_verbose)
        {
            cout << " to \"" << fixedCellValue << "\"" << endl;
        }
    } // for (all rows in category table)
} // End of CifCorrector::CorrectNumericList()


void CifCorrector::CorrectNotApplicableValues()
{
    // Get list of all categories in the CIF file
    Block& block = _cifFile.GetBlock(_cifFile.GetFirstBlockName());

    vector<string> catNames;
    block.GetTableNames(catNames);

    for (unsigned int catI = 0; catI < catNames.size(); ++catI)
    {
        const string& catName = catNames[catI];

        // Get category table pointer.
        ISTable* catTableP = block.GetTablePtr(catName);

        // Get category attributes.
        const vector<string>& attrNames = catTableP->GetColumnNames();

        for (unsigned int attrI = 0; attrI < attrNames.size(); ++attrI)
        {
            const string& attribName = attrNames[attrI];

            // Make a CIF item.
            string item;
            CifString::MakeCifItem(item, catName, attribName);

            if (!_dataInfo.IsItemDefined(item))
            {
#ifdef VLAD_DELETED
                cerr << "Warning: Item \"" << item << "\" is not defined in "\
                  "the internal dictionary." << endl;
#endif
            }

            // Check to see if column is present
            if (!catTableP->IsColumnPresent(attribName))
            {
#ifdef VLAD_DELETED
                cerr << "Warning: Table \"" << catName << "\" does not have "\
                  "column \"" << attribName << "\"." << endl;
#endif
                continue;
            }

            for (unsigned int rowI = 0; rowI < catTableP->GetNumRows(); ++rowI)
            {
                const string& currCellValue = (*catTableP)(rowI, attribName);

                string fixedCellValue;

                FixNotApplicable(fixedCellValue, currCellValue);

                if (fixedCellValue == currCellValue)
                {
                    // No fix is done
                    continue;
                }

                if (_verbose)
                {
                    cout << "Info: Changing item \"" << item <<
                      "\" value from \"" << currCellValue << "\"";
                }

                catTableP->UpdateCell(rowI, attribName, fixedCellValue);

                if (_verbose)
                {
                    cout << " to \"" << fixedCellValue << "\"" << endl;
                }
            } // for (all rows in category table)
        } // for (all attributes of a category)
    } // for (all categories in CIF file)
} // End of CifCorrector::CorrectNotApplicableValues()


void CifCorrector::CorrectMissingValues(const string& item,
  const string& refItem)
{
    if (!_dataInfo.IsItemDefined(item))
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Item \"" << item << "\" is not defined in the "\
          "internal dictionary." << endl;
#endif
    }

    if (!_dataInfo.IsItemDefined(refItem))
    {
        cerr << "Warning: Item \"" << refItem << "\" is not defined in the "\
          "internal dictionary." << endl;
        return;
    }

    string firstCat;
    CifString::GetCategoryFromCifItem(firstCat, item);

    Block& block = _cifFile.GetBlock(_cifFile.GetFirstBlockName());

    // Get a pointer to a category table
    ISTable* firstCatTableP = block.GetTablePtr(firstCat);
    if (firstCatTableP == NULL)
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Table \"" << firstCat << "\" not found." << endl;
#endif
        return;
    }

    string secondCat;
    CifString::GetCategoryFromCifItem(secondCat, refItem);

    // Get a pointer to a category table
    ISTable* secondCatTableP = block.GetTablePtr(secondCat);
    if (secondCatTableP == NULL)
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Table \"" << secondCat << "\" not found." << endl;
#endif
        return;
    }

    if ((firstCatTableP->GetNumRows() > 1) ||
      (secondCatTableP->GetNumRows() > 1))
    {
        cerr << "One of the tables has more than one row. Cannot fix "\
          "missing values." << endl;
        return;
    }

    // From item name, get attribute name
    string firstAttrib;
    CifString::GetItemFromCifItem(firstAttrib, item);

    // Check to see if column is present
    if (!firstCatTableP->IsColumnPresent(firstAttrib))
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Table \"" << firstCat << "\" does not have "\
          "column \"" << firstAttrib << "\"." << endl;
#endif
        return;
    }

    // From item name, get attribute name
    string secondAttrib;
    CifString::GetItemFromCifItem(secondAttrib, refItem);

    // Check to see if column is present
    if (!secondCatTableP->IsColumnPresent(secondAttrib))
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Table \"" << secondCat << "\" does not have "\
          "column \"" << secondAttrib << "\"." << endl;
#endif
        return;
    }

    const string& firstValue = (*firstCatTableP)(0, firstAttrib);
    const string& secondValue = (*secondCatTableP)(0, secondAttrib);

    bool firstValueUnknown = CifString::IsUnknownValue(firstValue);
    bool secondValueUnknown = CifString::IsUnknownValue(secondValue);

    if (firstValueUnknown && secondValueUnknown)
    {
        cerr << "Warning: Both values unknown. Cannot fix missing "\
          "values." << endl;
        return;
    }

    if (secondValueUnknown)
    {
        if (_verbose)
        {
            cout << "Info: Changing item \"" << refItem <<
              "\" value from \"" << secondValue << "\"";
        }

        secondCatTableP->UpdateCell(0, secondAttrib, firstValue);

        if (_verbose)
        {
            cout << " to \"" << firstValue << "\"" << endl;
        }
    }

    if (firstValueUnknown)
    {
        if (_verbose)
        {
            cout << "Info: Changing item \"" << item <<
              "\" value from \"" << firstValue << "\"";
        }

        firstCatTableP->UpdateCell(0, firstAttrib, secondValue);

        if (_verbose)
        {
            cout << " to \"" << secondValue << "\"" << endl;
        }
    }
} // End of CifCorrector::CorrectMissingValues()


void CifCorrector::CorrectLabeling(const string& item,
  const string& refItem)
{
    if (!_dataInfo.IsItemDefined(item))
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Item \"" << item << "\" is not defined in the "\
          "internal dictionary." << endl;
#endif
    }

    if (!_dataInfo.IsItemDefined(refItem))
    {
        cerr << "Warning: Item \"" << refItem << "\" is not defined in the "\
          "internal dictionary." << endl;
        return;
    }

    string catName;
    CifString::GetCategoryFromCifItem(catName, item);

    Block& block = _cifFile.GetBlock(_cifFile.GetFirstBlockName());

    // Get a pointer to a category table
    ISTable* catTableP = block.GetTablePtr(catName);
    if (catTableP == NULL)
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Table \"" << catName << "\" not found." << endl;
#endif
        return;
    }

    string refCat;
    CifString::GetCategoryFromCifItem(refCat, refItem);

    // Get a pointer to a category table
    ISTable* refCatTableP = block.GetTablePtr(refCat);
    if (refCatTableP == NULL)
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Table \"" << refCat << "\" not found." << endl;
#endif
        return;
    }

    // From item name, get attribute name
    string attribName;
    CifString::GetItemFromCifItem(attribName, item);

    // Check to see if column is present
    if (!catTableP->IsColumnPresent(attribName))
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Table \"" << catName << "\" does not have "\
          "column \"" << attribName << "\"." << endl;
#endif
        return;
    }

    // From item name, get attribute name
    string fromAttrib;
    CifString::GetItemFromCifItem(fromAttrib, refItem);

    // Check to see if column is present
    if (!refCatTableP->IsColumnPresent(fromAttrib))
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Table \"" << refCat << "\" does not have "\
          "column \"" << fromAttrib << "\"." << endl;
#endif
        return;
    }

    // VLAD-TODO Check if all values equal to 1
    for (unsigned int rowI = 0, fromRowI = 0; rowI < catTableP->GetNumRows();
      ++rowI)
    {
        const string& currCellValue = (*catTableP)(rowI, attribName);

        const string& itemValue = (*refCatTableP)(fromRowI, fromAttrib);

        if (currCellValue != itemValue)
        {
            if (_verbose)
            {
                cout << "Info: Changing item \"" << item <<
                  "\" value from \"" << currCellValue << "\"";
            }

            catTableP->UpdateCell(rowI, attribName, itemValue);

            if (_verbose)
            {
                cout << " to \"" << itemValue << "\"" << endl;
            }
        }

        ++fromRowI;
        if (fromRowI == refCatTableP->GetNumRows())
        {
            fromRowI = 0;
        }
    } // for (every row)
} // End of CifCorrector::CorrectLabeling()


void CifCorrector::CorrectBadSequence(const string& item,
  const string& refCondItem, const string& refCondItemValue)
{
    if (!_dataInfo.IsItemDefined(item))
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Item \"" << item << "\" is not defined in the "\
          "internal dictionary." << endl;
#endif
    }

    string catName;
    CifString::GetCategoryFromCifItem(catName, item);

    Block& block = _cifFile.GetBlock(_cifFile.GetFirstBlockName());

    // Get a pointer to a category table
    ISTable* catTableP = block.GetTablePtr(catName);
    if (catTableP == NULL)
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Table \"" << catName << "\" not found." << endl;
#endif
        return;
    }

    // From item name, get attribute name
    string attribName;
    CifString::GetItemFromCifItem(attribName, item);

    // Check to see if column is present
    if (!catTableP->IsColumnPresent(attribName))
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Table \"" << catName << "\" does not have "\
          "column \"" << attribName << "\"." << endl;
#endif
        return;
    }

    // Make this to operate on the whole column. Though you will not gain
    // performance
    vector<pair<vector<string>, vector<string> > > refMap;
    vector<vector<string> > refValues;

    CreateRefMap(refMap, refCondItem);
    ExtractRefValues(refValues, refCondItemValue);

    // IMPLEMENT
    // VerifyRefs(refMap, refValues);

    vector<string> fixedCellValues;
    FixBadSequence(fixedCellValues, refMap, refValues, block,
      catTableP->GetNumRows());

    if (fixedCellValues.size() != catTableP->GetNumRows())
    {
        // No fix is done. Something is wrong with references.
        return; 
    }

    for (unsigned int rowI = 0; rowI < catTableP->GetNumRows(); ++rowI)
    {
        const string& currCellValue = (*catTableP)(rowI, attribName);
        const string& fixedCellValue = fixedCellValues[rowI];

        if (fixedCellValue.empty())
        {
            // No fix is done. Something is wrong with references.
            continue; 
        }

        if (fixedCellValue == currCellValue)
        {
            // No fix is done
            continue;
        }

        if (_verbose)
        {
            cout << "Info: Changing item \"" << item <<
              "\" value from \"" << currCellValue << "\"";
        }

        catTableP->UpdateCell(rowI, attribName, fixedCellValue);

        if (_verbose)
        {
            cout << " to \"" << fixedCellValue << "\"" << endl;
        }
    } // for (all rows in category table)
} // End of CifCorrector::CorrectBadSequence()


void CifCorrector::Write(const string& outFileName)
{
    _cifFile.Write(outFileName);
}


void CifCorrector::CorrectEnumsSimple(CifFile& cifFile, DataInfo& dataInfo,
  const bool verbose)
{
    // Get list of all categories in the CIF file
    Block& block = cifFile.GetBlock(cifFile.GetFirstBlockName());

    vector<string> catNames;
    block.GetTableNames(catNames);

    for (unsigned int catI = 0; catI < catNames.size(); ++catI)
    {
        const string& catName = catNames[catI];

        // Get category table pointer.
        ISTable* catTableP = block.GetTablePtr(catName);

        // Get category attributes.
        const vector<string>& attrNames = catTableP->GetColumnNames();

        for (unsigned int attrI = 0; attrI < attrNames.size(); ++attrI)
        {
            const string& attribName = attrNames[attrI];

            // Make a CIF item.
            string item;
            CifString::MakeCifItem(item, catName, attribName);

            if (!dataInfo.IsItemDefined(item))
            {
#ifdef VLAD_DELETE
                cerr << "Warning: Item \"" << item << "\" is not defined in "\
                  "the internal dictionary." << endl;
#endif
            }

            // Check to see if column is present
            if (!catTableP->IsColumnPresent(attribName))
            {
#ifdef VLAD_DELETED
                cerr << "Warning: Table \"" << catName << "\" does not "\
                  "have column \"" << attribName << "\"." << endl;
#endif
                continue;
            }

            vector<string> enums = dataInfo.GetItemAttribute(item,
              CifString::CIF_DDL_CATEGORY_ITEM_ENUMERATION,
              CifString::CIF_DDL_ITEM_VALUE);

            // Check if the attribute is defined as an enumeration
            if (enums.empty())
            {
#ifdef VLAD_DELETE
                cerr << "Warning: Item \"" << item << "\" does not have "\
                  "enumerations defined." << endl;
#endif
                continue;
            }

            for (unsigned int rowI = 0; rowI < catTableP->GetNumRows(); ++rowI)
            {
                const string& currCellValue = (*catTableP)(rowI, attribName);

                // Find what enumeration it is to be replaced with by finding
                // enumInd within enumerations.

                unsigned int enumInd =
                  CifCorrector::FindEnumIndex(currCellValue, enums);

                if (enumInd == enums.size())
                {
                    // No enumerations fix is done
                    continue;
                }

                if (currCellValue == enums[enumInd])
                {
                    continue;
                }

                if (verbose)
                {
                    cout << "Info: Changing item \"" << item <<
                      "\" value from \"" << currCellValue << "\"";
                }

                catTableP->UpdateCell(rowI, attribName, enums[enumInd]);

                if (verbose)
                {
                    cout << " to \"" << enums[enumInd] << "\"" << endl;
                }
            } // for (all rows in category table)
        } // for (all attributes of a category)
    } // for (all categories in CIF file)
} // End of CifCorrector::CorrectEnumsSimple()


void CifCorrector::ValidateConfigTable()
{
    if (_configTableP == NULL)
    {
        return;
    }

    if (!_configTableP->IsColumnPresent("oper"))
    {
        throw NotFoundException("Table \"" + _configTableP->GetName() +
          "\" is missing attribute \"oper\".",
          "CifCorrector::ValidateConfigTable");
    }

    if (!_configTableP->IsColumnPresent("item"))
    {
        throw NotFoundException("Table \"" + _configTableP->GetName() +
          "\" is missing attribute \"item\".",
          "CifCorrector::ValidateConfigTable");
    }

    if (!_configTableP->IsColumnPresent("item_value"))
    {
        throw NotFoundException("Table \"" + _configTableP->GetName() +
          "\" is missing attribute \"item_value\".",
          "CifCorrector::ValidateConfigTable");
    }

    if (!_configTableP->IsColumnPresent("ref_item"))
    {
        throw NotFoundException("Table \"" + _configTableP->GetName() +
          "\" is missing attribute \"ref_item\".",
          "CifCorrector::ValidateConfigTable");
    }

    if (!_configTableP->IsColumnPresent("ref_item_value"))
    {
        throw NotFoundException("Table \"" + _configTableP->GetName() +
          "\" is missing attribute \"ref_item_value\".",
          "CifCorrector::ValidateConfigTable");
    }
}


unsigned int CifCorrector::FindEnumIndex(const string& value,
  const vector<string>& enums)
{
    for (unsigned int enumI = 0; enumI < enums.size(); ++enumI)
    {
        string modValue;
        String::RemoveWhiteSpace(value, modValue);
        String::LowerCase(modValue);

        string modEnum;
        String::RemoveWhiteSpace(enums[enumI], modEnum);
        String::LowerCase(modEnum);

        if (modValue == modEnum)
        {
            return (enumI);
        }
    }

    return (enums.size());
}


void CifCorrector::FixNumericList(string& outValue, const string& inValue)
{
    outValue = inValue;
    if (inValue.empty())
    {
        return;
    }

    // Replace semicolons with commas
    string::size_type semColInd = 0;
    while (semColInd != string::npos)
    {
        // VLAD - BUG semColInd + 1, but what to do with 0 index
        semColInd = outValue.find(';', semColInd);
        if (semColInd != string::npos)
        {
            outValue[semColInd] = ',';
        }
    }
}


void CifCorrector::FixNotApplicable(string& outValue, const string& inValue)
{
    if ((inValue == "N/A") || (inValue == "n/a") || (inValue == "N.A.") ||
      (inValue == "n.a."))
    {
        outValue = CifString::UnknownValue;
    }
    else
    {
        outValue = inValue;
    }
}


void CifCorrector::FixBadSequence(vector<string>& outValues,
  vector<pair<vector<string>, vector<string> > > & refMap,
  vector<vector<string> >& confRefValues, Block& block,
  const unsigned int numRows)
{
    outValues.clear();

    vector<string> prevRefValues;
    vector<string> currRefValues;

    for (unsigned int rowI = 0; rowI < numRows; ++rowI)
    {
        for (unsigned int refI = 1; refI < refMap.size(); ++refI)
        {
            const vector<string>& currSrcItems = refMap[refI].first;
            const vector<string>& currRefItems = refMap[refI].second;

            const vector<string>& prevRefItems = refMap[refI - 1].second;

            GetSrcValues(prevRefValues, currSrcItems, prevRefItems,
              confRefValues[refI], block, rowI);

            if (prevRefValues.empty())
            {
                return;
            }

            GetRefValues(currRefValues, currSrcItems, prevRefValues,
              currRefItems, block);

            if (currRefValues.empty())
            {
                return;
            }

            prevRefValues = currRefValues;

            if (refI == refMap.size() - 1)
            {
                outValues.push_back(currRefValues[0]);
            }
        }
    }
}


void CifCorrector::GetSrcValues(vector<string>& srcValues,
  const vector<string>& srcItems, const vector<string>& prevRefItems,
  const vector<string>& confValues, Block& block, const unsigned int rowIndex)
{
    srcValues.clear();

    if (srcItems.size() != confValues.size())
    {
        // VLAD - CORRECT THIS NOT TO THROW BUT TO INFORM USER
        throw InvalidStateException("Src items and conf values diff size.",
          "CifCorrector::GetSrcValues");
    }

    vector<string> tablePrevRefItemsValues;
    GetTableValues(tablePrevRefItemsValues, prevRefItems, block, rowIndex);

    if (tablePrevRefItemsValues.empty())
    {
        return;
    }

    for (unsigned int itemI = 0; itemI < srcItems.size(); ++itemI)
    {
        if (confValues[itemI] == CifString::InapplicableValue)
           srcValues.push_back(tablePrevRefItemsValues[itemI]);
        else
           srcValues.push_back(confValues[itemI]);
    }
}


void CifCorrector::GetRefValues(vector<string>& refValues,
  const vector<string>& srcItems, const vector<string>& srcValues,
  const vector<string>& refItems, Block& block)
{
    refValues.clear();

    // All source items must be in the same category. Here, the first
    // item is used to get the category name.
    string catName;
    CifString::GetCategoryFromCifItem(catName, srcItems[0]);

    ISTable* catTableP = block.GetTablePtr(catName);
    if (catTableP == NULL)
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Table \"" << catName << "\" not found." << endl;
#endif
        return;
    }

    // Prepare source attributes that will be searched
    vector<string> srcAttribs(srcItems.size());
    for (unsigned int itemI = 0; itemI < srcItems.size(); ++itemI)
    {
        CifString::GetItemFromCifItem(srcAttribs[itemI], srcItems[itemI]);
    }

    // Search for values
    unsigned int found = catTableP->FindFirst(srcValues, srcAttribs);
    if (found == catTableP->GetNumRows())
    {
        cerr << "ERROR: Table \"" << catName << "\" does not have these "\
          "values:" << endl;

        for (unsigned int itemI = 0; itemI < srcAttribs.size(); ++itemI)
        {
            cerr << "\"" << srcAttribs[itemI] << "\" = \"" <<
              srcValues[itemI] << "\"" << endl;
        }

        return;
    }

    vector<string> refAttribs(refItems.size());
    for (unsigned int refI = 0; refI < refItems.size(); ++refI)
    {
        CifString::GetItemFromCifItem(refAttribs[refI], refItems[refI]);
    }

    for (unsigned int refI = 0; refI < refAttribs.size(); ++refI)
    {
        refValues.push_back((*catTableP)(found, refAttribs[refI]));
    }
}


void CifCorrector::CreateRefMap(
  vector<pair<vector<string>, vector<string> > >& refMap,
  const string& refCondItem)
{
    refMap.clear();

    vector<string> refStrings;
    GetValuesFromSeparatedString(refStrings, refCondItem, "|");

    for (unsigned int refI = 0; refI < refStrings.size(); ++refI)
    {
        vector<string> srcDstString;
        GetValuesFromSeparatedString(srcDstString, refStrings[refI], ";");

        if (srcDstString.size() != 2)
        {
            // VLAD - CORRECT THIS NOT TO THROW BUT TO INFORM USER
            throw InvalidStateException("Wrong syntax in table.",
              "CifCorrector::CreateRefMap");
        }

        vector<string> srcString;
        GetValuesFromSeparatedString(srcString, srcDstString[0], ",");
        vector<string> dstString;
        GetValuesFromSeparatedString(dstString, srcDstString[1], ",");
 
        refMap.push_back(make_pair(srcString, dstString));
    }
}


void CifCorrector::ExtractRefValues(vector<vector<string> >& refValues,
  const string& refCondItemValue)
{
    refValues.clear();

    vector<string> refStrings;
    GetValuesFromSeparatedString(refStrings, refCondItemValue, ";");

    for (unsigned int refI = 0; refI < refStrings.size(); ++refI)
    {
        vector<string> valString;
        GetValuesFromSeparatedString(valString, refStrings[refI], ",");

        refValues.push_back(valString);
    }
}


void CifCorrector::GetTableValues(vector<string>& values,
  const vector<string>& items, Block& block, const unsigned int rowIndex)
{
    values.clear();

    // All items must be in the same category. Here, the first
    // item is used to get the category name.
    string catName;
    CifString::GetCategoryFromCifItem(catName, items[0]);

    ISTable* catTableP = block.GetTablePtr(catName);
    if (catTableP == NULL)
    {
#ifdef VLAD_DELETED
        cerr << "Warning: Table \"" << catName << "\" not found." << endl;
#endif
        return;
    }

    for (unsigned int itemI = 0; itemI < items.size(); ++itemI)
    {
        string attrib;
        CifString::GetItemFromCifItem(attrib, items[itemI]);
    
        if (!catTableP->IsColumnPresent(attrib))
        {
            cerr << "ERROR: Table \"" << catName << "\" does not "\
              "have column \"" << attrib << "\"." << endl;

            values.clear();

            return;
        }

        values.push_back((*catTableP)(rowIndex, attrib));
    }
}


void CifCorrector::GetValuesFromSeparatedString(vector<string>& values,
  const string& valueString, const string& sep)
{
    values.clear();

    if (valueString.empty())
    {
        return;
    }

    string::size_type sepInd = 0;
    string::size_type valueStartInd = 0;
    while (sepInd != string::npos)
    {
        sepInd = valueString.find(sep, sepInd + 1);
        if (sepInd != string::npos)
        {
            // VLAD-BUG !!! ".,." fails to be parsed
            //if (valueStartInd != (sepInd - 1))
            {
                values.push_back(valueString.substr(valueStartInd, sepInd -
                  valueStartInd));
                valueStartInd = sepInd + 1;
            }
        }
        else
        {
            values.push_back(valueString.substr(valueStartInd));
        }
    }
}


// VLAD - MOVE OR DELETE
void CifCorrector::NumberToLetter(char& letter, const unsigned int number)
{
    if (number > ('Z' - 'A'))
    {
        throw InvalidStateException("Invalid ASCII code offset \"" +
          String::IntToString(number) + "\"", "CifCorrector::NumberToLetter");
    }

    letter = 'A' + number;
}

