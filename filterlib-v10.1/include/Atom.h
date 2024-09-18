/*
FILE:     Atom.h
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
** \file Atom.h
**
** \brief Header file for Atom class.
*/

#ifndef _H_Atom_H_
#define _H_Atom_H_

#include <map>
#include <string>
#include <vector>

#include "Coord.h"

#define FORMAT_LONG  23
#define FORMAT_SIGATM 5
#define FORMAT_SIGUIJ 6

/**
**  \class Atom
** 
**  \brief Public class that respresents ATOM record.
** 
**  The class provides methods for construction, destruction, clearance, assignment
**  operator, data update and retrieval.
*/

namespace RCSB {

class Atom
{
   private:
        // atomic coordinates
        COORD  _orig;

        // ATOM record to store the following: 
        // "group_PDB",           index  0
        // "id",                  index  1
        // "label_atom_id",       index  2
        // "label_alt_id",        index  3
        // "label_comp_id",       index  4
        // "label_asym_id",       index  5
        // "label_seq_id",        index  6
        // "pdbx_PDB_ins_code",   index  7
        // "Cartn_x",             index  8
        // "Cartn_y",             index  9
        // "Cartn_z",             index 10
        // "occupancy",           index 11
        // "B_iso_or_equiv",      index 12
        // "type_symbol",         index 13
        // "pdbx_formal_charge",  index 14
        // "auth_seq_id",         index 15
        // "auth_comp_id",        index 16
        // "auth_asym_id",        index 17
        // "auth_atom_id",        index 18
        // "pdbx_auth_seq_id",    index 19
        // "pdbx_auth_comp_id",   index 20
        // "pdbx_auth_asym_id",   index 21
        // "pdbx_auth_atom_name", index 22
        // "label_entity_id"      index 23
        std::vector<std::string> _atom;

        // Same size as _atom, true value means the original value of _atom is empty.
        std::vector<bool> _empty_flag;

        // SIGATM record to store the following:
        // "Cartn_x_esd",         index  0
        // "Cartn_y_esd",         index  1
        // "Cartn_z_esd",         index  2
        // "occupancy_esd",       index  3
        // "B_iso_or_equiv_esd",  index  4
        std::vector<std::string> _sigatm;

        // ANISOU record to store the following:
        // "U[1][1]",             index  0
        // "U[2][2]",             index  1
        // "U[3][3]",             index  2
        // "U[1][2]",             index  3
        // "U[1][3]",             index  4
        // "U[2][3]",             index  5
        std::vector<std::string> _anisou;

        // SIGUIJ record to store the following:
        // "U[1][1]_esd",         index  0
        // "U[2][2]_esd",         index  1
        // "U[3][3]_esd",         index  2
        // "U[1][2]_esd",         index  3
        // "U[1][3]_esd",         index  4
        // "U[2][3]_esd",         index  5
        std::vector<std::string> _siguij;

        // Flag to indicate if the Atom is deleted
        // bool _deleted;

        // Flag to indicate if there is a TER card after this Atom
        // 0: no TER card after this atom
        // number: line no in PDB format file for TER record
        int _ter_flag;

        // line number in PDB format file for this atom or row number in atom_site category
        int _lineno;

        // Flag to indicate if the coordinate is 9999.000 9999.000 9999.000
        bool _is_9999_flag;

        // empty string data
        std::string _empty;

        // cif item vs. record position mapping for _atom
        std::map<std::string, short int> _atom_item_position_mapping;

        // atom_site_anisotrop cif item vs. record position mapping for _atom
        std::map<std::string, short int> _anisotrop_item_position_mapping;

        // cif item vs. record position mapping for _sigatm
        std::map<std::string, short int> _sigatm_item_position_mapping;

        // cif item vs. record position mapping for _anisou
        std::map<std::string, short int> _anisou_item_position_mapping;

        // cif item vs. record position mapping for _siguij
        std::map<std::string, short int> _siguij_item_position_mapping;

        std::map<std::string, std::string> _extra_item_value_mapping;

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
        void _clear_data();
   public:
        /**
        **  Constructs Atom
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
        Atom() { clear(); }

        /**
        **  Destructs Atom
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
        ~Atom() { _clear_data(); }

        /**
        **  Clear data storage, set up initial data storage
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
        void  clear();

        /**
        **  Copy a Atom to another Atom (assignment operator).
        **
        **  \param[in]: atom - reference to the source Atom
        **
        **  \return Reference to the destination Atom
        **
        **  \pre None
        **
        **  \post Constructed Atom is a clone as the Atom
        **        referenced by \e atom
        **
        **  \exception: None
        */
        Atom& operator=(const Atom& atom);

        /**
        **  Retrieves ATOM type
        **
        **  \param: None
        **
        **  \return Constant reference to a string that contains ATOM type
        **
        **  \pre None
        **
        **  \post None
        **
        **  \exception: None
        */
        const std::string& type() const;
        const std::string& atnum() const;
        const std::string& atmtype() const;
        const std::string& alt_loc() const;
        const std::string  alt_loc_char() const;
        const std::string& restype() const;
        const std::string& chnid() const;
        const std::string& resnum() const;
        const std::string& ins_code() const;
        const std::string  ins_code_char() const;
        const std::string& occ() const;
        const std::string& t_fct() const;
        const std::string& atom_type() const;
        const std::string& charge() const;
        const std::string& pdb_resnum() const;
        const std::string& pdb_resnam() const;
        const std::string& pdb_chnid() const;
        const std::string  pdb_chnid_char() const;
        const std::string& pdb_atmnam() const;
        const std::string& pub_resnum() const;
        const std::string& pub_resnam() const;
        const std::string& pub_chnid() const;
        const std::string& pub_atmnam() const;
        const COORD& orig() const;
        const double& x() const;
        const double& y() const;
        const double& z() const;
        const bool is_hydrogen() const;
        const bool has_sigatm() const; 
        const bool has_anisou() const; 
        const bool has_siguij() const; 
        // const bool& deleted() const;
        const int& ter_flag() const;
        const int& lineno() const;
        const bool& is_9999_flag() const;
        const std::string getAtomIndex() const;
        const std::string getAtomAltIndex() const;
        const std::string getAtomAllIndex(const bool& alt_flag = true) const;
        const std::string getAtomAllIndexWithoutAltId() const;
        const std::string& getValue(const int& pos) const;
        const std::string& getOrigValue(const int& pos) const;
        const std::string& getValue(const std::string&);
        const std::string& getAnisouValue(const std::string&);
        const std::vector<std::string>& getValue() const;
        void  getAuxiliaryValue(const std::string&, std::vector<std::string>&);
        void  getAtomIndex(std::vector<std::string>&, const bool& clear_flag = false, const bool& exclude_alt_loc_flag = false);
        void  getResidueIndex(std::vector<std::string>&, const bool& clear_flag = false);

        void set_type(const std::string&);
        void set_atnum(const int&);
        void set_atmtype(const std::string&);
        void set_alt_loc(const std::string&);
        void set_restype(const std::string&);
        void set_chnid(const std::string&);
        void set_resnum(const std::string&);
        void set_ins_code(const std::string&);
        void set_atom_type(const std::string&);
        void set_charge(const std::string&);
        void set_pdb_resnum(const std::string&);
        void set_pdb_resnam(const std::string&);
        void set_pdb_chnid(const std::string&);
        void set_pdb_atmnam(const std::string&);
        // void set_deletion();
        void set_ter_flag(const int&);
        void set_orig(const COORD&);
        void set_x(const double&);
        void set_y(const double&);
        void set_z(const double&);
        void setValue(const std::string&, const int&);
        void setValue(const std::vector<std::string>&, const int&);
        void setAuxiliaryValue(const std::vector<std::string>&, const bool convert_flag = false);
        void setExtraValue(const std::string& item, const std::string& value);
        void UpperCase();
        void writeAtomRecord(std::string& record);
};

}

extern double cal_distance(const COORD&, const COORD&);
extern double cal_distance(const RCSB::Atom*, const RCSB::Atom*);
extern bool is_same_atom(const RCSB::Atom*, const RCSB::Atom*);

extern const char *_format_atom_site_pdbx[FORMAT_LONG + 1];
extern const char *_format_sigatm_pdbx[FORMAT_SIGATM];
extern const char *_format_anisotrop_pdbx[FORMAT_LONG];
extern const char *_format_siguij_pdbx[FORMAT_SIGUIJ];

#endif

