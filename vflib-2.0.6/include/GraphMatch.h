/*
FILE:     GraphMatch.h
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
** \file GraphMatch.h
**
** \brief Header file for GraphMatch class.
*/

#ifndef _H_GRAPH_MATCH_H_
#define _H_GRAPH_MATCH_H_

#include <map>
#include <set>
#include <string>
#include <vector>

#include "argraph.h"

#define  REF_TO_TARGET   1
#define  TARGET_TO_REF   2

/**
**  \class _Node_
** 
**  \brief Public class that respresents a node in graph.
** 
*/

class _Node_ {
     public:
       std::string type, name;
     public:
       _Node_(const std::string& node_type,
              const std::string& node_name) {
             type = node_type;
             name = node_name;
       }
};

/**
**  \class _Edge_
** 
**  \brief Public class that respresents an edge in graph.
** 
*/

class _Edge_ {
     public:
       int type;
     public:
       _Edge_(const int edge_type) {
             type = edge_type;
       }
};

/**
**  \class _Node_Destroyer
** 
**  \brief Public class that delete a _Node_
** 
*/

class _Node_Destroyer: public AttrDestroyer {
     public:
       virtual void destroy(void *p) {
             delete ((_Node_ *) p);
       }
};

/**
**  \class _Edge_Destroyer
** 
**  \brief Public class that delete a _Edge_
** 
*/

class _Edge_Destroyer: public AttrDestroyer {
     public:
       virtual void destroy(void *p) {
             delete ((_Edge_ *) p);
       }
};

/**
**  \class _Node_Comparator
** 
**  \brief Public class that compares two _Node_s.
** 
*/

class _Node_Comparator: public AttrComparator {
     public:
       _Node_Comparator() {}
       virtual bool compatible(void *ptA, void *ptB) {
             _Node_ *a = (_Node_ *) ptA;
             _Node_ *b = (_Node_ *) ptB;
             return (a->type == b->type);
       }
};

/**
**  \class _Edge_Comparator
** 
**  \brief Public class that compares two _Edge_s.
** 
*/

class _Edge_Comparator: public AttrComparator {
     public:
       _Edge_Comparator() {}
       virtual bool compatible(void *ptA, void *ptB) {
             return (true);
       }
};

/**
**  \class GraphMatch
** 
**  \brief Public class that does graph matching using VFLib library.
** 
*/

class GraphMatch
{
  public:
       /* 
       **  Run substructure graph matching using VFLib library with all atoms defined
       **  in ref_atomlist & tgt_atomlist
       **
       **  \param[in]: ref_atomlist - reference atom type/name pair list
       **  \param[in]: ref_bondlist - reference bond list
       **  \param[in]: tgt_atomlist - target atom type/name pair list
       **  \param[in]: tgt_bondlist - target bond list
       **  \param[out]: matchlist   - matched atom name mapping between reference and target
       **  \param[in]: order        - value to define key/value order in matchlist:
       **                             REF_TO_TARGET: key is reference atom, value is target atom 
       **                             TARGET_TO_REF: key is target atom, value is reference atom
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       static void VF_GetMatch(const std::vector<std::pair<std::string, std::string> >& ref_atomlist,
                               const std::vector<std::pair<std::string, std::string> >& ref_bondlist,
                               const std::vector<std::pair<std::string, std::string> >& tgt_atomlist,
                               const std::vector<std::pair<std::string, std::string> >& tgt_bondlist,
                               std::map<std::string, std::string>& matchlist,
                               const int& order = TARGET_TO_REF);

       /* 
       **  Two steps processing:
       **  1. Run substructure graph matching using VFLib library with heavy atoms only defined
       **     in ref_atomlist & tgt_atomlist
       **  2. Add hydrogens mapping for found match(es)
       **  
       **  The method is same as VF_GetMatch if target has no hydrogen. It has much better 
       **  performance if hydrogens are presented.
       **
       **  \param[in]: ref_atomlist - reference atom type/name pair list
       **  \param[in]: ref_bondlist - reference bond list
       **  \param[in]: tgt_atomlist - target atom type/name pair list
       **  \param[in]: tgt_bondlist - target bond list
       **  \param[out]: matchlist   - matched atom name mapping between reference and target
       **  \param[in]: order        - value to define key/value order in matchlist:
       **                             REF_TO_TARGET: key is reference atom, value is target atom 
       **                             TARGET_TO_REF: key is target atom, value is reference atom
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       static void GetMatch(const std::vector<std::pair<std::string, std::string> >& ref_atomlist,
                            const std::vector<std::pair<std::string, std::string> >& ref_bondlist,
                            const std::vector<std::pair<std::string, std::string> >& tgt_atomlist,
                            const std::vector<std::pair<std::string, std::string> >& tgt_bondlist,
                            std::map<std::string, std::string>& matchlist,
                            bool& is_substructure_match, const int& order = TARGET_TO_REF,
                            const bool& heavy_atom_only = false, const bool& exclude_OP3_flag = false);

       /* 
       **  Add hydrogens mapping to existing heavy atom mapping
       **  
       **  \param[out]: atom_mapping_list - key/value atom mapping list
       **  \param[in]: linked_H_atoms1 - linked hydrogens for key atoms in atom_mapping_list
       **  \param[in]: linked_H_atoms2 - linked hydrogens for value atoms in atom_mapping_list
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       static void AddHydrogenMapping(std::map<std::string, std::string>& atom_mapping_list,
                             const std::map<std::string, std::set<std::string> >& linked_H_atoms1,
                             const std::map<std::string, std::set<std::string> >& linked_H_atoms2);
  private:
       /* 
       **  Build VFLib library's ARGraph
       **
       **  \param[in]: atomlist - atom type/name pair list
       **  \param[in]: bondlist - bond list
       **
       **  \return pointer to ARGraph
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       static ARGraph<_Node_*, _Edge_*> *getGraph(const std::vector<std::pair<std::string,
                                                  std::string> >& atomlist,
                                                  const std::vector<std::pair<std::string,
                                                  std::string> >& bondlist);

       /* 
       **  Run substructure graph matching using VFLib library with all atoms defined
       **  in ref_atomlist & tgt_atomlist
       **
       **  \param[in]: ref_atomlist - reference atom type/name pair list
       **  \param[in]: ref_bondlist - reference bond list
       **  \param[in]: tgt_atomlist - target atom type/name pair list
       **  \param[in]: tgt_bondlist - target bond list
       **  \param[out]: match_pair_lists - all possible matched atom name mappings between
       **                                  reference and target
       **  \param[in]: order        - value to define key/value order in matchlist:
       **                             REF_TO_TARGET: key is reference atom, value is target atom 
       **                             TARGET_TO_REF: key is target atom, value is reference atom
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       static void DoMatch(const std::vector<std::pair<std::string, std::string> >& ref_atomlist,
                           const std::vector<std::pair<std::string, std::string> >& ref_bondlist,
                           const std::vector<std::pair<std::string, std::string> >& tgt_atomlist,
                           const std::vector<std::pair<std::string, std::string> >& tgt_bondlist,
                           std::vector<std::vector<std::pair<std::string, std::string> > >&
                           match_pair_lists, const int& order);

       /* 
       **  Find best matched atom name mapping considering PDB atom name convention and connected
       **  hydrogens
       **
       **  \param[out]: matchlist   - best matched atom name mapping
       **  \param[in]: linked_H_atoms1 - linked hydrogens for pair.first atoms
       **  \param[in]: linked_H_atoms2 - linked hydrogens for pair.second atoms
       **  \param[in]: match_pair_lists - all possible matched atom name mappings between
       **                                  reference and target
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       static void GetBestMatch(std::map<std::string, std::string>& matchlist,
                           const std::map<std::string, std::set<std::string> >& linked_H_atoms1,
                           const std::map<std::string, std::set<std::string> >& linked_H_atoms2,
                           const std::vector<std::vector<std::pair<std::string, std::string> > >& match_pair_lists,
                           const bool& exclude_OP3_flag = false);

       /* 
       **  Find best matched atom name mapping index in matchlists
       **  hydrogens
       **
       **  \param[in]: linked_H_atoms1 - linked hydrogens for pair.first atoms
       **  \param[in]: linked_H_atoms2 - linked hydrogens for pair.second atoms
       **  \param[in]: linked_H_atoms1 - linked hydrogens for key atoms in matchlists
       **                                reference and target
       **
       **  \return index to matchlists vector
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       static int  GetBestMatchIndex(const std::map<std::string, std::set<std::string> >& linked_H_atoms1,
                                     const std::map<std::string, std::set<std::string> >& linked_H_atoms2,
                                     const std::vector<std::vector<std::pair<std::string, std::string> > >&
                                     matchlists, const bool& exclude_OP3_flag = false);

       /* 
       **  Reorder atom name character set alphabetically
       **
       **  \param[out]: name - atom name
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       static void reorderName(std::string& name);

       /* 
       **  Check if two atom names are close.
       **
       **  \param[in]: name1 - atom name
       **  \param[in]: name2 - atom name
       **
       **  \return true - if two names are close
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       static bool isCloseName(const std::string& name1, const std::string& name2);

       /* 
       **  Get two atom names similarity score
       **
       **  \param[in]: name1 - atom name
       **  \param[in]: name2 - atom name
       **
       **  \return 2  -  atom names are same
       **          1  -  atom names are close
       **          0  -  other
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       static int getAtomNameScore(const std::string& name1, const std::string& name2);

       /* 
       **  Generate heavy atoms only atom type/name pair list & bond list
       **
       **  \param[in]:  alist      - atom type/name pair list
       **  \param[out]: alist_no_H - atom type/name pair list without hydrogens
       **  \param[in]:  blist      - bond list
       **  \param[out]: blist_no_H - bond list without hydrogens
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       static void removeHydrogenAtoms(const std::vector<std::pair<std::string, std::string> >& alist,
                                       std::vector<std::pair<std::string, std::string> >& alist_no_H,
                                       const std::vector<std::pair<std::string, std::string> >& blist,
                                       std::vector<std::pair<std::string, std::string> >& blist_no_H);

       /* 
       **  Generate heavy atom vs linked hydrogen(s) mapping
       **
       **  \param[out]: atoms - heavy atom vs linked hydrogen(s) mapping
       **  \param[in]:  alist - atom type/name pair list
       **  \param[in]:  blist - bond list
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       static void getLinkedHydrogenAtoms(std::map<std::string, std::set<std::string> >& atoms,
                               const std::vector<std::pair<std::string, std::string> >& a_list,
                               const std::vector<std::pair<std::string, std::string> >& b_list);

       /* 
       **  Generate hydrogen pair list
       **
       **  \param[out]: atom_pair - hydrogen pair list
       **  \param[in]:  set1      - hydrogen atom set
       **  \param[in]:  set2      - hydrogen atom set
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       static void getHydrogenAtomPair(std::vector<std::pair<std::string, std::string> >& atom_pair,
                               const std::set<std::string>& set1, const std::set<std::string>& set2);

       /* 
       **  Convert hydrogen name set into list based on numbering
       **
       **  \param[in]:  atom_set  - hydrogen atom set
       **  \param[out]: atom_list - hydrogen atom list
       **
       **  \return true if all hydrogen names contain number
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       static bool getHydrogenList(const std::set<std::string>& atom_set, std::vector<std::string>& atom_list);
};

#endif
