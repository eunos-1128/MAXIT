/*
FILE:     BasePairInfo.h
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
#ifndef _H_BASE_PAIR_INFO_H_
#define _H_BASE_PAIR_INFO_H_

#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "BasePairParamDef.h"
#include "BaseParameter.h"
#include "BPIndex.h"
#include "ConnectDic.h"
#include "Contact.h"
#include "Link.h"
#include "Molecule.h"

#define DOUBLE_HELIX           1
#define PARALLEL_HELIX         2
#define A_DOUBLE_HELIX         4
#define B_DOUBLE_HELIX         8
#define Z_DOUBLE_HELIX        16
#define TRIPLE_HELIX          32
#define QUADRUPLE_HELIX       64
#define MISMATCH             128
#define HAIRPIN_LOOP         256
#define TETRA_LOOP           512
#define INTERNAL_LOOP       1024
#define BULGE_LOOP          2048
#define THREE_STEM          4096
#define FOUR_STEM           8192

class BasePairInfo
{
   public:
        BasePairInfo();
        ~BasePairInfo();
        void setCell(CrySymmetry* cell);
        void setCCDic(ConnectDic* ccdic);
        void setMolecule(RCSB::Molecule* mol);
        void calculateBasePairInfo();
        void getBasePairs(std::list<_BSPAIR>& bspairs);
        void getClassification(std::vector<std::string>& classifications);
        const std::list<_BASEBASE_PARAMS>& getBaseBaseParams() const;
        const std::list<_INTERBASE_PARAMS>& getInterBaseParams() const;
   private:
        typedef struct {
               int i, j;
               double dsum, dv;
               COORD dir, oave, zave;
               int type;
               int type_i;
               int type_j;
               int cis_or_trans;
               std::vector<CONTACT> bp_contact;
        } _BASEPAIR;

        CrySymmetry* _cell;
        ConnectDic* _ccDic;
        RCSB::Molecule* _mol;
        int _classification;
        std::map<std::string, std::set<int> > _paired_chain_id_map;
        std::multimap<BPIndex, CONTACT> _basepair_contacts;
        std::list<_BASEBASE_PARAMS> _basebase_params;
        std::list<_INTERBASE_PARAMS> _interbase_params;
        std::vector<_BASE> _baseframes;
        std::vector<std::vector<_BASEPAIR> > _basepairs;
        std::vector<_BASEPAIR> _best_basepairs;

        void _clear();
        int  _getBaseFrames(const std::vector<CONTACT>& contact, const std::string& selected_chainid, const int& op, const int& lx,
                            const int& ly, const int& lz);
        int  _getBasePairInformation(const std::vector<CONTACT>& contact, const std::string& selected_chainid, const bool& count_paired_flag);
        int  _refineBasePairs(const bool& count_paired_flag);
        void _getBasePairs(std::vector<_BASEPAIR>& base_pairs);
        void _updatePairedChainIdMap(const int& idx);
        void _getBestBasePair();
        void _getSymOperator(const std::vector<CONTACT>& contact, std::vector<std::pair<std::string, std::vector<int> > >& chnid_symop_list);
        bool _check_base_pair(const _BASE& base_i, const _BASE& base_j, _BASEPAIR& bs_pair);
        bool _is_overlap(const _BASE& base_i, const _BASE& base_j);
        double _get_overlap_area(const _BASE& base_i, const _BASE& base_j, const COORD& oave, const COORD& zave);
        std::vector<COORD> _get_polygon(const std::vector<COORD>& polygon_in, const COORD& oave, const COORD rot_mtx[3]);
        int _getBestPair(const std::vector<_BASEPAIR>&  base_pair, _BASEPAIR& bp_pair, const std::vector<int>& skip);
        void _re_ordering();
        int  _get_bp_context(std::vector<std::vector<int> >& bp_order, std::vector<std::vector<int> >& end_list);
        std::vector<std::vector<int> > _locate_helix(const int& num_ends, const std::vector<std::vector<int> >& bp_order,
                                               const std::vector<std::vector<int> >& end_list, std::vector<int>& bp_idx);
        void _five2three(const std::vector<std::vector<int> >& helix_idx, std::vector<int>& bp_idx);
        void _get_idx_and_min(const int& i, const int& start, const int& end, const std::vector<double>& dist_matrix,
                              std::vector<int>& ddidx, std::vector<double>& ddmin);
        bool _is_circle_helix(std::vector<std::vector<int> >& bp_order);
        bool _distance_ab(const _BASE& a, const std::string& ia, const _BASE& b, const std::string& ib, double& dist);
        void _first_step(const std::vector<int>& helix_idx, std::vector<int>& swapped, std::vector<int>& bp_idx);
        void _get_ij(const int& m, const int& swapped, int &i, int &j);
        bool _wc_bporien(const int& fst, const int& snd, const int& i1, const int& j1, const int& i2, const int& j2);
        bool _check_o3dist(const int& i1, const int& j1, const int& i2, const int& j2);
        void _reverse(const int& st, const int& n, std::vector<int>& vec);
        int _is_linked(const _BASE& base1, const _BASE& base2);
        bool _is_wc_geometry(const int& bs_idx);
        int _check_others(const int& i1, const int& j1, const int& i2, const int& j2);
        bool _check_direction(const std::vector<int>& helix_idx, std::vector<int>& bp_idx, std::vector<int>& swapped, int *direction);
        void _get_parameters(const std::vector<std::vector<int> >& pair_num);
        void _bpstep_par(COORD rot1[3], const COORD& org1, COORD rot2[3], const COORD& org2, double *pars, COORD mst_orien[3], COORD &mst_org);
        void _helical_par(COORD rot1[3], const COORD& org1, COORD rot2[3], const COORD& org2, double *pars, COORD mst_orien[3], COORD &mst_org);
        void _project_xyzP(const COORD& P_1_i_plus_1, const COORD& P_2_i, COORD mst_orien[3], const COORD& mst_org, COORD& aveP);
        void _find_additional_classification();
};

#endif
