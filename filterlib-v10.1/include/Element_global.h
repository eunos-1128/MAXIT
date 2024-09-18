/*
FILE:     Element_global.h
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
#ifndef _H_ELEMENT_GLOBAL_H_
#define _H_ELEMENT_GLOBAL_H_

#include <map>
#include <string>

#define NUM_PERIODIC_TABLE_ATOM 110

typedef struct {
        const char *Atom_Symbol;
        const char *Atom_Name;
        double Atom_Weight;
        double Covalent_Radii;
        double Ionic_Radii;
        int is_metal_element;
} ELEMENT_PERIODIC_TABLE;

static ELEMENT_PERIODIC_TABLE periodic_table[NUM_PERIODIC_TABLE_ATOM] = {
       { "H",   "Hydrogen",       1.0079,    0.32, 0.79, 0 },
       { "He",  "Helium",         4.0026,    0.93, 0.49, 0 },
       { "Li",  "Lithium",        6.941,     1.23, 2.05, 1 },
       { "Be",  "Beryllium",      9.01218,   0.90, 1.40, 1 },
       { "B",   "Boron",          10.81,     0.82, 1.17, 0 },
       { "C",   "Carbon",         12.011,    0.77, 0.91, 0 },
       { "N",   "Nitrogen",       14.0067,   0.75, 0.75, 0 },
       { "O",   "Oxygen",         15.9994,   0.73, 0.65, 0 },
       { "F",   "Fluorine",       18.998403, 0.72, 0.57, 0 },
       { "Ne",  "Neon",           20.179,    0.71, 0.51, 0 },
       { "Na",  "Sodium",         22.98977,  1.54, 2.23, 1 },
       { "Mg",  "Magnesium",      24.305,    1.36, 1.72, 1 },
       { "Al",  "Aluminum",       26.98154,  1.18, 1.82, 1 },
       { "Si",  "Silicon",        28.0855,   1.11, 1.46, 0 },
       { "P",   "Phosphorus",     30.97376,  1.06, 1.23, 0 },
       { "S",   "Sulfur",         32.06,     1.02, 1.09, 0 },
       { "Cl",  "Chlorine",       35.453,    0.99, 0.97, 0 },
       { "Ar",  "Argon",          39.948,    0.98, 0.88, 0 },
       { "K",   "Potassium",      39.0983,   2.03, 2.77, 1 },
       { "Ca",  "Calcium",        40.08,     1.91, 2.23, 1 },
       { "Sc",  "Scandium",       44.9559,   1.62, 2.09, 1 },
       { "Ti",  "Titanium",       47.9,      1.45, 2.00, 1 },
       { "V",   "Vanadium",       50.9415,   1.34, 1.92, 1 },
       { "Cr",  "Chromium",       51.996,    1.18, 1.85, 1 },
       { "Mn",  "Manganese",      54.938,    1.17, 1.79, 1 },
       { "Fe",  "Iron",           55.847,    1.17, 1.72, 1 },
       { "Co",  "Cobalt",         58.9332,   1.16, 1.67, 1 },
       { "Ni",  "Nickel",         58.7,      1.15, 1.62, 1 },
       { "Cu",  "Copper",         63.546,    1.17, 1.57, 1 },
       { "Zn",  "Zinc",           65.38,     1.25, 1.53, 1 },
       { "Ga",  "Gallium",        69.72,     1.26, 1.81, 1 },
       { "Ge",  "Germanium",      72.59,     1.22, 1.52, 1 },
       { "As",  "Arsenic",        74.9216,   1.20, 1.33, 1 },
       { "Se",  "Selenium",       78.96,     1.36, 1.22, 0 },
       { "Br",  "Bromine",        79.904,    1.14, 1.12, 0 },
       { "Kr",  "Krypton",        83.8,      1.12, 1.03, 0 },
       { "Rb",  "Rubidium",       85.4678,   2.16, 2.98, 1 },
       { "Sr",  "Strontium",      87.62,     1.91, 2.45, 1 },
       { "Y",   "Yttrium",        88.9059,   1.62, 2.27, 1 },
       { "Zr",  "Zirconium",      91.22,     1.45, 2.16, 1 },
       { "Nb",  "Niobium",        92.9064,   1.34, 2.09, 1 },
       { "Mo",  "Molybdenum",     95.94,     1.30, 2.01, 1 },
       { "Tc",  "Technetium",     98.0,      1.27, 1.95, 1 },
       { "Ru",  "Ruthenium",      101.07,    1.25, 1.89, 1 },
       { "Rh",  "Rhodium",        102.9055,  1.25, 1.83, 1 },
       { "Pd",  "Palladium",      106.4,     1.28, 1.79, 1 },
       { "Ag",  "Silver",         107.868,   1.34, 1.75, 1 },
       { "Cd",  "Cadmium",        112.41,    1.48, 1.71, 1 },
       { "In",  "Indium",         114.82,    1.44, 2.00, 1 },
       { "Sn",  "Tin",            118.69,    1.41, 1.72, 1 },
       { "Sb",  "Antimony",       121.75,    1.40, 1.53, 1 },
       { "Te",  "Tellurium",      127.6,     1.36, 1.42, 1 },
       { "I",   "Iodine",         126.9045,  1.33, 1.32, 0 },
       { "Xe",  "Xenon",          131.3,     1.31, 1.24, 0 },
       { "Cs",  "Cesium",         132.9054,  2.00, 3.00, 1 },
       { "Ba",  "Barium",         137.33,    1.98, 2.78, 1 },
       { "La",  "Lanthanum",      138.9055,  1.69, 2.74, 1 },
       { "Ce",  "Cerium",         140.12,    1.44, 2.16, 1 },
       { "Pr",  "Praseodymium",   140.9077,  1.34, 2.09, 1 },
       { "Nd",  "Neodymium",      144.24,    1.30, 2.02, 1 },
       { "Pm",  "Promethium",     145.0,     1.28, 1.97, 1 },
       { "Sm",  "Samarium",       150.4,     1.26, 1.92, 1 },
       { "Eu",  "Europium",       151.96,    1.27, 1.87, 1 },
       { "Gd",  "Gadolinium",     157.25,    1.30, 1.83, 1 },
       { "Tb",  "Terbium",        158.9254,  1.34, 1.79, 1 },
       { "Dy",  "Dysprosium",     162.5,     1.49, 1.76, 1 },
       { "Ho",  "Holmium",        164.9304,  1.48, 2.08, 1 },
       { "Er",  "Erbium",         167.26,    1.47, 1.81, 1 },
       { "Tm",  "Thulium",        168.9342,  1.46, 1.63, 1 },
       { "Yb",  "Ytterbium",      173.04,    1.46, 1.53, 1 },
       { "Lu",  "Lutetium",       174.967,   1.45, 1.43, 1 },
       { "Hf",  "Hafnium",        178.49,    1.43, 1.34, 1 },
       { "Ta",  "Tantalum",       180.9479,  1.50, 3.50, 1 },
       { "W",   "Tungsten",       183.85,    2.40, 3.00, 1 },
       { "Re",  "Rhenium",        186.207,   2.20, 3.20, 1 },
       { "Os",  "Osmium",         190.2,     1.65, 2.70, 1 },
       { "Ir",  "Iridium",        192.22,    1.65, 2.67, 1 },
       { "Pt",  "Platinum",       195.09,    1.64, 2.64, 1 },
       { "Au",  "Gold",           196.9665,  1.63, 2.62, 1 },
       { "Hg",  "Mercury",        200.59,    1.62, 2.59, 1 },
       { "Tl",  "Thallium",       204.37,    1.85, 2.56, 1 },
       { "Pb",  "Lead",           207.2,     1.61, 2.54, 1 },
       { "Bi",  "Bismuth",        208.9804,  1.59, 2.51, 1 },
       { "Po",  "Polonium",       209.0,     1.59, 2.49, 1 },
       { "At",  "Astatine",       210.0,     1.58, 2.47, 0 },
       { "Rn",  "Radon",          222.0,     1.57, 2.45, 0 },
       { "Fr",  "Francium",       223.0197,  1.56, 2.42, 1 },
       { "Ra",  "Radium",         226.0254,  1.74, 2.40, 1 },
       { "Ac",  "Actinium",       227.0278,  1.56, 2.25, 1 },
       { "Th",  "Thorium",        232.0381,  1.65, 3.16, 1 },
       { "Pa",  "Protactinium",   231.0359,  1.65, 3.14, 1 },
       { "U",   "Uranium",        238.029,   1.42, 3.11, 1 },
       { "Np",  "Neptunium",      237.0482,  1.65, 3.08, 1 },
       { "Pu",  "Plutonium",      244.0,     1.65, 3.05, 1 },
       { "Am",  "Americium",      243.0,     1.65, 3.02, 1 },
       { "Cm",  "Curium",         247.0,     1.65, 2.99, 1 },
       { "Bk",  "Berkelium",      247.0,     1.65, 2.97, 1 },
       { "Cf",  "Californium",    251.0,     1.65, 2.95, 1 },
       { "Es",  "Einsteinium",    252.0,     1.65, 2.92, 1 },
       { "Fm",  "Fermium",        257.0,     1.65, 2.90, 1 },
       { "Md",  "Mendelevium",    258.0,     1.65, 2.87, 1 },
       { "No",  "Nobelium",       259.0,     1.65, 2.85, 1 },
       { "Lr",  "Lawrencium",     260.0,     1.65, 2.82, 1 },
       { "Rf",  "Rutherfordium",  261.0,     0.32, 0.80, 1 },
       { "Db",  "Dubnuim",        262.0,     0.10, 0.10, 1 },
       { "Sg",  "Seaborgium",     263.0,     0.00, 0.00, 1 },
       { "Bh",  "Bohrium",        262.0,     0.00, 0.00, 1 },
       { "??",  "Hassium",        265.0,     0.00, 0.00, 1 },
       { "Mt",  "Meitnerium",     266.0,     0.00, 0.00, 1 },
       { "??",  "Element112",     277.0,     0.00, 0.00, 1 }
};

static std::map<std::string, int> _element_mapping;

#endif
