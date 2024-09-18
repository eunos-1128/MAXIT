/*
FILE:     MatthewsCoefficientSolventContentDistribution.h
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

#ifndef _H_MATTHEWS_COEFFICIENT_SOLVENT_CONTENT_DISTRIB_VALUE_H_
#define _H_MATTHEWS_COEFFICIENT_SOLVENT_CONTENT_DISTRIB_VALUE_H_

#define NUM_DISTRIB 100
#define NUM_VALUES  8

static const double __over_all_values[4] = { 1.410000, 6.890000, 14.9934, 82.07220 };

static const double __distrib_values[NUM_DISTRIB][NUM_VALUES] = {
       {      0.00,  1.88e+04,  0.00,    1.5,  1.410000,    3.5122,   14.9934,   64.7696  },
       {      0.00,  1.88e+04,   1.5,    1.7,  1.592700,    4.0792,   23.0567,   69.4569  },
       {      0.00,  1.88e+04,   1.7,   1.84,    1.5800,    4.1551,   25.1490,   70.3600  },
       {      0.00,  1.88e+04,  1.84,   1.95,    1.6200,    4.3650,  24.30000,   72.5100  },
       {      0.00,  1.88e+04,  1.95,   2.05,  1.683100,    4.3638,   26.6585,   71.7252  },
       {      0.00,  1.88e+04,  2.05,    2.2,    1.7381,    4.4676,   28.6804,   72.3919  },
       {      0.00,  1.88e+04,   2.2,   2.39,    1.7640,   4.95632,   30.4896,   75.1784  },
       {      0.00,  1.88e+04,  2.39,    2.6,    1.7417,    5.0166,   29.7238,   74.9796  },
       {      0.00,  1.88e+04,   2.6,   2.85,    1.6928,   5.70560,   33.5208,   78.3800  },
       {      0.00,  1.88e+04,  2.85,  100.0,    1.9030,  7.306973,   35.1970,  83.16597  },
       {  1.88e+04,  2.74e+04,  0.00,    1.5,  1.646200,    3.3690,   19.4878,   62.9738  },
       {  1.88e+04,  2.74e+04,   1.5,    1.7,  1.707412,    3.7900,   24.9000,   67.5325  },
       {  1.88e+04,  2.74e+04,   1.7,   1.84,    1.7593,    3.9107,   29.4100,   68.5363  },
       {  1.88e+04,  2.74e+04,  1.84,   1.95,    1.7600,    4.3800,  29.99100,   71.9090  },
       {  1.88e+04,  2.74e+04,  1.95,   2.05,  1.731655,    4.4095,   29.7786,   72.3814  },
       {  1.88e+04,  2.74e+04,  2.05,    2.2,    1.8123,    4.9116,   31.8207,   74.9126  },
       {  1.88e+04,  2.74e+04,   2.2,   2.39,    1.7953,   4.77760,   31.2932,   77.7668  },
       {  1.88e+04,  2.74e+04,  2.39,    2.6,    1.8072,    4.9156,   32.5164,   74.8536  },
       {  1.88e+04,  2.74e+04,   2.6,   2.85,    1.8900,   5.69720,   34.2404,   78.4988  },
       {  1.88e+04,  2.74e+04,  2.85,  100.0,    1.8696,  8.760800,   34.1142,  85.85300  },
       {  2.74e+04,  3.44e+04,  0.00,    1.5,  1.710600,    3.3300,   24.0000,   63.0729  },
       {  2.74e+04,  3.44e+04,   1.5,    1.7,  1.738500,    3.5660,   25.9400,   65.2095  },
       {  2.74e+04,  3.44e+04,   1.7,   1.84,    1.7400,    3.7758,   26.0000,   66.9902  },
       {  2.74e+04,  3.44e+04,  1.84,   1.95,    1.7761,    4.0634,  30.19602,   69.4300  },
       {  2.74e+04,  3.44e+04,  1.95,   2.05,  1.800000,    4.1365,   31.2788,   70.1575  },
       {  2.74e+04,  3.44e+04,  2.05,    2.2,    1.7890,    4.3900,   32.6145,   71.8655  },
       {  2.74e+04,  3.44e+04,   2.2,   2.39,    1.8400,   4.74270,   32.4424,   75.3111  },
       {  2.74e+04,  3.44e+04,  2.39,    2.6,    1.9065,    4.8100,   33.7970,   74.2445  },
       {  2.74e+04,  3.44e+04,   2.6,   2.85,    1.8800,   5.24840,   34.2832,   76.8384  },
       {  2.74e+04,  3.44e+04,  2.85,  100.0,    1.9028,  6.097200,   34.9460,  80.14040  },
       {  3.44e+04,   4.2e+04,  0.00,    1.5,  1.686024,    3.6700,   27.2464,   66.3094  },
       {  3.44e+04,   4.2e+04,   1.5,    1.7,  1.721800,    3.8282,   28.9360,   67.8164  },
       {  3.44e+04,   4.2e+04,   1.7,   1.84,    1.8100,    3.7514,   31.8559,   66.9570  },
       {  3.44e+04,   4.2e+04,  1.84,   1.95,    1.8564,    4.0028,  33.55480,   68.5468  },
       {  3.44e+04,   4.2e+04,  1.95,   2.05,  1.780000,    4.2300,   31.7022,   70.9715  },
       {  3.44e+04,   4.2e+04,  2.05,    2.2,    1.8300,    4.1627,   32.8963,   69.3542  },
       {  3.44e+04,   4.2e+04,   2.2,   2.39,    1.7325,   4.97000,   29.0075,   75.2060  },
       {  3.44e+04,   4.2e+04,  2.39,    2.6,    1.8788,    4.8880,   33.7588,   74.6956  },
       {  3.44e+04,   4.2e+04,   2.6,   2.85,    1.9000,   5.23860,   34.6986,   76.5244  },
       {  3.44e+04,   4.2e+04,  2.85,  100.0,    1.9330,  8.224500,   36.9315,  84.36500  },
       {   4.2e+04,  5.11e+04,  0.00,    1.5,  1.780000,    3.5413,   30.0000,   65.2878  },
       {   4.2e+04,  5.11e+04,   1.5,    1.7,  1.820000,    3.9700,   32.2000,   69.0080  },
       {   4.2e+04,  5.11e+04,   1.7,   1.84,    1.8115,    4.3635,   32.3015,   71.1800  },
       {   4.2e+04,  5.11e+04,  1.84,   1.95,    1.8300,    4.2600,  32.05800,   71.1768  },
       {   4.2e+04,  5.11e+04,  1.95,   2.05,  1.843500,    4.4195,   32.8500,   71.9090  },
       {   4.2e+04,  5.11e+04,  2.05,    2.2,    1.8100,    4.4372,   32.0028,   71.5700  },
       {   4.2e+04,  5.11e+04,   2.2,   2.39,    1.9200,   4.84480,   35.9152,   74.5596  },
       {   4.2e+04,  5.11e+04,  2.39,    2.6,    1.8845,    5.3155,   33.9685,   76.7660  },
       {   4.2e+04,  5.11e+04,   2.6,   2.85,    1.9020,   5.55600,   35.5280,   77.1000  },
       {   4.2e+04,  5.11e+04,  2.85,  100.0,    1.9415,  7.284013,   34.6020,  81.65500  },
       {  5.11e+04,  6.42e+04,  0.00,    1.5,  1.700000,    3.6156,   27.8480,   65.6900  },
       {  5.11e+04,  6.42e+04,   1.5,    1.7,  1.849500,    3.7800,   33.5090,   67.3855  },
       {  5.11e+04,  6.42e+04,   1.7,   1.84,    1.8200,    3.8257,   32.7446,   67.8713  },
       {  5.11e+04,  6.42e+04,  1.84,   1.95,    1.8716,    4.0528,  33.15760,   69.4152  },
       {  5.11e+04,  6.42e+04,  1.95,   2.05,  1.840000,    4.3980,   33.1240,   71.4880  },
       {  5.11e+04,  6.42e+04,  2.05,    2.2,    1.8800,    4.5000,   33.9714,   72.1404  },
       {  5.11e+04,  6.42e+04,   2.2,   2.39,    1.9000,   5.23900,   34.4998,   76.5116  },
       {  5.11e+04,  6.42e+04,  2.39,    2.6,    1.9000,    5.0493,   35.0136,   74.7660  },
       {  5.11e+04,  6.42e+04,   2.6,   2.85,    1.9236,   5.62880,   35.4900,   78.8108  },
       {  5.11e+04,  6.42e+04,  2.85,  100.0,    1.9586,  7.200000,   37.2134,  82.79360  },
       {  6.42e+04,  8.21e+04,  0.00,    1.5,  1.820000,    3.4360,   32.5140,   64.5964  },
       {  6.42e+04,  8.21e+04,   1.5,    1.7,  1.865500,    3.5590,   33.1650,   65.4520  },
       {  6.42e+04,  8.21e+04,   1.7,   1.84,    1.8400,    4.1050,   33.1175,   69.8100  },
       {  6.42e+04,  8.21e+04,  1.84,   1.95,    1.8700,    4.0378,  33.82170,   69.3137  },
       {  6.42e+04,  8.21e+04,  1.95,   2.05,  1.890000,    4.4169,   34.7257,   72.5607  },
       {  6.42e+04,  8.21e+04,  2.05,    2.2,    1.8900,    4.5100,   34.1900,   72.7700  },
       {  6.42e+04,  8.21e+04,   2.2,   2.39,    1.9000,   4.68540,   35.2306,   73.7454  },
       {  6.42e+04,  8.21e+04,  2.39,    2.6,    1.9100,    4.8172,   35.7166,   74.4903  },
       {  6.42e+04,  8.21e+04,   2.6,   2.85,    1.8900,   5.76525,   35.1765,   78.5360  },
       {  6.42e+04,  8.21e+04,  2.85,  100.0,    1.9850,  7.535000,   36.9050,  83.21500  },
       {  8.21e+04,  1.08e+05,  0.00,    1.5,  1.910000,    3.5149,   35.0754,   65.0292  },
       {  8.21e+04,  1.08e+05,   1.5,    1.7,  1.880000,    3.4740,   34.3864,   64.5440  },
       {  8.21e+04,  1.08e+05,   1.7,   1.84,    1.8860,    3.6180,   34.1500,   66.6400  },
       {  8.21e+04,  1.08e+05,  1.84,   1.95,    1.8700,    3.9480,  34.50550,   68.5495  },
       {  8.21e+04,  1.08e+05,  1.95,   2.05,  1.940000,    3.9684,   36.1160,   69.0188  },
       {  8.21e+04,  1.08e+05,  2.05,    2.2,    1.9000,    4.5251,   34.4200,   72.8006  },
       {  8.21e+04,  1.08e+05,   2.2,   2.39,    1.9228,   4.65680,   35.3936,   72.6636  },
       {  8.21e+04,  1.08e+05,  2.39,    2.6,    1.9738,    5.3662,   36.9628,   75.9848  },
       {  8.21e+04,  1.08e+05,   2.6,   2.85,    1.9400,   5.69980,   36.0980,   78.5082  },
       {  8.21e+04,  1.08e+05,  2.85,  100.0,    2.0400,  6.691400,   33.5047,  80.93520  },
       {  1.08e+05,  1.67e+05,  0.00,    1.5,  1.829900,    3.2014,   32.6000,   61.0213  },
       {  1.08e+05,  1.67e+05,   1.5,    1.7,  1.899200,    3.5104,   35.6696,   64.6960  },
       {  1.08e+05,  1.67e+05,   1.7,   1.84,    1.9029,    4.4000,   35.6379,   72.0400  },
       {  1.08e+05,  1.67e+05,  1.84,   1.95,    1.8704,    4.3992,  33.70860,   71.8052  },
       {  1.08e+05,  1.67e+05,  1.95,   2.05,  1.897500,    4.3925,   34.8900,   71.3600  },
       {  1.08e+05,  1.67e+05,  2.05,    2.2,    1.9100,    4.2988,   35.6012,   72.2694  },
       {  1.08e+05,  1.67e+05,   2.2,   2.39,    1.9318,   4.68740,   35.9260,   72.9028  },
       {  1.08e+05,  1.67e+05,  2.39,    2.6,    1.9176,    4.7434,   35.9182,   73.7240  },
       {  1.08e+05,  1.67e+05,   2.6,   2.85,    1.9582,   5.02880,   36.7668,   75.3236  },
       {  1.08e+05,  1.67e+05,  2.85,  100.0,    2.0900,  6.312400,   40.4518,  81.00000  },
       {  1.67e+05,  9.99e+10,  0.00,    1.5,  1.886300,    3.7614,   31.0000,   59.6822  },
       {  1.67e+05,  9.99e+10,   1.5,    1.7,  1.852200,    4.0452,   33.2568,   69.6251  },
       {  1.67e+05,  9.99e+10,   1.7,   1.84,    1.9196,    4.0804,   31.0000,   69.8354  },
       {  1.67e+05,  9.99e+10,  1.84,   1.95,    1.9569,    4.1655,  36.87380,   70.2727  },
       {  1.67e+05,  9.99e+10,  1.95,   2.05,  1.935800,    4.4042,   35.9066,   71.8756  },
       {  1.67e+05,  9.99e+10,  2.05,    2.2,    1.9900,    4.5898,   37.1896,   72.7968  },
       {  1.67e+05,  9.99e+10,   2.2,   2.39,    1.9913,   4.83740,   38.2817,   74.5870  },
       {  1.67e+05,  9.99e+10,  2.39,    2.6,    2.0000,    4.7601,   38.4896,   74.1303  },
       {  1.67e+05,  9.99e+10,   2.6,   2.85,    2.0600,   5.16950,   39.1830,   75.1746  },
       {  1.67e+05,  9.99e+10,  2.85,  100.0,    2.1000,  6.890000,   40.8130,  82.07220  }
};

#endif
