// GARLI Version 0.96b8 (May 2008)
// Copyright 2005-2008 Derrick J. Zwickl
// email: zwickl@nescent.org
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

THIS IS A BETA VERSION - CHECK RESULTS CAREFULLY!  
Please let me know of any problems, concerns or feedback (garli.support@gmail.com)

GARLI version 0.96 is much exhanced and expanded relative to the previous
version 0.951.  It should replace earlier versions, and should be backwards
compatible with all configuration files and datasets that were used with the
previous version.

->See the manual or support website (http://www.nescent.org/wg_garli)
 for detailed information on using the program.

->Example datasets and configuration files files can be found in the example folder

->For compilation help see the INSTALL file.  The program now comes with a full
 configure script for *nix configuration and installation.

New in version 0.96
1.Rigorous reading of Nexus datasets using Paul Lewis and Mark Holder's Nexus
 Class Library (NCL).
2.Ability to read Nexus starting trees using NCL.
3.Ability to perform inference under amino acid and codon-based models
 of sequence evolution (datatype = aminoacid, datatype = codon).
4.Ability to specify multiple search replicates in a single config file (searchreps = #).
5.Ability to specify outgroups for orientation of inferred trees (outgroup = # # #).
6.Ability to use backbone as well as normal topological constraints.
7.Ability to create fast likelihood stepwise addition starting trees
 (streefname = stepwise).
8.MPI version that spreads a specified number of serial runs across processors using
 a single config file, writing output to different output files (for example, to do
 25 bootstrap replicates simultaneously on each of 8 processors).
9.Ability to perform nucleotide inference using any sub-model of the General
 Time-Reversible model (GTR), in addition to all of the common “named” models (K2P, HKY, etc).
10.Speed increases for non-parametric bootstrapping

Condensed summary of new model settings (for more detailed descriptions and 
for unchanged settings see the manual or support webpage):
datatype = {nucleotide, aminoacid, codon, codon-aminoacid} - These set the type of model
 to be used.  Note that aminoacid and codon models are MUCH slower than nucleotide models.
-"aminoacid" is for datasets consisting of the 20 aminoacid single letter codes. 
-"codon" is for dna data (aligned in frame!) to be analyzed using a 60-62 state model
 that incorporates both the nucleotide substitution process and information on the
 genetic code.  This involves the estimation of at least one dN/dS ratio
 (aka nonsynonymous/synonymous rate ratio, or omega or w).  This is essentially
 the Goldman-Yang 1994 model and other related models.
-"codon-aminoacid" is for aligned dna sequences that are translated to aminoacids
 and analyzed under an aminoacid model

The different datatypes have different allowable model settings, listed here.
For nucleotide data:
ratematrix = {6rate, 2rate, 1rate, fixed, (a b c d e f) }
statefrequencies = {estimate, empirical, equal, fixed}
ratehetmodel = {none, gamma, gammafixed}
numratecategories = {#} (not including invariant site class, must be 1 for ratehetmodel = none)
invariantsites = {estimate, none}

For aminoacid or codon-aminoacid data:
ratematrix = {poisson, dayhoff, jones, wag, mtmam, mtrev}
statefequencies = {equal, dayhoff, jones, wag, mtmam, mtrev, empirical, estimate}
ratehetmodel = {none, gamma, gammafixed}
numratecategories = {#} (not including invariant site class, must be 1 for ratehetmodel = none)
invariantsites = {estimate, none, fixed}

For codon data:
ratematrix = {6rate, 2rate, 1rate, fixed, (a b c d e f) }
 (this is the nucleotide substitution process assumed by the codon model)
statefrequencies = {empirical, equal, F1x4, F3x4} (F1x4 and F3x4 are PAML's
 terminology, and calculate the codon frequencies as the product of the total
 nucleotide frequencies or the nucleotide frequencies at each codon position, respectively)
ratehetmodel = {none, nonsynonymous} {nonsynonymous estimates multiple dN/dS categories
 at the proportion of sites belonging to each.  This is the M3 model of PAML)
numratecategories = {#} (the number of dN/dS categories)
invariantsites = {none} (not allowed in codon model)

For codon or codon-aminoacid:
Geneticcode = {standard, vertmito, invertmito}

// GARLI Vesion 0.96b8
// Copyright 2005-2008 Derrick J. Zwickl
// email: zwickl@nescent.org
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
