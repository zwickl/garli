Garli version 0.95 - October 2006

Compilation help:
-to make the serial version type "make s" (or simply "make")
-If using gcc version 2.95, be sure to change which CC_FLAGS line you use in the Makefile.  The optimization in that version breaks the code.
-icc tends to make faster executables that gcc on intel hardware, so choose a compiler accordingly
-There have been reports of problems (crashes) when using gcc version 4.0 with the GARLI code, but I have not verified them.
-I may be able to provide compilation help if you've exhausted all other resources.

SAMPLE DATASET:
-The Rana dataset is from Hillis and Wilcox. 2005. Phylogeny of the New World true frogs (Rana). Mol Phylogenet Evol. 34(2):299-314. 	
-The 2 best trees for this dataset under the GTR+I+G model have log likelihood scores of -21812.66941 and -21812.64132. See if you can find them...  
-A sample starting tree is contained in ranastart.tre
-A sample constraint tree is contained in ranaconstraint.tre

Disclaimers:
This software is copyright Derrick Joel Zwickl, 2006.
It is made available without warranty:
BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY FOR THE PROGRAM,
TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN OTHERWISE STATED IN WRITING THE
COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY OF
ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS
TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE 
DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION. 

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL ANY COPYRIGHT
HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR REDISTRIBUTE THE PROGRAM AS PERMITTED
ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR
CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING
BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY
YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), 
EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
