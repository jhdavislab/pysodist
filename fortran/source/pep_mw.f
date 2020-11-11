C Copyright 2008 James R. Williamson and Michael T. Sykes
C
C http://williamson.scripps.edu/isodist/
C
C This file is part of isodist.
C 
C isodist is free software: you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published by
C the Free Software Foundation, either version 3 of the License, or
C (at your option) any later version.
C
C isodist is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C GNU General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with isodist.  If not, see <http://www.gnu.org/licenses/>.

        subroutine pep_mw

        include 'isodist.h'

        integer res_ptr

c       pep_length=lengthn(40,pepseq)

c       calculate molecular weight from peptide sequence

        do 330 ni=1,n_species
        
        mw=0.0

        do 310 i=1,natom_types
        mol_form(i)=0
310     continue

        do 300 i=1,pep_length

        j=match_residue(pepseq_res(i),n_res_lib,residue)
        
        do 302 k=1,natom_types
        kidx=atom_ptr_list(k)
        n=res_atm_comp(j,k,ni)
        mw=mw+float(n)*miso(kidx,1)
        mol_form(k)=mol_form(k)+n
302     continue

300     continue

        write(6,600)mw
600     format("Molecular Weight = ",f12.4)

        write(6,601)ni,(atom_symb(atom_ptr_list(i)),mol_form(i),i=1,natom_types)
601     format("species ",i2," molecular formula = ",10("(",a2,i3,")"))

330     continue

        return
        end
        

