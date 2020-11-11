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

        subroutine setup_residues

        include 'isodist.h'
        
        character*1 bl
        character*5 fixed_str
        
c       read in residue definition file

        open(unit=88, file=resfile, status="old", form="formatted")
        write(6,*)resfile
        ni=0
        
c       read in comment line
        read(88,880,end=800)rec
880     format(a200)

c       read in number of species
        read(88,880)rec
        call split(' ')
        n_species=i_ex(1)
        
c       read in species labels and initial amplitudes
        do 300 i=1,n_species
        read(88,880)rec
        call split(' ')
        species_lab(i)=rec(nf(1):nl(1))
        amp(i)=f_ex(2)
300     continue

c       read in atom type defs 
c       identify variable atoms
c           read in initial fraction
        n_var_atom=0
        read(88,880)rec
        call split(' ')
        natom_types=i_ex(1)
        
        do 301 j=1,natom_types
        read(88,880)rec
        call split(' ')
        nat=match_atom(rec(nf(1):nl(1)),nelements,atom_symb)
        atom_ptr_list(j)=nat
        if(rec(nf(2):nl(2)).eq."default")atom_fix(nat)=0
        if(rec(nf(2):nl(2)).eq."fixed") then
                atom_fix(nat)=0
                fiso(nat,niso(nat))=f_ex(3)
                fiso(nat,1)=1.0-f_ex(3)
        endif                           
        if(rec(nf(2):nl(2)).eq."variable")then
                atom_fix(nat)=1
                n_var_atom=n_var_atom+1
                var_atom(n_var_atom)=nat
                fiso(var_atom(n_var_atom),niso(var_atom(n_var_atom)))=f_ex(3)
                fiso(var_atom(n_var_atom),1)=1.0-f_ex(3)
        endif
        
301     continue

C       Calculate mu-domain function for each atom that is used

        do 310 i=1,natom_types  
        call iso_ft(atom_ptr_list(i))
310     continue

c       read in residue definitions

        nr=0
        n_var_res=0
        n_fix_res=0
900     read(88,880,end=800)rec
        call split(' ')
        nr=nr+1
        residue(nr)=rec(nf(1):nl(1))

        if(rec(nf(2):nl(2)).eq."default")then
                fixed_var_res(nr)=0
        endif
        if(rec(nf(2):nl(2)).eq."fixed")then
                fixed_var_res(nr)=1
                n_fix_res=n_fix_res+1
                fix_res(n_fix_res)=nr
                phi_fix(n_fix_res)=f_ex(3)
        endif
        if(rec(nf(2):nl(2)).eq."variable")then
                fixed_var_res(nr)=2
                n_var_res=n_var_res+1
                var_res(n_var_res)=nr
                phi(n_var_res)=f_ex(3)
        endif

c       read in atom composition of the residue for each species        
        do 302 j=1,n_species
        read(88,880)rec
        call split(' ')
        do 303 k=1,natom_types
        res_atm_comp(nr,k,j)=i_ex(k)
        fixed_var_atom(k,nr,j)=atom_fix(atom_ptr_list(k))
303     continue
302     continue

        goto 900
800     continue

        return
        end
        
