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

        subroutine calc_residues(ifix)
                
        include 'isodist.h'
        
        integer ifix
        
c       ifix = 0  fixed part of residues, called at startup
c       ifix = 1  variable part of residues, called at each iteration

c       first calculate variable atoms

        do 300 ni=1,n_species
        do 301 j=1,n_res_lib

c       initialize fixed part of all residues on first pass
        
        if(ifix.eq.0) call init_c(ceft_fix(1,j,ni),ncp,(1.0,0.0))

c       skip inactive residues
        if(ifix.eq.1.and.active_res(j).eq.0) goto 301

c       initialize varialble part of all residues       
        if(ifix.eq.1) call init_c(ceft_var(1,j,ni),ncp,(1.0,0.0))

        do 302 k=1,natom_types
                kidx=atom_ptr_list(k)
                n=res_atm_comp(j,k,ni)  
                do 303 jn=1,n
                do 305 jc=1,ncp
                        if(ifix.eq.0.and.fixed_var_atom(k,j,ni).eq.0)then
                                ceft_fix(jc,j,ni)=ceft_fix(jc,j,ni)*cel_ft(jc,kidx)
                        endif
                        if(ifix.eq.1.and.fixed_var_atom(k,j,ni).eq.1) then
                                ceft_var(jc,j,ni)=ceft_var(jc,j,ni)*cel_ft(jc,kidx)
                        endif
305             continue
303             continue
302     continue

301     continue
300     continue

c       calculate residue spectra

        if(ifix.eq.1.)then
                do 320 j=1,n_res_lib
                if(active_res(j).eq.0)goto 320
                do 321 ni=1,n_species           
                do 322 jc=1,ncp
                ceft_res(jc,j,ni)=ceft_fix(jc,j,ni)*ceft_var(jc,j,ni)
322             continue
321             continue
320             continue
        endif

c       calculate fixed residues
        do 330 i=1,n_fix_res
        j=fix_res(i)
        do 331 jc=1,ncp
        ceft_res(jc,j,2)=(1.0-phi_fix(i))*ceft_res(jc,j,1)+phi_fix(i)*ceft_res(jc,j,2)
331     continue
330     continue
        
c       calculate variable residues
        
        do 310 i=1,n_var_res
        j=var_res(i)
        do 311 jc=1,ncp
        ceft_res(jc,j,2)=(1.0-phi(i))*ceft_res(jc,j,1)+phi(i)*ceft_res(jc,j,2)
311     continue
310     continue
        
        return
        end
