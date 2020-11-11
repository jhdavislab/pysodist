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

        subroutine setup_model
                
        include 'isodist.h'

        character*4 alab,flab,plab
        
c       Universal fitting parameters:
c               x(1)=b
c               x(2)=off
c               x(3)=gw

        parlabel(1)="B"
        parlabel(2)="OFF"
        parlabel(3)="GW"
        
        nfitpar=3

c       amplitude parameters

        alab="AMP_"
        do 300 i=1,n_species
        nfitpar=nfitpar+1
        amp_idx(i)=nfitpar
        parlabel(nfitpar)=alab//species_lab(i)
        x_init(amp_idx(i))=amp(i)
300     continue

c       fractional atom parameters

        flab="FRC_"
        do 301 i=1,n_var_atom
        nfitpar=nfitpar+1
        frac_idx(i)=nfitpar
        parlabel(nfitpar)=flab//atom_symb(var_atom(i))
        x_init(frac_idx(i))=fiso(var_atom(i),niso(var_atom(i)))
301     continue

c       fractional residue parameters

        plab="PHI_"
        do 302 i=1,n_var_res
        nfitpar=nfitpar+1
        phi_idx(i)=nfitpar
        parlabel(nfitpar)=plab//residue(var_res(i))
        x_init(phi_idx(i))=phi(i)
302     continue

c       assemble fit schedule

c       flags for fixed universal parameters
        nb=fixed_pars(1)
        no=fixed_pars(2)
        ng=fixed_pars(3)

c       fit baseline and amplitudes
        nround=1
        call fit_sched(nround,nb,0,0,1,0,0)

c       fit baseline and gaussian       
        if(ng.ne.0)then
                nround=nround+1
                call fit_sched(nround,nb,0,1,0,0,0)
        endif

c       fit 3 universal params
        if(no.ne.0)then
                nround=nround+1
                call fit_sched(nround,nb,1,ng,0,0,0)
        endif

c       fit 3 universal, and fractional residues
        if(n_var_res.gt.0)then
                nround=nround+1
                call fit_sched(nround,nb,no,ng,0,0,1)
        endif

c       fit 3 universal and fractional atoms    
        if(n_var_atom.gt.0)then
                nround=nround+1
                call fit_sched(nround,nb,no,ng,0,1,0)
        endif

c       fit all params
        nround=nround+1
        call fit_sched(nround,nb,no,ng,1,1,1)

c       write out model header for csv output file
        call model_header
        return
        end
