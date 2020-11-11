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

        subroutine fit_sched(j,nbase,noff,nwid,namp,nfrac,nphi)

c       routine to fill up the fit_schedule array based on a set of flags
c               0 = parameter held constant
c               1 = parameter varied
c       these values are passed to marqmin to fix or free parameters

c               j = round #
c               nbase = baseline flag
c               noff = offset flag
c               nwid = gaussian width flag
c               namp = amplitude parameters flag
c               nfrac = fractional atoms flag
c               nphi = fractional residues flag

        include 'isodist.h'
                
        fit_schedule(j,1)=nbase
        fit_schedule(j,2)=noff
        fit_schedule(j,3)=nwid
        
        do 310 i=4,nfitpar
        fit_schedule(j,i)=0
310     continue
        
        if(namp.eq.1)then
                do 300 i=1,n_species
                fit_schedule(j,amp_idx(i))=1
300             continue
        endif
        
        if(nfrac.eq.1)then
                do 301 i=1,n_var_atom
                fit_schedule(j,frac_idx(i))=1
301             continue
        endif

        if(nphi.eq.1)then
                do 302 i=1,n_var_res
                fit_schedule(j,phi_idx(i))=1
302             continue
        endif

        return
        end
