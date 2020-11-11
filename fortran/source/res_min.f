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

        subroutine res_min(nmz2,na,afit,ycalc,iflag)

c       chi squared calculation for minimization

c       this is routine "fcn" called by lmdif1 from minpack

        include 'isodist.h'

        integer iflag,na
        real afit(na)
        real ycalc(nmz2)
                
c       handle fixed/free parameters

        nna=1
        do 350 k=1,nfitpar
        if(fit_schedule(nctr,k).eq.1) then
                x(k)=afit(nna)
                nna=nna+1
        endif
350     continue

c       fix parameters gone amok

c       baseline
        if(x(1).lt.0)then
                x(1)=-x(1)
        endif
        
c       offset
        if(x(2).gt.1.0.or.x(2).lt.-1.0)then
                x(2)=0.0
        endif

c       gaussian width
        if(x(3).lt.0)then
                x(3)=-x(3)
        endif

c       amplitude parameters
        do 330 i=1,n_species
        if(x(amp_idx(i)).lt.0.0)then
                x(amp_idx(i))=-x(amp_idx(i))
        endif
330     continue

c       frac parameters
        do 331 i=1,n_var_atom
        if(x(frac_idx(i)).lt.0.0)then
                x(frac_idx(i))=0.1
        endif
        if(x(frac_idx(i)).gt.1.0)then
                x(frac_idx(i))=0.9
        endif
331     continue

c       phi parameters
        do 332 i=1,n_var_res
        if(x(phi_idx(i)).lt.0.0)then
                x(phi_idx(i))=0.1
        endif
        if(x(phi_idx(i)).gt.1.0)then
                x(phi_idx(i))=0.9
        endif
332     continue

        call calc_iso_dist

c       calculate residuals and chisq

        chisq=0.0
        do 300 i=1,nmz2
        ii=mz_ptr(i)
        ycalc(i)=eft_obs(ii)-eft(ii)
        chisq=chisq+ycalc(i)*ycalc(i)
300     continue

        if(iflag.eq.0) write(6,678)chisq,(x(j),j=1,nfitpar)
678     format(4x,e12.4,1x,4(f10.3,1x),f12.6,1x,f10.3)

        chisq_end=chisq

        return
        end
        
        
