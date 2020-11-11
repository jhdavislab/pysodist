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

        subroutine calc_model

        include 'isodist.h'

c       calculate the isotope distribution for the species
c       and sum them up

c       initialize spectrum array
        call init_c(ceft(1),ncp,(0.0,0.0))
        do 360 i=1,nfitpar
360     continue
        
c       apply the gaussian and offset

        do 320 i=1,n_species
        do 321 j=1,ncp
        cseft(j,i)=cseft(j,i)*cgft(j)*cpft(j)           
321     continue
320     continue


c       multiply species by the amplitudes
c       sum up the spectra for the species

        do 322 i=1,n_species
        do 323 j=1,ncp
        ceft(j)=ceft(j)+amp(i)*cseft(j,i)
323     continue
322     continue

c       IFT final spectrum 

        call bogus_ft(eft,np,-1)
        
c       add the baseline
        
        do 342 j=1,np
        eft(j)=eft(j)+b
342     continue

        return
        end
