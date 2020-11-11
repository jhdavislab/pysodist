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

        subroutine calc_species
        
        include 'isodist.h'
        
c       calculate the final residue-based spectra
c               multiplying the u-domain functino for each residue

        do 320 ni=1,n_species

        call init_c(cseft(1,ni),ncp,(1.0,0.0))  

        do 322 n=1,pep_length
        do 323 j=1,n_res_lib
        if(residue(j).eq.pepseq(n:n))then
                do 300 nn=1,ncp
                cseft(nn,ni)=cseft(nn,ni)*ceft_res(nn,j,ni)
300             continue
        endif
                
323     continue
322     continue
320     continue

        return
        end
        
        
