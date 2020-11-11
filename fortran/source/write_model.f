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

        subroutine write_model
        
        include 'isodist.h'
        
c       find rightmost dot

        do 300 i=70,1,-1
        if(infile(i:i).eq.".") then
                nd=i-1
                goto 900
        endif
300     continue
900     continue



        write(18,1807)infile(1:nd),protlab,pepseq(1:pep_len),mw,nz,ptim,chisq_end,psymb,pmz,(x(i),i=1,nfitpar)
1807    format(a,",",a2,",",a,",",f8.3,",",i2,",",a,",",e12.5,",",i2,",",f8.3,8(",",e12.5))
        
        
        return
        end
