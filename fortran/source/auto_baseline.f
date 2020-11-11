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

        subroutine auto_baseline
        
        include 'isodist.h'
        
c       autobasline routine
c       ad libitum algorithm by jrw
c       sort data, then average 1st and 2nd quartile values.
c       seems to work ok, but is very heuristic
        
        do 300 i=1,nmz
        ysort(i)=yobs(i)
300     continue

        call sort(ysort,ysort2,nmz)
        
c       sort gives increasing order
c       get 1st and 2nd quartile values
        
        nq1=int(float(nmz)/4.0+0.5)
        nq2=int(float(nmz)/2.0+0.5)
        
        x_init(1)=(ysort(nq2)+ysort(nq1))/2.0
        
        return
        end


                
