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

        subroutine iso_ft(i)
        
        include 'isodist.h'

c       i is the atom index

c       initialize the element array
                
        do 301 j=1,np
        el_ft(j,i)=0.0
301     continue
        
        do 302 j=1,niso(i)
        n=int(miso(i,j)*scale_mz+0.5)+1
        el_ft(n,i)=fiso(i,j)
302     continue

        call bogus_ft(el_ft(1,i),np,1)
        
c       set imaginary part of 1st point to 0
c               this gives baseline offsets when raised to higher powers
c       Note:  this is due to paciing values for F(0) and F(n/2) into the 
c               first complex point. Since this gets multiplied by a gaussian
c               later, it doesn't matter...(I think).

        el_ft(2,i)=0.0

        return
        end
        
        
