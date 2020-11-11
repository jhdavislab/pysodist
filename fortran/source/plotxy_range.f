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

        subroutine plotxy_range(dat,n,nfr,nto)
        
        include 'isodist.h'
        
        real dat(n)

        open(unit=7,file=outfile,status="unknown",form="formatted")
        
        write(7,600)(real(i)/scale_mz+mz_hd,dat(i),i=nfr,nto)
600     format(e16.8,",",e16.8)

        close(unit=7)

        return
        end
        
        

