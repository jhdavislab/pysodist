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

        subroutine model_header

c       write out the header containing column labels for the output file
        
        include 'isodist.h'
                
        outlabels="file,protein,pep,mw,z_charge,tim,chisq,symb,mz"
                
        do 300 i=1,nfitpar
        n=lengthn(200,outlabels)        
        outlabels=outlabels(1:n)//","//parlabel(i)
300     continue

        write(18,180)outlabels
180     format(a200)

        return
        end
        
        
