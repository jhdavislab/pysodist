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

        subroutine calc_atoms
        
        include 'isodist.h'
        
c       recalculate u-domain distribution for variable atoms

        do 300 i=1,n_var_atom
        fiso(var_atom(i),1)=1.0-frac(i)
        fiso(var_atom(i),2)=frac(i)
        call iso_ft(var_atom(i))
300     continue

        return
        end
