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

        subroutine setup_atoms

c       masses and natural abundances taken from 
c               iMass v. 1.0 by Urs Roethlisberger,
c               who references CRC Handbook, 74th Ed. (1993)

        include 'isodist.h'

        do 300 i=1,nelements
        atom_symb(i)=" "
300     continue

c       read in atom definition file
c       field 1 = atomic number
c       field 2 = atomic symbol
c       field 3 = # isotopes
c       field 4 = comment (read, not used)

        open(unit=88, file=atomfile, status="old", form="formatted")
        write(6,*)atomfile
        
        na=0
900     read(88,880,end=800)rec
880     format(a200)
        na=na+1
        call split(' ')
        nat=i_ex(1)

        atom_num(na)=nat
        atom_symb(nat)=rec(nf(2):nl(2))
        niso(nat)=i_ex(3)

        atom_comment(nat)=rec(nf(4):nl(nstr))

        do 301 i=1,niso(nat)
        read(88,880)rec
        call split(' ')
        miso(nat,i)=f_ex(1)
        fiso(nat,i)=f_ex(2)
301     continue

        goto 900
        
800     continue

        if(na.gt.natoms)then
        write(6,600)natoms,na
600     format("ERROR: expecting ",i3," atoms, read in ",i3,".  Check atom_definitions.txt and header file")
        stop
        endif

        return
        end
        
