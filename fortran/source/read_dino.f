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

        subroutine read_dino
        
        include 'isodist.h'
        
        write(6,*) infile
        open(unit=2,file=infile,status='old',form='formatted')
        
        nmz=0
        mzfr=np
        mzto=0
        ifirst=0
        
800     read(2,100,end=900)rec
100     format(a200)
        call split(' ')
        rmz=f_ex(1)
        ri=f_ex(2)

c       use first point as heterodyne mass

        if(ifirst.eq.0)then
                mz_hd=float(nz)*rmz
                ifirst=1
        endif

        n=int((rmz*float(nz)-mz_hd)*scale_mz+0.5)

        if(n.gt.np) then
                write(6,600)n,np,scale_mz
                goto 900
        endif
        
600     format(' PROBLEM WITH MZ DATA:  n = ',i8,'   np = ',i8,' scale_mz = ',f9.3)

        if(n.eq.0) goto 800
        
        nmz=nmz+1

        if(n.lt.mzfr) mzfr=n
        if(n.gt.mzto) mzto=n
        mz(n)=rmz*float(nz)-mz_hd
        eft_obs(n)=ri
        mz_ptr(nmz)=n
        xobs(nmz)=(rmz*float(nz)-mz_hd)*scale_mz
        yobs(nmz)=ri
        goto 800
900     continue

        write(6,601)mz_hd
601     format(' heterodyne mass = ',f9.3)

c       this section reads in the whole peak file in the case where the calc mass does not observe the input mass
c               this is usually because the peptide is modified.  A crude hack to prevent crashing in batch mode.
        if(mz_hd.eq.0.0)then
                bogusfile = 'bogus.dat'
                open(unit=22,file=bogusfile,status='unknown',form='formatted')
                write(22,660)
                write(6,660)
660             format(' BOGUS PEPTIDE...possible modification???')
                close(unit=22)
                stop
        endif

        return
        end
        
