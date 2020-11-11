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

        subroutine bogus_ft(xd,nd,idir)

c       perform real fourier transform, both forward and back
c       idir = 0 initializes working array
c       idir = 1 forward transform
c       idir = -1 back transform
c
c       uses a subset of routines from fftpack5:
c               rfft1i,rfft1f, rfft1b and dependencies

c       routine is supposed to substitute for realft from Numerical Recipes
c               using public domain source code

c       fftpack5 realft differs from realft by a (nearly) circular permutation,
c               and also how it deals with the first and last points
c               so bogus shuffling was introduced to allow substitution
c               of this routine.
c
c       The isodist calculations multiply u-domain functions as complex numbers xc(i),i=1,ncp
c       which is equivalent to the real array x(j),j=1,np, where np=2*ncp.  The fftpack5 routine
c       RFFT1F returns values where x(1) is a real value, (x(2*i),x(2*i+1)), i=1,ncp-1 are complex numbers,
c       and x(np) is a real value.  The shuffling orders/reorders the elements between the u-domain order
c       and the fftpack order.  Yuk. Highly inelegant.

c       Life is short, computers are fast, and this code is free.

        include 'isodist.h'

        integer nd,idir
        real xd(nd)
        real fn2

        fn2=float(nd/2)

c       initialize working array
        
        if(idir.eq.0) call rfft1i(nd,wsave,nwork2,ier)  

c       FORWARD XFORM

        if(idir.eq.1)then

        call rfft1f(nd,1,xd,nd,wsave,nwork2,work,nwork2,ier)
        
c       shuffle elements

        temp=xd(nd)*fn2
        
        do 300 j=nd,3,-1
        xd(j)=xd(j-1)*fn2
300     continue

        xd(1)=xd(1)*2.0*fn2
        xd(2)=temp
        
        endif
        
c       BACK XFORM
        
        if(idir.eq.-1)then
        
c       unshuffle

        temp=xd(2)

        do 311 j=2,nd-1
        xd(j)=xd(j+1)
311     continue
        
        xd(nd)=temp/2.0
        xd(1)=xd(1)/2.0
        
        call rfft1b(nd,1,xd,nd,wsave,nwork2,work,nwork2,ier)
        
        endif
        
        return  
        end
