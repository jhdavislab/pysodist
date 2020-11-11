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

        subroutine char_num(c,i)

c       subroutine to decide if a character is alpha or numeric
c               for alphanumeric string handling
c          ...grim fortran hack
        
        character*1 c
        character*52 alphabet
        character*10 numbers
        
        i=0

        alphabet="abcdefghijklmnnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
        numbers="0123456789"
        
        do 300 j=1,52
        if(alphabet(j:j).eq.c)i=1
300     continue
        do 301 j=1,10
        if(numbers(j:j).eq.c)i=2
301     continue

        return
        end
        
        real function f_ex(nfld)
        
        include 'isodist.h'

        character*1 c
        character*3 cxp
        real mult

c       floating point extraction
c       works in conjunction with the split routine
c       after calling split, 
c       extract the nth field from the current record
c               and convert it to real
c       takes a string and converts it to floating point number
c          ..this goes with SPLIT routine, which is jrw's
c          fortran hack to deal with field free format

c       modified 8/8/00 by jrw
c               to include input as integer returned as float
        
        n=nl(nfld)-nf(nfld)+1

        mult=1.0
        zneg=1.0
        nd=0
        ne=0
        e=1.0

c       check for negative sign
        if(rec(nf(nfld):nf(nfld)).eq.'-')zneg=-1.0
        
c       find decimal point and exponent
        
        do 300 i=1,n
        c=rec(nf(nfld)+i-1:nf(nfld)+i-1)
        if(c.eq.'.')nd=i
        if(c.eq.'E')ne=i
        if(c.eq.'e')ne=i
300     continue

c       ***modified 8/8/00  for the case of integer
        if(nd.eq.0)nd=n+1
        
c       handle left of decimal

        f=0.0
        do 302 i=nd-1,1,-1
        c=rec(nf(nfld)+i-1:nf(nfld)+i-1)
        if(c.eq.'-')goto 302
        k=ifc(c)
        f=f+mult*float(k)
        mult=mult*10.0
302     continue

c       handle exponent

        if(ne.eq.0)then
                nt=n
        else
                nt=ne-1
                nmult=1
                j=0
                ineg=0
                do 310 i=ne+3,ne+1,-1
                c=rec(nf(nfld)+i-1:nf(nfld)+i-1)
                if(c.eq.'-')then
                        ineg=1
                        goto 310
                endif
                if(c.eq.'+')then
                        ineg=0
                        goto 310
                endif
                k=ifc(c)
                j=j+k*nmult
                nmult=nmult*10
310             continue
                if(ineg.eq.1)j=-j
                e=10.0**j
        endif

c       handle right of decimal

        mult=0.1
        do 301 i=nd+1,nt
        c=rec(nf(nfld)+i-1:nf(nfld)+i-1)
        k=ifc(c)
        f=f+mult*float(k)
        mult=mult/10.0
301     continue

        f_ex=f*e*zneg

        return
        end

        integer function i_ex(nfld)

        include 'isodist.h'

        character*1 c

c       integer extraction
c       extract the nth field from the current record
c               and convert it to an integer
c               used in conjunction with SPLIT routine to implement
c               field free input
        
        n=nl(nfld)-nf(nfld)+1

        mult=1
        j=0
        ineg=0

        do 300 i=n,1,-1
        k=-1

c	JRW 4/2/17        
c       c=rec(nf(nfld)+i-1:nf(nfld)+1-1)
        c=rec(nf(nfld)+i-1:nf(nfld)+i-1)

        
        if(c.eq.'-')then
                ineg=1
                goto 300
        endif
        
        k=ifc(c)                
        j=j+k*mult
        mult=mult*10
        
300     continue
        if(ineg.eq.1)j=-j
        
        i_ex=j
        
        return
        end
        integer function ifc(c)
        character*1 c

c       "10 for efficiency, 2 for elegance"
c           not my proudest hour, but got fed up with
c           stupid fortran nonsense trying to read simple files
c               that any scripting language handles with ease
c               resulting in SPLIT routine, and amazingly
c               insightful pieces of code like the following:


        if(c.eq.'0')k=0
        if(c.eq.'1')k=1
        if(c.eq.'2')k=2
        if(c.eq.'3')k=3
        if(c.eq.'4')k=4
        if(c.eq.'5')k=5
        if(c.eq.'6')k=6
        if(c.eq.'7')k=7
        if(c.eq.'8')k=8
        if(c.eq.'9')k=9

        ifc=k
        return
        end
        integer function lengthn(n,str)
        character*1 str(n)

c       returns length of a string, ignoring trailing blanks
c               why would you need to know that?

        do 300 i=n,1,-1
        if(str(i).ne.' ')then
                lengthn=i
                return
        endif
300     continue 
        return
        end
        subroutine split(c)

c       generalized split function
c
c       takes contents of "rec" and split at each instance of character "c"
c       used for field free input with c=" ", or for string handling
c       with other chars, such as ",", ".", "_", "/"
c               ...written as a crude fortran hack by jrw
c               however, it solves a lot of problems
c
c       nstr gives the number of strings
c       arrays nf and nl give the position of the first and last chars
c               of each string

        include 'isodist.h'
        
        character*1 c
        
        do 300 i=1,20
        nf(i)=0
        nl(i)=0
300     continue

        nstr=0
        istr=0
        ipflag=0
        
c       find last nonblank character

        do 302 i=200,1,-1
        call char_num(rec(i:i),j)
        if(j.gt.0)then
                lastchar=i
                goto 800
        endif
302     continue

800     continue

c       char(9) is the tab character

        do 301 i=1,lastchar
        
        if(rec(i:i).eq.'(')ipflag=1
        if(rec(i:i).eq.')')ipflag=0
        if(istr.eq.0)then
                if(rec(i:i).ne.c.and.rec(i:i).ne.char(9))then
                        nstr=nstr+1
                        nf(nstr)=i
                        istr=1
                        go to 301
                endif
        else
                if(rec(i:i).eq.c.or.rec(i:i).eq.char(9))then
                        if(ipflag.eq.1)goto 301
                        nl(nstr)=i-1
                        istr=0
                        goto 301
                endif
        endif
301     continue
        if(nl(nstr).eq.0) nl(nstr)=lastchar

        return
        end
