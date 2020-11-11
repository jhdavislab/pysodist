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

        subroutine marq_lmdif1(niter,res_min,idone)
                
        include 'isodist.h'

        real xadj(ndim)
        real wa(nwork)
        integer iwa(ndim)
        real ycalc(np)
        
c       driver routine for marquardt minimization

        chisq=0.0
        nctr=0
        
        do 374 nj=1,nround
        nctr=nctr+1

c       fix or free parameters according to fit schedule

        nadj=0
        do 375 k=1,nfitpar
        if(fit_schedule(nj,k).eq.0)fixfree(k)=fix
        if(fit_schedule(nj,k).eq.1)then
                fixfree(k)=free
                nadj=nadj+1
                xadj(nadj)=x(k)
        endif
375     continue

        write(6,660)nj,(fixfree(k),parlabel(k),k=1,nfitpar)
660     format('ROUND ',i2,' CHISQ    ',6(a4,1x,a6,1x))


        tol=sqrt(spmpar(1))
        tol=1.0e-7
        if(nj.eq.nround)tol=1.0e-7
        lwa=ndim*np
        
        call lmdif1(res_min,nmz,nadj,xadj,ycalc,tol,info,iwa,wa,nwork)

374     continue

c       FINISH UP
        
        write(6,*)
        write(6,683)
683     format('FINAL FIT PARAMETERS')
        write(6,*)
        write(6,682)(parlabel(i),x(i),i=1,nfitpar)
682     format(a6,1x,f14.6)
        write(6,*)

        chisq_end=chisq

992     continue
        return
        end
        
