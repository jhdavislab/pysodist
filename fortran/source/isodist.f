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

        program isodist
                
        include 'isodist.h'
        
        real xx(ndim),x0(ndim),xdel(ndim)
        external res_min

C       Write out license information
          write(6,*)
          write(6,600)
600   format('isodist Copyright (C) 2008 James R. Williamson ')
          write(6,601)
601   format('and Michael T. Sykes\n')
          write(6,602)
602   format('This program comes with ABSOLUTELY NO WARRANTY\n')
          write(6,603)
603   format('This is free software, and you are welcome to ')
          write(6,604)
604   format('redistribute it under the terms of the')
          write(6,605)
605   format('GNU General Public License\n')
      write(6,606)
606   format('isodist_release v1.0 - May 14th, 2008\n')
        
c       read in control data

        do 355 i=1,ndim
        fixed_pars(i)=1
355     continue

        call getarg(1,infile)
        open(unit=1,file=infile,form="formatted",status="old")

        read(1,112) rec
        call split(' ')
        opt=rec(nf(1):nl(1))

        read(1,112) rec
        call split(' ')
        batchfile=rec(nf(1):nl(1))

        read(1,112) rec
        call split(' ')
        atomfile=rec(nf(1):nl(1))

        read(1,112) rec
        call split(' ')
        resfile=rec(nf(1):nl(1))

        read(1,112) rec
        call split(' ')
        niter=i_ex(1)

        read(1,112) rec
        call split(' ')
        sig_global=f_ex(1)

        read(1,112) rec
        call split(' ')
        x_init(1)=f_ex(1)
        if(rec(nf(2):nl(2)).eq."fixed")fixed_pars(1)=0
        if(rec(nf(2):nl(2)).eq."auto")then
                fixed_pars(1)=0
                auto_base_flag=1
        endif

        read(1,112) rec
        call split(' ')
        x_init(2)=f_ex(1)
        if(rec(nf(2):nl(2)).eq."fixed")fixed_pars(2)=0

        read(1,112) rec
        call split(' ')
        x_init(3)=f_ex(1)
        if(rec(nf(2):nl(2)).eq."fixed")fixed_pars(3)=0
                
        close(unit=1)

112     format(a80)
113     format(a200)
c       set up input and output files

        open(unit=17,file=batchfile,status="old",form="formatted")
        write(6,*)batchfile

        rec=batchfile
        call split('.')
        batchout=rec(1:nl(nstr))//".csv"
        open(unit=18,file=batchout,status="unknown",form="formatted")
        write(6,*)batchout
        
c       calculate mz axis
        do 300 i=1,np
        mz(i)=float(i)/scale_mz
300     continue

c       calculate mu axis
        do 301 i=1,ncp
        cmz(i)=float(i-1)/(float(np))
301     continue

c       initialize error array  --- not used
        do 360 i=1,np
        sig(i)=sig_global
360     continue

c       initialize working array for fftpack5
        call bogus_ft(eft,np,0)
        
c       setup atoms and calculate atomic mu-domain functions
        call setup_atoms
        
c       set up residue library
        call setup_residues
        
c       set up species and labeling model
        call setup_model

c       precalculate fixed part of residues
        call calc_residues(0)

c       misc constants
        mfslug="XX"
        psymb=19
        fix="FIX "
        free="FREE"
        npep_ctr=0

c       MAIN LOOP

900     continue

c       read line from batch file and extract filename, charge, and sequence

        read(17,113,end=901)rec
        call split(' ')
        nz=i_ex(2)
        infile=rec(nf(3):nl(3))

c       add terminal residues and residues for charge
c               NEED TO MODIFY THIS SO TERMINAL RESIDUES INCLUDE CHARGE 

        pepseq=rec(nf(1):nl(1))//"Z"
        pep_length=nl(1)-nf(1)+2
        pep_len=nl(1)-nf(1)+1
        
        npep_ctr=npep_ctr+1
        
        do 340 i=1,nz
        pepseq=pepseq(1:pep_length)//"X"
        pep_length=pep_length+1
340     continue

        write(6,*)
        write(6,665)npep_ctr,pepseq, nz, pep_length
665     format("PEPTIDE # ",i4,1x,a,1x,i2,1x,i2)

        rec=infile
        call split('.')
        outfile=rec(1:nl(nstr-1))//".fit"
        datfile=rec(1:nl(nstr-1))//".dat"

c       extract info from filename:  timepoint, protein #, mz, retention time
        
        rec=outfile
        call split('/')
        froot=rec(nf(nstr):nl(nstr))

        rec=froot
        call split('_')
        protlab=rec(nf(1)+1:nl(1))

        lc=nf(2)
        do 341 i=nf(2),nl(2)
        call char_num(froot(i:i),ic)
        if(ic.eq.2)lc=i
341     continue

        ptim=froot(nf(2):lc)
        pmz=f_ex(3)
        prt=f_ex(4)
        
c       set active_residue flags
        
        do 320 i=1,n_res_lib
        active_res(i)=0
320     continue

        var_res_num=0
        do 321 i=1,pep_length
        do 323 j=1,n_res_lib
        if(residue(j).eq.pepseq(i:i)) active_res(j)=1
c       if(j.eq.var_res_idx) var_res_num=var_res_num+1
323     continue
321     continue

        call pep_mw

c       read in peak to be fit
        call read_dino
        
        if(auto_base_flag.eq.1)call auto_baseline

c       write out experimental peak file   DO WE REALLY NEED THIS?

        open(unit=4,file=datfile,status="unknown",form="formatted")
        write(4,611)(mz(mz_ptr(i))+mz_hd,eft_obs(mz_ptr(i)),i=1,nmz)
611     format(e16.8,",",e16.8)
        close(unit=4)
        
        idone=0
        nit=2

c       set parameters to initial values
        
        do 375 j=1,nfitpar
        x(j)=x_init(j)
375     continue

        if(opt.eq."tryit")goto 990

c       do minimization 

        call marq_lmdif1(niter,res_min,idone)

374     continue

990     continue

c       evaluate with final fit params
        call calc_iso_dist

        chisq=0.0
        do 388 i=1,nmz
        ii=mz_ptr(i)
        chisq=chisq+(eft_obs(ii)-eft(ii))*(eft_obs(ii)-eft(ii))
388     continue
        chisq_end=chisq

c       write out the fit
        call plotxy_range(eft,np,mzfr,mzto)
c       call plotxy_range(eft,np,1,np)

c       output final results to batch csv file
        call write_model
        
        goto 900

901     continue
        
        end
        
        
