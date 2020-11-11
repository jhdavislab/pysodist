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

c       global parameters

        parameter (np=65536)            ! np is the number of points calculated
        parameter (ncp=np/2)            ! ncp is the number of complex points                           
        parameter (natoms=24)           ! natoms is the maximum number of atom types defined                            
        parameter (niso_max=4)          ! niso_max is the maximum number of isotopes defined                            
        parameter (nelements=103)       ! nelements is the maximum number of elements                           
        parameter (scale_mz=1000.0)     ! scale_mz determines the number of points calculated per dalton                                
        parameter (nres_max=100)        ! nres_max determines the max peptide length                            
        parameter (n_species_max=10)    ! n_species_max defines the max total number of species                         
        parameter (n_var_atom_max=20)   ! n_var_atom_max  defines the max number of variable atom parameters                            
        parameter (n_var_res_max=20)    ! n_var_res_max defines the max number of variable residue parameters                           
        parameter (n_res_lib=60)        ! n_res_lib defines the max number of residues in the residue library                           
        parameter (ndim=16)             ! ndim defines the max number of fitting parameters: physical array dimension                           
        parameter (n_round_max=10)      ! n_round_max defines the maximum length of the fitting schedule
        parameter (nwork=np)            ! nwork defines the size of the working array for minpack
        parameter (nwork2=4*np+15)      ! nwork2 defines the size of the working arrays for fftpack5
        
        integer niso(nelements)
        real miso(nelements,niso_max)
        real fiso(nelements,niso_max)
        common /iso/ niso,miso,fiso

        character*2 atom_symb(nelements)
        integer atom_num(nelements),atom_fix(nelements)
        character*40 atom_comment(nelements)
        common /periodic/ atom_num,atom_fix,atom_comment,atom_symb

c       arrays for ft of element
        real el_ft(np,nelements)
        complex cel_ft(ncp,nelements)
        equivalence(el_ft,cel_ft)
        common /elements_ft/ el_ft
        
c       arrays for various spectra      
        real eft(np),eft_obs(np),eftd(np)
        complex ceft(ncp),ceft_obs(ncp)
        equivalence (eft,ceft),(eft_obs,ceft_obs)
        common /spectra/ eft,eft_obs,eftd
        
c       arrays for species spectra
        real seft(np,n_species_max)
        complex cseft(ncp,n_species_max)
        equivalence (seft,cseft)
        common /species/ seft

c       arrays for residue spectra
        real eft_fix(np,n_res_lib,n_species_max)
        complex ceft_fix(ncp,n_res_lib,n_species_max)
        real eft_var(np,n_res_lib,n_species_max)
        complex ceft_var(ncp,n_res_lib,n_species_max)
        real eft_res(np,n_res_lib,n_species_max)
        complex ceft_res(ncp,n_res_lib,n_species_max)
        equivalence(eft_fix,ceft_fix)
        equivalence(eft_var,ceft_var)
        equivalence(eft_res,ceft_res)
        common /res_spectra/ ceft_fix,ceft_var,ceft_res

c       arrays for intermediate calc. functions
        real pft(np),gft(np)
        complex cpft(ncp),cgft(ncp)
        equivalence (pft,cpft),(gft,cgft)
        common /intermed/pft,gft

        real pi,pi2,gw3
        common /constants/ pi,pi2,gw3

        real mz(np)
        real cmz(ncp)
c	JRW 4/2/17
c       real mz_ptr(np)
	integer mz_ptr(np)
        real mz_hd
        real sig_global
        integer mzfr,mzto,nmz,nz,pep_length
        character*100 pepseq
        character*1 pepseq_res(100)
        equivalence(pepseq,pepseq_res)
        integer mol_form(natoms)
        common /ms/ mz,cmz,mz_ptr,mzfr,mzto,nmz
        common /ms2/ nz,mz_hd,sig_global,pepseq,pep_length
        common /ms3/ mol_form
        
        real xobs(np),yobs(np),sig(np),ysort(np),ysort2(np)
        common /obsdata/xobs,yobs,sig,ysort,ysort2
        
        character*5 opt
        character*200 infile,outfile,datfile,bogusfile,batchfile,batchout,froot
        character*70 atomfile,resfile
        character*200 outlabels
        character*6 parlabel(ndim)
        character*4 fix,free
        character*4 fixfree(ndim)
        common /iofiles/ infile,outfile,datfile,opt,fix,free,parlabel,fixfree
        common /iofiles2/ atomfile,resfile
        common /bogus/ bogusfile,batchout,outlabels,froot
        
        character*2 protlab,mfslug,species_lab(n_species_max)
        character*3 ptim
        common /misc_strings/ protlab,mfslug,ptim, species_lab
        
        real pmz,prt,frac_lab,frac_lab_err,chisq_end
        integer psymb,pep_len
        common /misc_vars/ pmz,prt,psymb,frac_lab,frac_lab_err,chisq_end,pep_len
        
        integer fit_schedule(n_round_max,ndim)
        real mw,b,gw,ftol,off,ul_amp,x(ndim),x_init(ndim)
        real frac1,frac2
        integer nround,nfitpar,nround_ctr
        integer fixed_pars(ndim)
        common /fitparams/ mw,b,gw,ftol,off,ul_amp,x,x_init
        common /fitparams2/ nround,nfitpar,nround_ctr,fit_schedule,frac1,frac2,fixed_pars
        
        character*1 residue(nres_max)
        common /res/ residue
        
        integer natom_types,n_var_atom
        integer res_atm_comp(nres_max,natoms,n_species_max)
        integer atom_ptr_list(natoms)
        common /reslib/ natom_types, n_var_atom, res_atm_comp, atom_ptr_list

        integer fixed_var_res(nres_max)
        integer fixed_var_atom(natoms,nres_max,n_species_max)
        integer n_var_res
        integer var_atom(natoms),var_res(n_res_lib)
        integer active_res(n_res_lib)
        integer fix_res(n_var_res_max)
        integer n_fix_res
        common /reslib2/ fixed_var_res, fixed_var_atom, n_var_res, var_atom, var_res, active_res
        common /reslib4/ fix_res,n_fix_res
        
        integer amp_idx(n_species_max),frac_idx(n_var_atom_max),phi_idx(n_var_res_max)
        real amp(n_species_max),frac(n_var_atom_max),phi(n_var_res_max), phi_fix(n_var_res_max)
        common/reslib3/ amp_idx,frac_idx,phi_idx,amp,frac,phi,phi_fix

        integer err_flag,auto_base_flag,nctr
        common /flags/ err_flag,auto_base_flag,nctr

        integer n_species,n_var_atoms,n_var_residues,var_res_idx,var_res_num
        integer species_array(n_species_max),var_atom_array(n_var_atom_max),var_res_array(n_var_res_max)
        common /species/ n_species,species_array,n_var_atoms,n_var_residues,var_atom_array,var_res_array
        common /species2/ var_res_idx,var_res_num

        character*200 rec
        integer nf(20),nl(20),nstr
        common /parser/ nf,nl,nstr,rec
        
        real wsave(nwork2),work(nwork2)
        common /work/ wsave,work
