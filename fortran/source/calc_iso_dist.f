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

        subroutine calc_iso_dist
                
        include 'isodist.h'
        
        complex efac,gfac,cgw

c       calculate the isotope distribution in the u-domain

        pi=4.0*atan(1.0)
        pi2=2.0*pi

        call model_params

c       calculate the offset: h(u)

        do 303 j=1,ncp
        efac=complex(0.0,cmz(j)*pi2*((off+mz_hd)*scale_mz))
        cpft(j)=cexp(-efac)
303     continue

c       calculate the gaussian: g(u)

        cgw=complex(gw,0.0)*sqrt(2.0)
        gw3=gw*gw*gw
        gw2p=gw*sqrt(pi2)

        do 304 j=1,ncp
        gfac=(complex(cmz(j),0.0)/(cgw))**2
        cgft(j)=cexp(-gfac)/gw2p
304     continue

        call calc_atoms
        call calc_residues(1)
        call calc_species
        call calc_model

        return
        end
