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

        subroutine model_params

c       the minimization routine uses the variable x to hold the fitting parameters.
c       Because of the dynamic assignment of fitting parameters based on the
c       input model, we have to keep track of the variables in the arrays
c       amp, frac, and  phi, and pack them into x for each call.  The
c       indices amp_idx, frac_idx, and phi_idx are set up in setup_model.f

        include 'isodist.h'

        b=x(1)
        off=x(2)
        gw=x(3)

        do 300 i=1,n_species
        amp(i)=x(amp_idx(i))
300     continue

        do 301 i=1,n_var_atom
        frac(i)=x(frac_idx(i))
301     continue

        do 302 i=1,n_var_res
        phi(i)=x(phi_idx(i))
302     continue

        return
        end
