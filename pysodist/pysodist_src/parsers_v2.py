import numpy as np
from numpy.fft import fft,ifft,rfft,irfft
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from math import pi,e
import ast

class MainInfoParser(): #parses the .in file to hold relevant variables
    def __init__(self, in_file):
        self.file=in_file
        f=open(in_file,'r')
        self.opt = f.readline().split()[0]

        self.batchfile = f.readline().split()[0]
        self.atomfile = f.readline().split()[0]
        self.resfile = f.readline().split()[0]
        self.fit_iter = int(f.readline().split()[0])
        self.sig_global = float(f.readline().split()[0])

        baseline_line = f.readline()
        self.b_init = baseline_line.split()[0]
        self.b_tag = baseline_line.split()[1]

        self.m_off_init = float(f.readline().split()[0])
        self.gw_init = float(f.readline().split()[0])

        f.close()

class AtomInfoParser(): #parses the atom model file to hold relevant data - masses and frequencies
    def __init__(self, file='atom_defs.txt'):
        self.file = file
        self.atom_masses = dict()
        self.atom_freqs = dict()

        reading=True
        current_atom=None
        f = open(self.file,'r')
        while reading:
            new_line=f.readline()
            if new_line=='':
                reading=False
                break
            elif new_line.split()[1][0].isalpha(): #check if new atom
                current_atom=new_line.split()[1]
                self.atom_masses[current_atom]=[]
                self.atom_freqs[current_atom]=[]
            else: #otherwise add mass,freq to data
                split_line=new_line.split()
                self.atom_masses[current_atom].append(float(split_line[0]))
                self.atom_freqs[current_atom].append(float(split_line[1]))

class ResInfoParser(): #parses the residue model file. Gives the variable residues, the atomic composition of each residue, as well as amplitudes and variable atoms
    def __init__(self, res_file, atom_parser, corr_atoms = None):
        f = open(res_file,'r')
        f.readline()#skip initial line
        ##################
        #collect all info
        ##################

        #Parse the species data - number of species and initial amplitudes
        self.num_species = int(f.readline().split()[0]) #number of species
        self.species_names=[] #U, L, F, etc. Name of the species (labeling condition)
        self.species_amps= [] #Amplitudes of each species
        for i in range(self.num_species):
            line=f.readline().split()
            self.species_names.append(line[0])
            self.species_amps.append(float(line[1]))

        #Parse the variable atoms and their initial values
        num_atoms = int(f.readline().split()[0]) #number of atom types in model: C, H, N, O, S, and variable atoms, etc.
        self.atom_names = [] #Elemental symbol
        self.atom_init_values = [] #Initial frequency of special atoms (variable or fixed)
        self.atom_modes = [] #How to treat the special atoms
        self.corr_atoms = corr_atoms

        if corr_atoms is not None:
            self.residue_corrs = dict()#[] #Correlated carbons for each amino acid
            self.max_corr_group = 0
            calc_corr = True
        else:
            calc_corr = False

        for i in range(num_atoms):
            line=f.readline().split()
            self.atom_names.append(line[0])
            self.atom_modes.append(line[1])
            try:
                self.atom_init_values.append(float(line[2]))
            except:
                self.atom_init_values.append(None)

        #residues multiplicity of atoms
        #keys are symbols, values lists of atom multiplicities, one for each species
        self.residue_composition=dict()
        reading=True
        self.residue_names = [] #One letter identifier for amino acids
        self.residue_modes = [] #Modes for each amino acid, for SILAC
        self.residue_init_values = [] #Initial frequencies for special residues

        self.res_freqs = dict() #Dictionary for frequencies overall. Frequency is just 1
        #for most residues, but are set to be the values in residue_init_values in
        #pysodist

        while reading:
            line = f.readline()
            if line=='':
                reading=False
                break

            current_res = line.split()[0]
            current_mode = line.split()[1]
            if calc_corr:
                try:
                    if line.split()[-1] == "[]":
                        self.residue_corrs[current_res] = []
                    else:
                        self.residue_corrs[current_res] = ast.literal_eval(line.split()[-1])

                    if max(self.residue_corrs[current_res], default = 0) > self.max_corr_group:
                        self.max_corr_group = max(self.residue_corrs[current_res])
                except:
                    print("Correlation definitions incorrect")
            self.residue_names.append(current_res)
            self.residue_modes.append(current_mode)


            init_values = []
            for freq in line.split()[2:]:
                try:
                    init_values.append(float(freq))
                except:
                    pass

            self.residue_init_values.append(init_values)
            self.residue_composition[current_res] = []
            self.res_freqs[current_res] = []

            for i in range(self.num_species):
                a = f.readline()
                line=list(map(int,a.split()[:num_atoms]))
                self.residue_composition[current_res].append(line)
                self.res_freqs[current_res].append(1)

        #If we use correlated atoms, we change the "residue composition" to reflect this
        #using artificial "atoms" that are doublets, triplets, etc.
        if calc_corr:
            corr_atom_names = []
            corr_atom_modes = []
            corr_atom_init_values = []
            for k in range(num_atoms):
                atom = self.atom_names[k]
                corr_atom_names.append(atom)
                corr_atom_modes.append(self.atom_modes[k])
                corr_atom_init_values.append(self.atom_init_values[k])
                if atom == self.corr_atoms:
                    for num in range(2,self.max_corr_group+1):
                        corr_atom_names.append(atom + str(num))
                        corr_atom_modes.append("Correlated")
                        corr_atom_init_values.append(None)


            for i in range(self.num_species):
                for residue in self.residue_composition:
                    correlation_sets = list(set(self.residue_corrs[residue]))
                    correlation_mults = list(map(lambda n:self.residue_corrs[residue].count(n), correlation_sets))
                    zero_padded_mults = [] #Zero padded mults is a vector [#2, #3, #4, #5], symbolizing
                    #the number of doublets to quintets

                    for k in range(2,self.max_corr_group+1):
                        if k in correlation_sets:
                            zero_padded_mults.append(correlation_mults[correlation_sets.index(k)])
                        else:
                            zero_padded_mults.append(0)

                    corr_residue_composition = []

                    for k in range(num_atoms):
                        atom = self.atom_names[k]
                        if atom == self.corr_atoms:
                            #print(self.corr_atoms)
                            if self.residue_composition[residue][i][k] > 0: #To determine if this is a species containing the variable atom
                                corr_residue_composition.append(self.residue_composition[residue][i][k] - sum(self.residue_corrs[residue]))
                                corr_residue_composition += zero_padded_mults
                            else:
                                corr_residue_composition += (len(zero_padded_mults)+1)*[0]
                        else:
                            corr_residue_composition.append(self.residue_composition[residue][i][k])
                    self.residue_composition[residue][i] = corr_residue_composition
                    print(self.atom_names)
            self.atom_names = corr_atom_names
            self.atom_modes = corr_atom_modes
            self.atom_init_values = corr_atom_init_values
        self.original_species_amps = self.species_amps.copy()


class BatchInfoParser(): #Parses the batch file in the input file. Gives the sequence, charges, and spectrum file
    def __init__(self, batch_file):
        f=open(batch_file,'r')
        #f.readline() #skip header line
        batch_syms = [] #List of symbols for each peptide. Used as index for multiplicities
        batch_mults = [] #Multiplicities for number of each residue in a peptide
        charges=[] #Charges for each peptide
        data_files=[] #File path for each peptide's spectrum file
        pep_names=[] #Sequence of the peptide on paper, without added protons/water
        reading=True
        while reading:
            line = f.readline()
            if line=='':
                break
            line=line.split()
            pepseq=line[0]

            charge=int(line[1])
            charges.append(charge)
            data_files.append(line[2]) #skip over retentiontime
            #build batch_model
            pep_names.append(pepseq)
            pepseq+='Z'+charge*'X'#adds on a water and protons

            res_syms = list(set(pepseq))
            batch_syms.append(res_syms)
            res_mults = list(map(lambda sym:str.count(pepseq,sym),res_syms))
            batch_mults.append(res_mults)

        ########################
        #Externally useful data
        ########################

        self.num_batches=len(batch_mults)
        self.charges=charges
        self.data_files=data_files
        self.pep_names=pep_names
        self.batch_syms = batch_syms
        self.batch_mults = batch_mults


class ExpSpectrumParser(): #parses the experimental spectrum data file
    def __init__(self,data_file,charge):
        lines=open(data_file,'r').read().split('\n')[:-1] #cut off empty line at end
        self.masses=[]
        self.raw_intensities=[]
        try:
            for line in lines:
                splitline=line.split(' ')
                self.masses.append(float(splitline[0])*charge)
                self.raw_intensities.append(float(splitline[1]))
            self.m_hd=self.masses[0]
            self.largest_mass=self.masses[-1]
        except:
            self.m_hd=None
    def get_unbinned_target_array(self):#returns a tuple of target masses, values
        #subtract off baseline offset
        min_intensity=min(self.raw_intensities)
        baselined_intensities=list(map(lambda x:x-min_intensity,self.raw_intensities))
        self.vert_shift=min_intensity
        return self.masses,baselined_intensities

    def get_target_array_v1(self,N,dm,exp_box_size,m_hd):
        #subtract off baseline offset
        min_intensity=min(self.raw_intensities)
        baselined_intensities=list(map(lambda x:x-min_intensity,self.raw_intensities))
        self.vert_shift=min_intensity
        #bin into N,dm array
        target_array=np.zeros(N,'float')
        for i,mass in enumerate(self.masses):
            index=round(exp_box_size/dm*round((mass-m_hd)/exp_box_size))%N
            #index=round((mass-m_hd)/dm)%N
            target_array[index]+=baselined_intensities[i]#normalized_intensities[i]
        return self.cubic_interpolate_zeros(target_array)
    def cubic_interpolate_zeros(self,array):
        non_zero_indices=[]
        non_zero_values=[]
        for i,val in enumerate(array):
            if val>0:
                non_zero_indices.append(i)
                non_zero_values.append(val)
        f_interp=interp1d(non_zero_indices,non_zero_values,kind='cubic')
        out=array.copy()
        for i in range(len(out)):
            if out[i]==0:
                try:
                    out[i]=max(f_interp(i),0)
                except: #i out of interpolation range, just leave
                    pass
        return out
