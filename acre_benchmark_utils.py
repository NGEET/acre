import code   # For development: code.interact(local=locals())
import numpy as np
import xml.etree.ElementTree as et
from scipy.io import netcdf
import matplotlib.dates as mdates
from datetime import datetime, timedelta


# BA_SCPF    (SIZE x PFT)
# DDBH_SCPF  (SIZE x PFT)
# (M1_SCPF + M2_SCPF + M3_SCPF + M4_SCPF + M5_SCPF + M6_SCPF + M7_SCPF + M8_SCPF) / NPLANT_SCPF
# RECRUITMENT (PFT)


# This object is bound to each site
# It should contain a list of viable benchmarks

class benchmark_obj:

    def __init__(self,census_filename):
        
        self.bvarlist = []

        # Lets check through the census file and see if any of these variables
        # are in the file.  We will later look through the model output
        # and pop off list entries that are not there.

        if(census_filename.strip() != ''):
            
            print(census_filename)
            
            fp = netcdf.netcdf_file(census_filename, 'r', mmap=False)

            cens_var_name = 'basal_area_by_size_census'
            if (fp.variables.has_key(cens_var_name)):
                self.bvarlist.append(benchmark_vars('Basal Area','BA_SCPF',cens_var_name,'scpf','size-class','m2/ha'))
                self.bvarlist[-1].load_census(fp)
            else:
                print('Census variable: '+cens_var_name+', was not found in the census file')

            cens_var_name = 'growth_increment_by_size_census'
            if (fp.variables.has_key(cens_var_name)):
                self.bvarlist.append(benchmark_vars('Growth Increment','DDBH_SCPF',cens_var_name,'scpf','size-class','cm/yr'))
                self.bvarlist[-1].load_census(fp)
            else:
                print('Census variable: '+cens_var_name+', was not found in the census file')




# =======================================================================================

class benchmark_vars:

    def __init__(self,name,mod_symbol,obs_symbol,mod_dimclass,obs_dimclass,unit):

        self.mod_symbol = mod_symbol
        self.obs_symbol = obs_symbol
        self.name   = name
        self.scv    = False    # Size-class variable
        self.sv     = False    # Scalar variable
        self.unit   = unit
        
        self.hfile_id = -9
        

        self.mod_dimclass = mod_dimclass
        self.obs_dimclass = obs_dimclass

        # This is a size-dimensioned variable
        if(self.mod_dimclass == 'scpf'):

            # Initialize size-class array
            self.scv = True
#            self.scv_mod_ar = np.zeros((hdims.nsize,n_htypes))
#            self.scv_mod_n  = np.zeros((hdims.nsize,n_htypes))
            
#            self.scv_x      = np.zeros((hdims.nsize,1))
            self.scv_x_unit = 'DBH [cm]'

            # Initialize scalar array
            self.sv        = True
#            self.sv_mod_ar = np.zeros((n_htypes,1))
#            self.sv_mod_n  = np.zeros((n_htypes,1))
            self.sv_obs_ar = np.zeros((1,1))


        if(self.mod_dimclass == 'scalar'):

            # Initialize scalar array
            self.sv        = True
#            self.sv_mod_ar = np.zeros((n_htypes,1))
#            self.sv_mod_n  = np.zeros((n_htypes,1))


    def load_census(self,fp):

        code.interact(local=locals())

        # The 'size-class' type census variable is actually a
        # fp.variables[self.obs_symbol].shape
        # (7, 10, 3)
        # fp.variables[self.obs_symbol].dimensions
        # ('cens', 'dclass', 'confidence')

        fp.variables[self.obs_symbol].data

        d_sizes = fp.variables[self.obs_symbol].shape
        dim_names = fp.variables[self.obs_symbol].shape

        if(self.obs_dimclass == 'size-class'):

            if(dim_names[0] != 'cens'):
                print('expected census data to have cens as first dimension')
                print('exiting')
                exit(2)

            # Condense the census dimension into 1 size
            self.scv_obs_ar = np.zeros((d_sizes[1],d_sizes[2]))
            self.scv_obs_ar[:,:] = np.mean(fp.variables[self.obs_symbol].data, axis=0 )
        
        else:
            print('The census variable: {},   with dim-type: {}'.(self.obs_symbol,self.obs_dimclass))
            print(' does not have dimensions set-up yet. Exiting')
            exit(2)
