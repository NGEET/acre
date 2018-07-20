import code   # For development: code.interact(local=locals())
import numpy as np
import xml.etree.ElementTree as et
from scipy.io import netcdf
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import matplotlib.pyplot as plt


# BA_SCPF    (SIZE x PFT)
# DDBH_SCPF  (SIZE x PFT)
# (M1_SCPF + M2_SCPF + M3_SCPF + M4_SCPF + M5_SCPF + M6_SCPF + M7_SCPF + M8_SCPF) / NPLANT_SCPF
# RECRUITMENT (PFT)


# This object is bound to each site
# It should contain a list of viable benchmarks

class benchmark_obj:

    def __init__(self,census_filename):
        
        self.bvarlist = []
        self.census_filename = census_filename

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

            fp.close()

    # ===================================================================================
    ## Check the first history file in the list to see which benchmarking variables
    #  are available.

    def init_history(self,hist_file0,n_htypes):


        fp = netcdf.netcdf_file(hist_file0, 'r', mmap=False)

        for bvar in self.bvarlist:
            
            if( fp.variables.has_key(bvar.mod_symbol) ):

                bvar.active = True

                if( bvar.mod_dimclass == 'scpf'):
                    
                    dims = fp.variables[bvar.mod_symbol].dimensions
                    if(dims[1] != 'fates_levscpf'):
                        print('An SCPF benchmark variable: {} does not actually'.format(bvar.mod_symbol))
                        print(' have the correct dimensions: {}... exiting'.format(dims))
                        exit(2)

                    fates_levscls = fp.variables['fates_levscls'].data

                    if (fates_levscls[0] == 0.0):
                        bvar.offset0 = True

                        if(fates_levscls.size-1 != bvar.scv_obs_ar.shape[0]):
                            print('Dimensions of model output size-classes dont match observations')
                        for isc,scvar in enumerate(fates_levscls[1:]):
                            if( np.abs(scvar- bvar.scv_x[isc])>1.0e-10 ):
                                print('Dimensions of model output size-classes dont match observations')
                                print('Observed classes: {}'.format(bvar.scv_x))
                                print('Modeled (0 is ignored): {}',format(fates_levscls))

                    else:
                        bvar.offset0 = False
                        if(fates_levscls.size != bvar.scv_obs_ar.shape[0]):
                            print('Dimensions of model output size-classes dont match observations')
                        for isc,scvar in enumerate(fates_levscls[:]):
                            if( np.abs(scvar- bvar.scv_x[isc])>1.0e-10 ):
                                print('Dimensions of model output size-classes dont match observations')
                                print('Observed classes: {}'.format(bvar.scv_x))
                                print('Modeled (0 is ignored): {}',format(fates_levscls))


                    d_sizes = bvar.scv_obs_ar.shape
                    bvar.modlist = []
                    for imod in range(n_htypes):
                        bvar.modlist.append(mod_scv_array(d_sizes[0]))

        fp.close()
    
    # ===================================================================================

    def load_history(self,filename,h_index,site_index):

        # Objective is to push new estimates of the benchmark variables 

        fp = netcdf.netcdf_file(filename, 'r', mmap=False)

        for bvar in self.bvarlist:

            if(bvar.active):

                if( (bvar.obs_dimclass=='size-class') and (bvar.mod_dimclass=='scpf') ):

                    # Create a mapping between FATES size-classes and
                    fates_levscls       = fp.variables['fates_levscls'].data
                    fates_scmap_levscpf = fp.variables['fates_scmap_levscpf'].data
                    
                    d_sizes = fp.variables[bvar.mod_symbol].shape
                    # fates_scmap_levscpf
                    # These are the expected dimensions
                    # ('time', 'fates_levscpf', 'lndgrid')

                    for itime in range(d_sizes[0]):
                        for isc,isc_map0 in enumerate(fates_scmap_levscpf):
                            
                            if ( (bvar.offset0 == False) or ((bvar.offset0 == True)and(isc_map0-1 != 0 ))   ):

                                isc_map = isc_map0-1

                                if(bvar.offset0 == False): arr_map = isc_map
                                if(bvar.offset0 == True):  arr_map = isc_map-1

                                bvar.modlist[h_index].var_ar[arr_map] = \
                                                                    ( bvar.modlist[h_index].var_ar[arr_map] * bvar.modlist[h_index].var_n[arr_map] \
                                                                    + fp.variables[bvar.mod_symbol].data[itime,isc,site_index] ) \
                                                                    / (bvar.modlist[h_index].var_n[arr_map] + 1.0)
                                bvar.modlist[h_index].var_n[arr_map] = bvar.modlist[h_index].var_n[arr_map] + 1.0

                
                else:

                    print('Only scpf to sc is available now. Exiting.')
                    print(0)

        fp.close()

# =======================================================================================

class benchmark_vars:

    def __init__(self,name,mod_symbol,obs_symbol,mod_dimclass,obs_dimclass,unit):

        self.active     = False
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
            # Initialize scalar array
            self.sv  = True

        if(self.mod_dimclass == 'scalar'):

            # Initialize scalar array
            self.sv        = True


    def load_census(self,fp):


        # The 'size-class' type census variable is actually a
        # fp.variables[self.obs_symbol].shape
        # (7, 10, 3)
        # fp.variables[self.obs_symbol].dimensions
        # ('cens', 'dclass', 'confidence')

        #        fp.variables[self.obs_symbol].data

        d_sizes = fp.variables[self.obs_symbol].shape
        dim_names = fp.variables[self.obs_symbol].dimensions

        if(self.obs_dimclass == 'size-class'):

            if(dim_names[0] != 'cens'):
                print('expected census data to have cens as first dimension: {}'.format(dim_names))
                print('exiting')
                exit(2)
            
            # Condense the census dimension into 1 size
            self.scv_obs_ar      = np.zeros((d_sizes[1],d_sizes[2]))
            self.scv_obs_ar[:,:] = np.mean(fp.variables[self.obs_symbol].data, axis=0 )

            # Note that the dimensions in the census dictate the output dimension
            self.scv_x           = np.zeros((d_sizes[1],1))
            fp.variables['dclass'].data.resize(self.scv_x.shape)
            self.scv_x[:]        = fp.variables['dclass'].data
            self.scv_x_unit      = 'DBH [cm]'

        else:
            print('The census variable: {},   with dim-type: {}'.format(self.obs_symbol,self.obs_dimclass))
            print(' does not have dimensions set-up yet. Exiting')
            exit(2)


class mod_scv_array:

    def __init__(self,d_size):
    
        self.var_ar =  np.zeros((d_size,1))
        self.var_n  =  np.zeros((d_size,1))



def plot_bmarks(site,pdf):


    # Plot size structured benchmarks
    rline_types = ["k-","b--"]
    moline_types = ["ko-","b^--"]
    marksize = 10

    for bvar in bvarlist:

        plt.ioff()

        if( bvar.obs_dimclass == 'size-class' ):

            
            plt.plot(bvar.scv_x,bvar.scv_obs_ar[:,2])
        
            for imod, mod in enumerate(bvar.modlist):
                plt.plot(bvar.scv_x,mod.var_ar, \
                         "{}".format(moline_types[imod]),markersize=marksize, \
                         label="SAUCE{}")
