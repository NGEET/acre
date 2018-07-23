import code   # For development: code.interact(local = dict(globals(), **locals()))
import numpy as np
import xml.etree.ElementTree as et
from scipy.io import netcdf
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import math

# The CTFS processed files use an invalid flag of -9e+30
# We will consider anything very large negative invalid

invalid_flag = -9.9e10


# Anything that is a rate, needs to be normalized by the number of plants
# This is a restriction on 

nplant_scpf_name = 'NPLANT_SCPF'

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
                self.bvarlist.append(benchmark_vars(  name         = 'Basal Area', \
                                                      mod_symbols  = 'BA_SCPF', \
                                                      obs_symbol   = cens_var_name, \
                                                      mod_dimclass = 'scpf',      \
                                                      obs_dimclass = 'size-class', \
                                                      unit         = 'm2/ha', \
                                                      vartype      = 'quantity'))
                self.bvarlist[-1].load_census(fp)
            else:
                print('Census variable: '+cens_var_name+', was not found in the census file')

            cens_var_name = 'growth_increment_by_size_census'
            if (fp.variables.has_key(cens_var_name)):
                self.bvarlist.append(benchmark_vars(  name         = 'Growth Increment', \
                                                      mod_symbols  = 'DDBH_SCPF',  \
                                                      obs_symbol   = cens_var_name, \
                                                      mod_dimclass = 'scpf', \
                                                      obs_dimclass = 'size-class', \
                                                      unit         = 'cm/yr', \
                                                      vartype      = 'rate'))
                self.bvarlist[-1].load_census(fp)
            else:
                print('Census variable: '+cens_var_name+', was not found in the census file')


            cens_var_name = 'mortality_rate_by_size_census'
            if (fp.variables.has_key(cens_var_name)):
                self.bvarlist.append(benchmark_vars(  name         = 'Mortality Rate', \
                                                      mod_symbols  = 'M1_SCPF,M2_SCPF,M3_SCPF,M4_SCPF,M5_SCPF,M6_SCPF,M7_SCPF,M8_SCPF', \
                                                      obs_symbol   = cens_var_name, \
                                                      mod_dimclass = 'scpf',      \
                                                      obs_dimclass = 'size-class', \
                                                      unit         = '/yr', \
                                                      vartype      = 'rate'))
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
            
            all_symbols_found = True
            for mod_symbol in bvar.mod_symbols:
                if( not fp.variables.has_key(mod_symbol) ):
                    all_symbols_found = False

            if( all_symbols_found ):

                bvar.active = True

                if( bvar.mod_dimclass == 'scpf'):
                    
                    dims = fp.variables[bvar.mod_symbols[0]].dimensions
                    if(dims[1] != 'fates_levscpf'):
                        print('An SCPF benchmark variable: {} does not actually'.format(bvar.mod_symbols[0]))
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

        #code.interact(local = dict(globals(), **locals()))

        for bvar in self.bvarlist:

            if(bvar.active):

                if( (bvar.obs_dimclass=='size-class') and (bvar.mod_dimclass=='scpf') ):

                    # Create a mapping between FATES size-classes and the SCPF map
                    # ------------------------------------------------------------
                    fates_levscls       = fp.variables['fates_levscls'].data
                    fates_scmap_levscpf = fp.variables['fates_scmap_levscpf'].data
                    
                    
                    # fates_scmap_levscpf
                    # These are the expected dimensions
                    # ('time', 'fates_levscpf', 'lndgrid')
                    d_sizes = fp.variables[bvar.mod_symbols[0]].shape

                    hist_arr = np.ma.zeros(fp.variables[bvar.mod_symbols[0]].shape)

                    
                    if(bvar.vartype == 'quantity'):

                        for mod_symbol in bvar.mod_symbols:
                            hist_arr = hist_arr + fp.variables[mod_symbol].data

                    elif(bvar.vartype == 'rate'):
                        
                        #                        code.interact(local = dict(globals(), **locals()))
                        for mod_symbol in bvar.mod_symbols:
                            hist_arr = hist_arr + \
                                        fp.variables[mod_symbol].data / fp.variables[nplant_scpf_name].data
                            
                        # Mask out when there are no cohort counts
                        hist_arr[fp.variables[nplant_scpf_name].data <= 0.] = np.ma.masked
                        
                        # For rates... we are encountering some weird behavior where some files are not reporting
                        hist_arr[hist_arr.data <= 0.] = np.ma.masked

                    else:
                        print("Unhandled variable type submitted to registry: {}".format(bvar.vartype))
                        print("Must be one of: quantity or rate.  Exiting")
                        exit(2)

                    # Mask again if the variable has a no data flag
                    hist_arr[hist_arr.data<invalid_flag] = np.ma.masked

                    for itime in range(d_sizes[0]):

                        # Loop Sizes and then PFTs
                        # For quantities, add them
                        # For rates, take the mean
                        # code.interact(local = dict(globals(), **locals()))
                        # for isc,isc_map0 in enumerate(fates_scmap_levscpf):

                        for isc,isc_val in enumerate(fates_levscls):

                            if ( (bvar.offset0 == False) or ((bvar.offset0 == True)and(isc != 0 ))   ):

                                if(bvar.offset0==True):
                                    isc0 = isc-1
                                else:
                                    isc0 = isc

                                sc_maps = [i for i, x in enumerate(fates_scmap_levscpf) if x == isc+1 ]
                                
                                local_vars = hist_arr[itime,sc_maps,site_index]

                                if ( np.ma.count(local_vars)>0 ):

                                    if(bvar.vartype == 'quantity'):
                                        local_var = local_vars.sum()
                                    elif(bvar.vartype == 'rate'):
                                        local_var = local_vars.mean()
                                        print('entered {}'.format(local_vars))
                                    else:
                                        print('Unknown vartype')
                                        exit(2)

                                    bvar.modlist[h_index].var_ar[isc0] = ( bvar.modlist[h_index].var_ar[isc0] \
                                                                           * bvar.modlist[h_index].var_n[isc0] \
                                                                           + local_var) / (bvar.modlist[h_index].var_n[isc0] + 1.0)
                                    bvar.modlist[h_index].var_n[isc0] = bvar.modlist[h_index].var_n[isc0] + 1.0

                
                else:

                    print('Only scpf to sc is available now. Exiting.')
                    print(0)

        fp.close()

# =======================================================================================

class benchmark_vars:

    def __init__(self,name,mod_symbols,obs_symbol,mod_dimclass,obs_dimclass,unit,vartype):

        self.active     = False
        self.obs_symbol = obs_symbol
        self.name   = name
        self.scv    = False    # Size-class variable
        self.sv     = False    # Scalar variable
        self.unit   = unit
        self.vartype = vartype  # Is this a simple quantity, or a rate?
        
        self.hfile_id = -9

        # This will convert mod_symbols into a list
        self.mod_symbols = mod_symbols.split(',')
        

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

            # This is the mean across census intervals
            # AND... the lowest lower bound CI across census, 
            # and... the highest upper bound CI across census

            
            # Mask out bad data (for rates of change, probably missing
            # first census, or perhaps lowest
            masked_data = np.ma.array(fp.variables[self.obs_symbol].data, \
                                      mask=fp.variables[self.obs_symbol].data<invalid_flag)


            self.scv_obs_ar[:,0] = masked_data[:,:,0].min(axis=0).data
            self.scv_obs_ar[:,1] = masked_data[:,:,1].mean(axis=0).data
            self.scv_obs_ar[:,2] = masked_data[:,:,2].max(axis=0).data


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
    obs_mean_color = (0.6,0.6,0.6)
    obs_ci_color   = (0.8,0.8,0.8)

    moline_types = ["ko-","b^--"]
    marksize = 6

    #    code.interact(local = dict(globals(), **locals()))
    for bvar in site.benchmarks.bvarlist:

        fig, ax = plt.subplots()
        

        if( bvar.obs_dimclass == 'size-class' ):

            n_x         = bvar.scv_x.size
            obs_x       = np.resize(bvar.scv_x,(n_x))
            obs_mean    = np.resize(bvar.scv_obs_ar[:,1],(n_x))
            obs_ci_low  = np.resize(bvar.scv_obs_ar[:,0],(n_x))
            obs_ci_high = np.resize(bvar.scv_obs_ar[:,2],(n_x))

            # Plot out the observations, including a shaded CI
            ax.plot(obs_x,obs_mean,color = obs_mean_color,label="Census")
            ax.fill_between(obs_x, obs_ci_low , obs_ci_high,facecolor = obs_ci_color )

            print("OBS SUM:{}".format(np.sum(obs_mean)))

            for imod, mod in enumerate(bvar.modlist):
                ax.plot(bvar.scv_x,mod.var_ar, \
                         "{}".format(moline_types[imod]),markersize=marksize, \
                         label="Model")
                print("MOD SUM: {}".format(np.sum(mod.var_ar)))


            #            code.interact(local = dict(globals(), **locals()))
            ax.set_xlabel(bvar.scv_x_unit)
            ax.set_ylabel(bvar.unit)
            ax.set_title(bvar.name)
            ax.grid(True)

            ax.legend(loc='upper left')
        
        fig.suptitle("{}".format(site.name),fontsize=14,horizontalalignment='center')
        pdf.savefig(fig)
        plt.close(fig)
