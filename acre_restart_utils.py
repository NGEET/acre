import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join 
import sys
import code  # For development: code.interact(local=locals())
import datetime

## Define the restart variables (rvars) class type
#  rvars is a class that contains diagnostic variables
#  and also methods to fill those diagnostics from the netcdf data files
class rvars:

    def __init__(self,n_rtypes,n_r_files,test_name,base_name):

        ## Total dead biomass within the restart file's 
        #  grid-cell of interest (GOI) [MgC/ha]
        self.bdead_gc      = np.zeros((n_rtypes,n_r_files))
        
        ## Total AGB within the restart file's GOI [MgC/ha]
        self.agb_gc        = np.zeros((n_rtypes,n_r_files))
    
        ## Total leaf biomass within the restart file's GOI [MgC/ha]
        self.bleaf_gc      = np.zeros((n_rtypes,n_r_files))
        
        ## Total basal area within the restart file's GOI [MgC/ha]
        self.ba_gc         = np.zeros((n_rtypes,n_r_files))
        
        ## Total dead biomass among trees less than 5cm diameter,
        # within the restart file's GOI [MgC/ha]
        self.bdead_dbh5_gc = np.zeros((n_rtypes,n_r_files))
        
        ## Total dead biomass among trees less than 5cm diameter,
        # within the restart file's GOI [MgC/ha]
        self.bdead_dbh2_gc = np.zeros((n_rtypes,n_r_files))

        ## The file index and array counter for the different filetypes
        self.fid           = np.zeros((n_rtypes))

        ## The mean number of cohorts per patch
        self.nc_per_pa_gc  = np.zeros((n_rtypes,n_r_files))

        ## The number of patches per column
        self.npa_gc        = np.zeros((n_rtypes,n_r_files))

        self.test_name_str = [test_name,base_name]

    ## load_restart reads and processes data in individual restart files. 
    # This includes opening the netcdf file, retrieving a list of variables, 
    # incrementing and averaging thos variables for this time-stamp, closing the
    # file and returning variables.
    # @param file The filename of the current restart file, without path
    # @param file_prefix The path to the filename
    # @param igr The grid index containing the SOI
    # @return rtot_bdead,rtot_bleaf,rtot_ba,rtot_bdead_recr5,rtot_bdead_recr2
    def load_restart(self,file,igr,rtype):
        from scipy.io import netcdf


        ## Temporary Flag
        gridstruct = False

        ## A no data flag consistent with CESM/ACME restart no data flags
        nodataflag = -9999

        ## Conversion multiplier for square meters to hectares
        m2_to_ha = 10000.0  
        
        ## Conversion multiplier for kilograms to Megagrams
        kg_to_Mg  = 0.001   
    
        ## Conversion multiplier for square centimeters to square meters
        cm2_to_m2 = 0.0001

        # Alias the index counter for arrays files
        fid = self.fid[rtype]
   
        # Load the file
        print('Loading: '+file)
        fp = netcdf.netcdf_file(file, 'r')
 
        co_ar_size = fp.dimensions['cohort']

        # Test to see if this is the old or new method

        if ('fates_PatchesPerSite' in fp.variables.keys()):

            patchpcell = fp.variables['fates_PatchesPerSite'].data
            if( patchpcell[igr]<=0 ):
                print('There arent any ed patches at the specified coordinates')
                sys.exit()

            gid_b = igr*co_ar_size   # Cohort vector first index for gridcell
            gid_e = gid_b+co_ar_size     # Cohort vector last index for gridcell

            # In this call, we identify the vector of cohort and patch data
            # (which currrently use the same space) associated with a grid cell
            patchmap  = fp.variables['fates_CohortsPerPatch'].data[gid_b:gid_e]

            # array indices
            pa_ids = np.where(np.greater(patchmap,nodataflag))  

            # total area of all FATES patches
            max_area = sum(fp.variables['fates_area'].data[gid_b + pa_ids[0][:]])

            # Patch loop
            for pid in gid_b + pa_ids[0][:]:     # This is a tuple I think

                # Retrieve patch level variables
                p_area  = fp.variables['fates_area'].data[pid]
                p_age   = fp.variables['fates_age'].data[pid]
        
                # Retrieve cohort level variables
                n_cohorts = patchmap[pid]
                cid_b = pid
                cid_e = pid+n_cohorts  # python uses the last index as a greater than type

                c_bleaf  = fp.variables['fates_bl'].data[cid_b:cid_e]
                c_bdead  = fp.variables['fates_bdead'].data[cid_b:cid_e] # [kgC/plant]
                c_dbh    = fp.variables['fates_dbh'].data[cid_b:cid_e]
                c_heigt  = fp.variables['fates_height'].data[cid_b:cid_e]
                c_nplant = fp.variables['fates_nplant'].data[cid_b:cid_e]    # [plant/patch]
                c_pft    = fp.variables['fates_pft'].data[gid_b:gid_e]

                # If any of the cohorts are zero, then something is wrong
                if (any(np.less_equal(c_nplant,0.0))):
                    print('cohort map believes a zero density cohort is valid, exiting');
                    sys.exit()
    
                # Since N is plants per patch, we don't even need to do area weighting
                # we just simply add the stuff up

                self.bdead_gc[rtype,fid] += kg_to_Mg * sum(c_bdead*c_nplant)
                self.bleaf_gc[rtype,fid] += sum(c_bdead*c_nplant)
                self.ba_gc[rtype,fid]    += sum(0.25*np.pi*cm2_to_m2*c_dbh*c_dbh*c_nplant) 

                # Filter those cohorts smaller than 5cm dbh
                idrecr5 = np.less_equal(c_dbh,5.0)
                self.bdead_dbh5_gc[rtype,fid] += kg_to_Mg * sum(c_bdead[idrecr5]*c_nplant[idrecr5])

                idrecr2 = np.less_equal(c_dbh,2.0)
                self.bdead_dbh2_gc[rtype,fid] += kg_to_Mg * sum(c_bdead[idrecr2]*c_nplant[idrecr2])

                self.nc_per_pa_gc[rtype,fid] += n_cohorts*(p_area/max_area)

                self.npa_gc[rtype,fid] += 1


        else:

            # Array of gridcell indices for each column
            col2gr = fp.variables['cols1d_gridcell_index'].data

            # Array of number of patches in those grid-cells
            patchpcol = fp.variables['fates_PatchesPerSite'].data

            # Column weight on that grid-cell
            cols1d_wtxy = fp.variables['cols1d_wtxy'].data

            # cols1d_active
            cols1d_active = fp.variables['cols1d_active'].data


            # We want to filter columns that both have patches
            # and also match the grid-cell of interest
            
            icrs = np.where( (patchpcol>nodataflag) & (col2gr==(igr+1)) & (cols1d_wtxy>0.0 ) & (cols1d_active==1))

            sumwt = sum(cols1d_wtxy[icrs])

            for idx, icr in enumerate(icrs):

                colwt = cols1d_wtxy[icr]/sumwt
                gid_b = icr*co_ar_size
                gid_e = gid_b + co_ar_size

                # In this call, we identify the vector of cohort and patch data
                # (which currrently use the same space) associated with a grid cell
                patchmap  = fp.variables['fates_CohortsPerPatch'].data[gid_b:gid_e]

                # array indices
                pa_ids = np.where(np.greater(patchmap,nodataflag))  

                max_area = sum(fp.variables['fates_area'].data[gid_b + pa_ids[0][:]])
                
                # Patch loop
                for pid in gid_b + pa_ids[0][:]:     # This is a tuple I think

                    # Retrieve patch level variables
                    p_area  = fp.variables['fates_area'].data[pid]
                    p_age   = fp.variables['fates_age'].data[pid]
        
                    # Retrieve cohort level variables
                    n_cohorts = patchmap[pid]
                    cid_b = pid
                    cid_e = pid+n_cohorts  # python uses the last index as a greater than type

                    c_bleaf  = fp.variables['fates_bl'].data[cid_b:cid_e]
                    c_bdead  = fp.variables['fates_bdead'].data[cid_b:cid_e] # [kgC/plant]
                    c_dbh    = fp.variables['fates_dbh'].data[cid_b:cid_e]
                    c_heigt  = fp.variables['fates_height'].data[cid_b:cid_e]
                    c_nplant = fp.variables['fates_nplant'].data[cid_b:cid_e]    # [plant/patch]
                    c_pft    = fp.variables['fates_pft'].data[cid_b:cid_e]

                    # If any of the cohorts are zero, then something is wrong
                    if (any(np.less_equal(c_nplant,0.0))):
                        print('cohort map believes a zero density cohort is valid, exiting');
                        sys.exit()
    
                    # Since N is plants per patch, we don't even need to do area weighting
                    # we just simply add the stuff up

                    self.bdead_gc[rtype,fid] += kg_to_Mg * sum(c_bdead*c_nplant) * colwt
                    self.bleaf_gc[rtype,fid] += sum(c_bdead*c_nplant) * colwt
                    self.ba_gc[rtype,fid]    += sum(0.25*np.pi*cm2_to_m2*c_dbh*c_dbh*c_nplant) * colwt

                    # Filter those cohorts smaller than 5cm dbh
                    idrecr5 = np.less_equal(c_dbh,5.0)
                    self.bdead_dbh5_gc[rtype,fid] += kg_to_Mg * sum(c_bdead[idrecr5]*c_nplant[idrecr5]) * colwt

                    idrecr2 = np.less_equal(c_dbh,2.0)
                    self.bdead_dbh2_gc[rtype,fid] += kg_to_Mg * sum(c_bdead[idrecr2]*c_nplant[idrecr2]) * colwt

                    self.nc_per_pa_gc[rtype,fid] += n_cohorts*(p_area/max_area)
                    
                    self.npa_gc[rtype,fid] += 1



        

        fp.close()   # Close the netcdf file
        
        # Increment the file array counter
        self.fid[rtype] += 1




