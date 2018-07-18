# =======================================================================================
## @package acre_history_utils a bunch of scripts that help interperate


import code   # For development: code.interact(local=locals())
import numpy as np
import xml.etree.ElementTree as et
from scipy.io import netcdf
import matplotlib.dates as mdates
from datetime import datetime, timedelta

# =======================================================================================
## simple function that converts clm history file timing stamps, integer YYYYMMDD and 
# floating point seconds into integer year, month, date and hour of day.  Note that 
# the hour is returned on a floor function.  Time-stamps past the hour will refer
# backwards.  ***NOT SURE WHAT THE STAMPING CONVENTION FOR HIST FILES ARE
# @param yyyymmdd (integer) code where June 27, 1979 would be 19790627
# @param sec (float) number of seconds in the day
# @param date_bounds (float) bounds in days since reference of beg-end of time point
# @return year (integer) gregorian year
# @return moy (integer) month of year
# @return dom (integer) day of month
# @return hod (integer) hour of day (0-23)
def hist_dateint_to_num(yyyymmdd,sec,date_bounds):

    # Calculate the backwards adjustment
    adjustment = (date_bounds[0]-date_bounds[1])/2.0
    
    year = np.int16(np.floor(yyyymmdd/10000.0))          
    moy  = np.int16(np.floor((yyyymmdd-year*10000)/100.0))
    dom  = np.int16(yyyymmdd - (year*10000 + moy*100))
    hod  = np.int16(np.floor(sec/3600.0))

    adjusted = mdates.num2date(mdates.date2num(datetime(year,moy,dom,hod))+adjustment)

    year = adjusted.year
    moy  = adjusted.month
    dom  = adjusted.day
    hod  = adjusted.hour

    return(year,moy,dom,hod)


# =======================================================================================
## loop through the XML file (currently acre_constants.xml) entries containing definitions 
# of the history variables we are interested in analysing.  These definitions
# contain the name of the netcdf variable (as the tag) as well as the a list of codes 
# for the type of evaluation to perform.  The codes are the only entry for each tag, and 
# are strings with the codes comma delimited
# @param xmlfile (string) the name of the XML file
# @param hdims (class) data structure containing dimension info for the netcdf files
# @param file0 (string) the name of a netcdf file from which some info about the variable
# of interest can be found
# @return hvarlist (list) a list of variables of interest in the history files, stored
# as an hvar class
def define_histvars(xmlfile,hdims,file0,n_htypes,test_name,base_name):

    hvarlist = []
    siteroot = et.parse(xmlfile).getroot()

    ## ASSUME H0 for now
    hfile_id = 0
    
    for elem in siteroot.iter('hvar'):
        name = elem.attrib['tag']
        
        atypes = []
        units  = []
        mults  = []
        offs   = []
        for atype in elem.iter('atype'):
            atypes.append(atype.attrib['tag'])
            units.append(atype.find('unit').text)
            try:
                atype.find('mult').text
            except AttributeError:
                mults.append("1.0")
            else:
                mults.append(atype.find('mult').text)

            try:
                atype.find('offset').text
            except AttributeError:
                offs.append("0.0")
            else:
                offs.append(atype.find('offset').text)


        fp = netcdf.netcdf_file(file0, 'r', mmap=False)

        if (fp.variables.has_key(name)):

            # Determine dimension info
            dimnames = fp.variables[name].dimensions
#            units    = fp.variables[name]._attributes['units']
            hvarlist.append(hist_vars(name,atypes,units,mults,offs, \
                                      hdims,dimnames,hfile_id,n_htypes,test_name,base_name))
        else:
            print('History variable: '+name+', was not found in the history files')

        fp.close()

    return(hvarlist)

# =======================================================================================
## class containing several scratch space variables that will be used and re-used.  The 
#  space in the vectors needs to be large enough to accomodate the time dimension in
#  each history file. We allocate 125% of the space necessary in the first file.

class scratch_space:

    ## the constructor
    # @param file0 (string) a netcdf file from which time information can be diagnosed
    def __init__(self,file0):
        
        fp = netcdf.netcdf_file(file0, 'r',mmap=False)
        ntimemax = int(1.25*fp.variables['mcdate'].data.__len__())
        self.rawdata = np.zeros((ntimemax))
        self.movec   = np.zeros((ntimemax),dtype=np.int16)
        self.hrvec   = np.zeros((ntimemax),dtype=np.int16)
        self.yrvec   = np.zeros((ntimemax),dtype=np.int16)
        fp.close()

    ## the timing information in history files may be vectors, an alternative
    # formulation to read them in uses the scratch class to avoid memory allocation
    # this was supposed to be fast, but for some reason it deosn't save time
    # @param yyyymmdds (integer) year-month-day code, the "s" is plural
    # @param secs (float) seconds of the current date
    def hist_dateints_to_nums(self,yyyymmdds,secs,year_offset):

        nlen = secs.__len__()
        self.yrvec[:nlen] = (np.floor(yyyymmdds/10000.0)-year_offset).astype(int)
        self.movec[:nlen] = np.floor( (yyyymmdds - np.floor(yyyymmdds/10000.0)*10000.0) \
                                     /100.0).astype(int)
        self.hrvec[:nlen] = np.floor(secs/3600.0).astype(int)

# =======================================================================================
## Define the history variables (hvars) class type.  hvars is a class that contains 
# all the diagnostic variables and also some methods to fill those diagnostics from the 
# netcdf data files

class hist_vars:

    ## the constructor of the hist_vars class, and initialization of whatever components
    #  are available at the time of construction.
    # @param name (string) the short name of the variable as defining the netcdf entry
    # @param atypes (list) a list of the types of averageing available for this var
    # @param units (list) a list of the units associated with the different atypes
    # @param mults (list) a list of the multiplier used to achieve the units in atype
    # @param offs (list) a list of the offset used to achieve the units of atype
    # @param hdims (class) contains information on the dimension info of the history
    #                      files
    # @param dimnames (list of strings) the names, in order, of the dimensions of which
    #        the current history variable contains
    # @param hfile_id (integer) [NOT SET CURRENTLY] the history file identifier
    #        this variable is found in, ie H0, H1, H2, up to H6
    # @param n_htypes (integer) defines the size of allocated arrays, where 0 is for
    #        non-regression evaluations, and when regression comparisons against a
    #        baseline are ordered, n_htypes = 1
    def __init__(self,name,atypes,units,mults,offs,hdims,dimnames, \
                 hfile_id,n_htypes,test_name,base_name):

        ## Name in the netcdf file
        self.name     = name
        self.mmv      = False
        self.dmv      = False
        self.amv      = False
        self.amds     = False
        ## history file that contains the variable (h0,h1,etc)
        #  only h0 is availble right now
        self.hfile_id = hfile_id
        self.test_name_str = [test_name,base_name]
        # time dimension
        if('time' in dimnames):
            self.timedim = [id for id in range(dimnames.__len__())\
                            if dimnames[id] in 'time'][0]
        else:
            print('History Variable: '+name+'has no time dimension?')

        # Classify the variable

        # 2D Classes:
        # Current allowable classes:  2DLND  ('time','lndgrd')
        #                             3DLND  ('time','lat','lon')

        self.dimclass = None
        
        if(dimnames.__len__() == 2):
            if(dimnames[0] == 'time' and dimnames[1] == 'lndgrid'):
                self.dimclass = '2dlnd'

        if(dimnames.__len__() == 3):
            if(dimnames[0] == 'time' and dimnames[1] == 'lat' and dimnames[2] == 'lon'):
                self.dimclass = '3dlnd'

        if(dimnames.__len__() == 3):
            if(dimnames[2] == 'lndgrid' and dimnames[1] == 'fates_levscpf'):
                self.dimclass = '3dlndscpf'

        if(self.dimclass is None):
            print('History Variable: '+name+' does not have a registered dimensionality')
            print('in this toolset.  ')
            if(dimnames.__len__()>0):
                for dimname in dimnames:
                    print(dimname)
            exit()

        # Determine what types of output formats are requested
        # for this variable.  
        # Only allow monthly analysis if the history files contain monthly or finer data
        # Only allow diurnal analysis if the history files contain 6 hourly or finer data

        for atype,unit,mult,offset in zip(atypes,units,mults,offs):
            if(atype=="MMV"):
                self.mmv = True
                self.mmv_ar = np.zeros((12,n_htypes))
                self.mmv_n  = np.zeros((12,n_htypes))
                if(n_htypes==1):
                    self.mmv_x  = np.array([range(1,13)]).transpose()
                else:
                    self.mmv_x  = np.array([range(1,13),range(1,13)]).transpose()
                self.mmv_unit = unit
                self.mmv_mult = float(mult)
                self.mmv_off  = float(offset)

            if(atype=="DMV"):
                self.dmv = True
                self.dmv_ar = np.zeros((24,n_htypes))
                self.dmv_n  = np.zeros((24,n_htypes))
                if(n_htypes==1):
                    self.dmv_x  = np.array([range(24)]).transpose()
                else:
                    self.dmv_x  = np.array([range(24),range(24)]).transpose()
                self.dmv_unit = unit
                self.dmv_mult = float(mult)
                self.dmv_off  = float(offset)
            if(atype=="AMV"):
                self.amv = True
                self.amv_ar   = np.zeros((hdims.nyears,n_htypes))
                self.amv_n    = np.zeros((hdims.nyears,n_htypes))
                if(n_htypes==1):
                    self.amv_x  = np.array([range(hdims.yeara,hdims.yearz+1)]).transpose()
                else:
                    self.amv_x  = np.array([range(hdims.yeara,hdims.yearz+1), \
                                            range(hdims.yeara,hdims.yearz+1)]).transpose()
                self.amv_unit = unit
                self.amv_mult = float(mult)
                self.amv_off  = float(offset)


    ## push the values of a given time-stamp (rawdata) into the MONTHLY average
    # @param rawdata (float 1D-vector) the raw data in the file
    # @param movec (integer 1D-vector) the months associated with the raw data
    # @param htype (integer) specifies if this is test or baseline data
    def push_mmv(self,rawdata,movec,htype):
        for idx, val in enumerate(movec):
            self.mmv_ar[int(val),htype] = self.mmv_ar[int(val),htype] + \
                                            rawdata[idx]*self.mmv_mult + self.mmv_off
            self.mmv_n[int(val),htype]  = self.mmv_n[int(val),htype]  + 1.0

    ## push the values of a given time-stamp (rawdata) into the DIURNAL average
    # @param rawdata (float 1D-vector) the raw data in the file
    # @param hrvec (integer 1D-vector) the hours associated with the raw data
    # @param htype (integer) specifies if this is test or baseline data
    def push_dmv(self,rawdata,hrvec,htype):
        for idx, val in enumerate(hrvec):
            self.dmv_ar[int(val),htype] = self.dmv_ar[int(val),htype] + \
                                          rawdata[idx]*self.dmv_mult + self.dmv_off
            self.dmv_n[int(val),htype]  = self.dmv_n[int(val),htype]  + 1.0

    ## push the values of a given time-stamp (rawdata) into the ANNUAL average
    # @param rawdata (float 1D-vector) the raw data in the file
    # @param anvec (integer 1D-vector) the years associated with the raw data
    # @param htype (integer) specifies if this is test or baseline data   
    def push_amv(self,rawdata,yrvec,htype):
        for idx, val in enumerate(yrvec):
            self.amv_ar[int(val),htype] = self.amv_ar[int(val),htype] + \
                                          rawdata[idx]*self.amv_mult + self.amv_off
            self.amv_n[int(val),htype]  = self.amv_n[int(val),htype]  + 1.0

    ## like push_mmv but does vectorized processes across unique month entries
    def push_mmvvec(self,rawdata,movec,htype):
        uns = np.unique(movec)
        for idx,val in enumerate(uns):
            mask = movec == val
            self.mmv_ar[int(val),htype] = self.mmv_ar[int(val),htype] + \
                np.sum(rawdata[mask])*self.mmv_mult + self.mmv_off
            self.mmv_n[int(val),htype]  = self.mmv_n[int(val),htype]  + np.sum(mask)

    ## like push_dmv but does vectorized processes across unique day entries
    def push_dmvvec(self,rawdata,hrvec,htype):
        uns = np.unique(hrvec)
        for idx,val in enumerate(uns):
            mask = hrvec == val
            self.dmv_ar[int(val),htype] = self.dmv_ar[int(val),htype] + \
                np.sum(rawdata[mask])*self.dmv_mult + self.dmv_off
            self.dmv_n[int(val),htype]  = self.dmv_n[int(val),htype]  + np.sum(mask)

    ## like push_amv but does vectorized processes across unique month year entries
    def push_amvvec(self,rawdata,yrvec,htype):
        uns = np.unique(yrvec)
        for idx,val in enumerate(uns):
            mask = yrvec == val
            self.amv_ar[int(val),htype] = self.amv_ar[int(val),htype] + \
                np.sum(rawdata[mask])*self.amv_mult + self.amv_off
            self.amv_n[int(val),htype]  = self.amv_n[int(val),htype]  + np.sum(mask)

    ## after all of the raw data has been pushed into their respective averaging arrays,
    #  the sums need to be converted into averages by simple normalization.
    #  Note that there may be cases where no-data was pushed to a certain index.  This
    #  is fine, and a NaN will show up in the plots with no-data.  To allow this
    #  without error, the numpy error flag is by-passed and re-set.
    def normalize_diagnostics(self):

        old_div0 = np.seterr(divide='ignore',invalid='ignore')
        if(self.mmv):
            self.mmv_ar = np.divide(self.mmv_ar,self.mmv_n)

        if(self.dmv):
            self.dmv_ar = np.divide(self.dmv_ar,self.dmv_n)
            
        if(self.amv):
            self.amv_ar = np.divide(self.amv_ar,self.amv_n)
            mask = self.amv_n < 2 
            self.amv_ar[mask] = np.nan
            
            # Mask out values that are incomplete
            maxmask = max(self.amv_n[:,0])
            mask = self.amv_n < np.floor(0.75*maxmask)
            self.amv_ar[mask] = np.nan

        np.seterr(**old_div0)

# =======================================================================================
## This is the main call to process data for any given history file.  A single file is
#  passed as an argument, and this file cycles through the different variables it has
#  in its database, and based on their definitions, will process data for each of them
#  accordingly.  
# @param file (string) the netcdf file name currently being opened
# @param site (class) the site class object (for spatial indexing)
# @param hvarlist (list) the list of history variables of class (hvar)
# @param htype (integer) 0 for test version, 1 for baseline version
# @param scr (class) the scratch space for mathing things hard
# @param hdims (class) information about the file dimensions
def load_history(file,site,hvarlist,htype,scr,hdims):

    vectorizedates = False

    # Load up a file to retrieve dimension info
    fp = netcdf.netcdf_file(file, 'r', mmap=False)
    
    yyyymmdd    = fp.variables['mcdate'].data
    sec         = fp.variables['mcsec'].data
    date_bounds = fp.variables['time_bounds'].data    #[time,bounds_interval]
    ntimes      = int(yyyymmdd.__len__())

    print('Loading: '+file)

    if(~vectorizedates):
        # Strangely, when the time operations were vectorized, the code
        # took longer to complete.  I will leave the vectorized calls
        # in place and perhaps I will figure out why some day
        for it in range(ntimes):
            year,moy,dom,hod = hist_dateint_to_num(yyyymmdd[it],sec[it],date_bounds[it,:])
            scr.yrvec[it] = np.int16(year)-hdims.yeara
            scr.movec[it] = np.int16(moy)-1
            scr.hrvec[it] = np.int16(hod)-1

    else:
        scr.hist_dateints_to_nums(yyyymmdd,sec,hdims.yeara)

    for hvar in hvarlist:
        # IF THIS IS A 1D:
        if(hvar.dimclass=='2dlnd'):
            rawshape = fp.variables[hvar.name].shape
            if(rawshape[1] == 1):
                scr.rawdata[:ntimes] = \
                    fp.variables[hvar.name].data[:ntimes].reshape(-1)
            else:
                scr.rawdata[:ntimes] = \
                    fp.variables[hvar.name].data[:ntimes,site.igh].reshape(-1)

            if(hvar.mmv):
                hvar.push_mmvvec(scr.rawdata[:ntimes],scr.movec[:ntimes],htype)
                
            if(hvar.dmv):
                hvar.push_dmv(scr.rawdata[:ntimes],scr.hrvec[:ntimes],htype)

            if(hvar.amv):
                hvar.push_amvvec(scr.rawdata[:ntimes],scr.yrvec[:ntimes],htype)

        if(hvar.dimclass=='3dlnd'):
            rawshape = fp.variables[hvar.name].shape
            if(rawshape[1] == 1):
                scr.rawdata[:ntimes] = \
                    fp.variables[hvar.name].data[:ntimes].reshape(-1)
            else:
                scr.rawdata[:ntimes] = \
                    fp.variables[hvar.name].data[:ntimes,site.ilath,site.ilonh].reshape(-1)

            if(hvar.mmv):
                hvar.push_mmvvec(scr.rawdata[:ntimes],scr.movec[:ntimes],htype)

            if(hvar.dmv):
                hvar.push_dmv(scr.rawdata[:ntimes],scr.hrvec[:ntimes],htype)

            if(hvar.amv):
                hvar.push_amvvec(scr.rawdata[:ntimes],scr.yrvec[:ntimes],htype)

        if(hvar.dimclass=='3dlndscpf'):
            rawshape = fp.variables[hvar.name].shape
            scpf_dim = 1   # 0 is time, and 2 is space

            # ----------------------------------------------------------------
            # In some cases we want to just condense the second dimension
            # the "scpf" dimension into a single value.
            # This is a 1D output, it is acceptable to condense the 
            # second [1] dimension.  The only question is whether or not the 
            # variable is a mean, normalized (like per plant), or a total
            # 1 = sum
            # 2 = mean
            # ----------------------------------------------------------------

            if(fp.variables[hvar.name].units=='kgC/m2/yr'):
                scale_type=1
            if(fp.variables[hvar.name].units=='m2/ha'):
                scale_type=1
            else:
                print('Rescaling a 3d variables requires some understanding')
                print('of its units. Type {} for variable {} is unknown'.format( \
                    fp.variables[hvar.name].units,hvar.name))
                exit(0)

            if(hvar.mmv ):
                if(scale_type==1):
                    scr.rawdata[:ntimes] = np.sum(fp.variables[hvar.name].data[:ntimes,:,site.igh],axis=scpf_dim).reshape(-1)
                    hvar.push_mmvvec(scr.rawdata[:ntimes],scr.movec[:ntimes],htype)

            if(hvar.dmv):
                if(scale_type==1):
                    scr.rawdata[:ntimes] = np.sum(fp.variables[hvar.name].data[:ntimes,:,site.igh],axis=scpf_dim).reshape(-1)
                    hvar.push_dmv(scr.rawdata[:ntimes],scr.hrvec[:ntimes],htype)

            if(hvar.amv):
                if(scale_type==1):
                    scr.rawdata[:ntimes] = np.sum(fp.variables[hvar.name].data[:ntimes,:,site.igh],axis=scpf_dim).reshape(-1)
                    hvar.push_amvvec(scr.rawdata[:ntimes],scr.yrvec[:ntimes],htype)


    fp.close()

    


# =======================================================================================        
## history file dimensions
class hist_dims:

    # the constructor of hist_dims
    ## @param file0 (string) an arbitrary file name used to extract dimension data
    def __init__(self,file0):

        fp = netcdf.netcdf_file(file0, 'r', mmap=False)
        
        print('Dimension information for history files\n')
 
        ## This is the list of all dates
        #  when fully constructed, it will have a nested
        #  list with sub-file dates, within the outer list
        #  of file indices.  It presents as a 2d list

        self.datelist = []
        self.hperiod  = -9
        self.ntimes   = -9
        self.nyears   = -9

        ## Size dimension
        self.dbh = fp.variables['fates_levscls'].data

        ## A mapping from the scpf vector to pft index
        self.pft_map     = fp.variables['fates_pftmap_levscpf'].data
        
        ## A mapping from the scpf vector to size index
        self.scls_map    = fp.variables['fates_scmap_levscpf'].data
        
        fp.close()

    ## process timing information for a specific filetype, where the history file list
    # may have up to 3 files.  The first file and possibly the second (if each file has 
    # 1 time-stamp) are used to determine the temporal frequency of the time dimensions.
    # The last file is used in conjunction with the first file to determine the total
    # range in years of output, this is necessary for defining array sizes.
    # @param filelist (list) a list of netcdf file names (probably first, second last)
    def timing(self,filelist):

        if(filelist.__len__()==3):
            filea = filelist[0]
            fileb = filelist[1]
            filez = filelist[2]
        elif(filelist.__len__()==2):
            filea = filelist[0]
            fileb = filelist[1]
            filez = filelist[1]
        elif(filelist.__len__()==1):
            filea = filelist[0]
            fileb = filelist[0]
            filez = filelist[0]
        else:
            print('Retrieving file timing info from history')
            print('files requires first, second and last files.')
            print('You provided {} number of files'.format(threefilelist.__len__()))
            sys.exit(2)

        # Open the timing info on the first file
        fpa = netcdf.netcdf_file(filea, 'r', mmap=False)
        yyyymmdd_a = fpa.variables['mcdate'].data
        sec_a      = fpa.variables['mcsec'].data
        db_a       = fpa.variables['time_bounds'].data
        ntimes_a   = (yyyymmdd_a.__len__());
        
        # Open the timing info on the second file
        fpb = netcdf.netcdf_file(fileb, 'r', mmap=False)
        yyyymmdd_b = fpb.variables['mcdate'].data
        sec_b      = fpb.variables['mcsec'].data
        db_b       = fpb.variables['time_bounds'].data
        
        # Open the year info for just the last file
        fpz        = netcdf.netcdf_file(filez, 'r', mmap=False)
        yyyymmdd_z = fpz.variables['mcdate'].data
        sec_z      = fpz.variables['mcsec'].data
        db_z       = fpz.variables['time_bounds'].data
        
        if(ntimes_a>1):  # Multiple time-stamps per file
            
            yr1,moy1,dom1,hod1 = hist_dateint_to_num(yyyymmdd_a[0],sec_a[0],db_a[0,:])
            yr2,moy2,dom2,hod2 = hist_dateint_to_num(yyyymmdd_a[1],sec_a[1],db_a[1,:])

        else:
            
            yr1,moy1,dom1,hod1 = hist_dateint_to_num(yyyymmdd_a[0],sec_a[0],db_a[0,:])
            yr2,moy2,dom2,hod2 = hist_dateint_to_num(yyyymmdd_b[0],sec_b[0],db_b[0,:])

        # This does not need to be perfect at large numbers
        # hourly period between history time-stamps
        self.hperiod = int(np.ceil( (yr2-yr1)*365*24 + \
                                      (moy2-moy1)*30*24 + \
                                      (dom2-dom1)*24 + 
                                      (hod2-hod1) ))

                                    
        yrz,moyz,domz,hodz = hist_dateint_to_num(yyyymmdd_z[-1],sec_z[-1],db_z[-1,:])

#        yrz  = int(np.floor(yyyymmdd_z[-1]/10000.0))  

        self.nyears = yrz-yr1+1
        self.yeara = yr1
        self.yearz = yrz

        fpa.close()   # Close the netcdf file
        fpb.close()   # Close the netcdf file
        fpz.close()   # Close the netcdf file
