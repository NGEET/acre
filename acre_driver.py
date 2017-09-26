# =======================================================================================
## @package acre_driver a driver for the Answer Changing Regression Evaluator
#
# For usage: $python acre_driver.py -h
#
#
#
# This script is intended to diagnose the output of several single site runs, as a 
# rapid visual and textual analysis on ecosystem response over a period of 5 - 50 years.  
#
# As a default, the user must provide the file-path to restart files. As a default,
# plotting is turned off.
#
# Options include 1) the ability to do a regression against another set of files (base).
#                 2) analysis of history files (which have a 1-hour suggested output)
#                 3) plotting (most of the analysis right now is visual, not much
#                              really happens right now without plotting)
#
# Next capabilities that are in queue:
#     - add more restart output diagnostics
#     - come up with a list of PASS/FAIL diagnostics
#     - h1,h2,h3 files
#     - add more checks
#     - optimize math and averaging
#
#
#  FYI: REALLY HELPFULL DEBUG TOOL
#  UNcomment the imported "code" library and add the following call
#  in the code where you want to add a debug stop:
#        code.interact(local=locals())
#
#
#
# =======================================================================================

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import sys
import getopt
import code  # For development: code.interact(local=locals())
import time
import acre_history_utils as hutils
import acre_restart_utils as rutils
import acre_plot_utils as putils
import acre_table_utils as tutils

# =======================================================================================
# Some Global Parmaeters

## The name of the xml file containing site data (should not change)
xmlfile = 'acre_controls.xml'

## A proximity threshold for matching the locations of Sites Of Interest (SOIs)
# with locations with the model output grids, units: [degrees] 
geo_thresh = 2.0
        


# =======================================================================================

## Define the sitetype class structure
#  Sitetype is a class used to describe Sites Of Interest (SOIs).
#  SOIs are associated with geographic locations on the globe
#  where we want to closely evaluate model output.  These sites may
#  be chosen because they present unusual ecology, or because these
#  location contains instrumentation or a monitoring network.
# \callgraph
class sitetype:
    
    ## The constructor only initializes name and location.
    # @param self The object
    # @param name The name (text string) of the SOI
    # @param lat  The latitude of the SOI in decimal degrees
    # @param lon  The longitude of the SOI in decimal degrees
    def __init__(self,name,lat,lon):
        
        ## The name (text string) of the SOI
        self.name = name  
        ## The latitude of the SOI in decimal degrees
        self.lat  = lat
        ## The longitude of the SOI in decimal degrees
        self.lon  = lon    
        ## The associated grid-index of the restart file (index grid restart = igr)
        self.igr  = -9     
        ## The associated grid-index of the history file (index grid history = igh)
        self.igh  = -9
        ## If a gridded simulation, instead of lndgrid we have 2d geographic coordinates
        self.ilath = -9
        self.ilonh = -9

        ## This index specifies what type of history grid we have
        self.hgrid = -9

        
# =======================================================================================


## get a list of files given a directory prefix and specification on the type
# @param file_prefix a string with the full or relative path to the data
# @return filelist a list of strings, each of which a file
def getnclist(file_prefix,filetype):

    from os.path import isfile, join 
    from os import listdir

    file_prefix = file_prefix+'/'

    if(filetype == 'restart'):
        filelist = sorted([file_prefix+f for f in listdir(file_prefix) \
                               if ( isfile(join(file_prefix, f)) & \
                                        ('clm2.r.' in f) & ('.nc' in f ))  ])

    elif(filetype == 'h0'):
        filelist = sorted([file_prefix+f for f in listdir(file_prefix) \
                               if ( isfile(join(file_prefix, f)) & \
                                        ('clm2.h0.' in f) & ('.nc' in f ))  ])

    elif(filetype == 'h1'):
        filelist = sorted([file_prefix+f for f in listdir(file_prefix) \
                               if ( isfile(join(file_prefix, f)) & \
                                        ('clm2.h1.' in f) & ('.nc' in f ))  ])        

    elif(filetype == 'h2'):
        filelist = sorted([file_prefix+f for f in listdir(file_prefix) \
                               if ( isfile(join(file_prefix, f)) & \
                                        ('clm2.h2.' in f) & ('.nc' in f ))  ])

    else:

        print("filetype: {}, is not supported".format(filetype))
        sys.exit(2)

    return(filelist)


# ========================================================================================

## load_sites interprets the Site of Interest (SOI) XML file.
# It does this by loading the XML entries automatically into an element tree.
# The generated tree is fairly simple, and the important tags are stored.  
# There is only one parent entry class "site". Within which is the name 
# (its tag), and the lat/lon coordinates.
# @param xmlfile the xml filename including path
# @param sitetype the class type
# @return sites a list of class sitetype

def load_sites(xmlfile,sitetype):
    
    import xml.etree.ElementTree as et
    sites = []

    siteroot = et.parse(xmlfile).getroot()
    for elem in siteroot.iter('site'):
        
        name = elem.attrib['tag']
        lat  = float(elem.find('lat').text)
        lon  = float(elem.find('lon').text)
        if(lon<0):
            lon = 360.0+lon
        sites.append(sitetype(name,lat,lon))

    return(sites)

# ========================================================================================

## print the help message for argument passing

def usage():
     print('')
     print('=======================================================================')
     print('')
     print(' python acre_driver.py -h --plotmode --regressmode --restartmode')
     print('                         --test-hist-pref=<path> --base-hist-pref=<path>')
     print('                         --test-rest-pref=<path> --base-rest-pref=<path>')
     print('                         --test-name=<test> --base-name=<text>')
     print('')
     print('  This script is intended to diagnose the output of several single site')
     print('  runs, as a rapid pass/fail visual analysis on ecosystem response over')
     print('  time.  Simulations have on been tested on single site runs, but were')
     print('  designed to accomodate gridded output as well.  The key requisit is to ')
     print('  have available, multiple time points of history output.')
     print('')
     print('')
     print(' -h --help ')
     print('     print this help message')
     print('')
     print(' --plotmode')
     print('     [Optional] logical switch, turns on plotting')
     print('     default is False')
     print('')
     print(' --regressmode')
     print('     [Optional] logical switch, turns on regression tests')
     print('     against a baseline. Requires user to also set --base-rest-pref')
     print('     default is False')
     print('')
     print(' --restartmode')
     print('     [Optional] logical switch, turns on evaluations of restart type output')
     print('     in this mode, the code will search for .r. files and key FATES demogrpahics.')
     print('     --------------------------------------- ')
     print('     DO NOT USE THIS MODE ON NON-FATES OUTPUT')
     print('     --------------------------------------- ')
     print('')
     print(' --eval-id=<id-string>')
     print('     a string that gives a name, or some id-tag associated with the')
     print('     evaluation being conducted. This will be used in output file naming.')
     print('     Any spaces detected in string will be trimmed.')
     print('')
     print(' --test-hist-pref=<path>')
     print('     the full path to the folder with history files of the test')
     print('     version of output')
     print('') 
     print(' --base-hist-pref=<path>')
     print('     [Optional]  the full path to history files of a baseline')
     print('     version of output')
     print('')
     print(' --test-rest-pref=<path>')
     print('     [Optional] the full path to the folder with the restart files of ')
     print('     the test version of output')
     print('') 
     print(' --base-rest-pref=<path>')
     print('     [Optional] the full path to the folder with the restart files of ')
     print('     the base version of output')
     print('') 
     print(' --test-name=<text>')
     print('     [Optional] a short descriptor for the test case that will be used')
     print('     for labeling plots. The default for the test case is "test".')
     print('')
     print(' --base-name=<text>')
     print('     [Optional] a short descriptor for the base case that will be used')
     print('     for labeling plots. The default for the base case is "base".')
     print('')
     print('')
     print('=======================================================================')

# =======================================================================================

## interpret timing information from restart files.
# This information is appended to a list of "datelist" objects that are defined in the 
# datetime library.
# @param file The filename of the current restart file, without path
# @param file_prefix The path to the filename
# @param datelist A special date class defined by the datetime library.
# @return datelist 

def load_restart_dates(file,datelist):
    from scipy.io import netcdf
    import math
    import datetime

    fp = netcdf.netcdf_file(file, 'r', mmap=False)
 
    # Retrieve the timing information from the file
    # ymd should be an integer
    ymd   = fp.variables['timemgr_rst_curr_ymd'].data
    dom   = int(ymd - 100*math.floor(ymd/100.0))
    moy   = int(math.floor(ymd/100) - 100*math.floor(ymd/10000))
    yr    = int(math.floor(ymd/10000))
    sod   = int(fp.variables['timemgr_rst_curr_tod'].data)
    datelist.append(datetime.datetime(yr,moy,dom,sod))
    fp.close()   # Close the netcdf file
    return(datelist)


# =======================================================================================

## compare the list of all known SOIs against the first restart file (which should have 
# the same grid structure as all restarts). If a grid-cell exists that is within a 
# distance tolerance, then the site is categorized as available, and the index of 
# the grid-cell is saved.
# @param restart_file  the name of the first test restart file
# @param sites_known   a list of class sitetype 
# @param geo_thresh    the maximum acceptable distance 
# between the SOI and gridcell (degrees)
#
def filter_rest_hist_sites(file,filetype,sites_known,geo_thresh):

    from scipy.io import netcdf
    import numpy as np

    ## netcdf object, think this is some type of pointer map to the file contents
    fp = netcdf.netcdf_file(file, 'r', mmap=False)

    if(filetype=='restart'):
    ## longitude data from the 1st restart file
        lons = fp.variables['grid1d_lon'].data
    ## latitude data from the 1st restart file
        lats = fp.variables['grid1d_lat'].data

    elif(filetype=='history'):

        ## Test if this is a 2d or irregular
        if ('lat' in fp.dimensions):
            hgrid = 2
        else:
            hgrid = 1

        ## longitude data from the 1st history file
        lons = fp.variables['lon'].data
        
        ## latitude data from the 1st history file
        lats = fp.variables['lat'].data

    else:
        print('Bad file type passed to the site filter')
        sys.exit(2)


    ## list of sites that have confirmed coverage within the file
    #  this is a list of class sitetype
    sites_avail = []
    
    for site in sites_known:

        ## the gridcell index with the smallest squared difference
        # in geoposition with the current SOI
        
        if(filetype=='restart' or (hgrid==1) ):
            igr = np.argmin( (lons-site.lon)**2.0 + (lats-site.lat)**2.0 )

            if ( (((lons[igr]-site.lon)**2.0 + (lats[igr]-site.lat)**2.0)**0.5) < geo_thresh):
                print('Site: '+site.name+' was located in the restart grid')
                sites_avail.append( sitetype(site.name,site.lat,site.lon)   )
                if(filetype=='restart'):
                    sites_avail[-1].igr = igr
                else:
                    sites_avail[-1].igh = igr
                    sites_avail[-1].hgrid = hgrid
            else:
                print('Site: '+site.name+' was NOT located in the restart grid')
                print('glon: '+str(lons[igr])+' slon: '+str(site.lon)+' glat: '+str(lats[igr])+' slat: '+str(site.lat))

        elif(filetype=='history' and (hgrid==2) ):
            ilat = np.argmin( (lats-site.lat)**2.0 )
            ilon = np.argmin( (lons-site.lon)**2.0 )
            if ( (((lons[ilon]-site.lon)**2.0 + (lats[ilat]-site.lat)**2.0)**0.5) < geo_thresh):
                print('Site: '+site.name+' was located in the restart grid')
                sites_avail.append( sitetype(site.name,site.lat,site.lon)   )
                sites_avail[-1].ilath=ilat
                sites_avail[-1].ilonh=ilon
                sites_avail[-1].hgrid=hgrid
            else:
                print('Site: '+site.name+' was NOT located in the restart grid')
                print('glon: '+str(lons[ilon])+' slon: '+str(site.lon)+' glat: '+str(lats[ilat])+' slat: '+str(site.lat))
        else:
            print('Uknown type (restart/history) and grid combination, exiting')
            exit()
            

        
    fp.close()
    return sites_avail




# ========================================================================================

## interp_args processes the arguments passed to the python script at executions.
#  The options are parsed and logical checks and comparisons are made.
# @param argv
# @return plotmode
# @return regressionmode
# @return restartmode
# @return test_r_prefix
# @return base_r_prefix
# @return test_h_prefix
# @return base_h_prefix
# @return test_name
# @return base_name

def interp_args(argv):

    argv.pop(0)  # The script itself is the first argument, forget it

    ## Binary flag that turns on and off plotting
    plotmode = False
    ## Binary flag that turns on and off regression tests against a baseline run
    regressionmode = False
    ## Binary flag that turns on and off the use of restart files for evaluation and comparison
    restartmode = False
    ## File path to the directory containing restart files from that baseline simulation
    base_r_prefix  = ''
    ## File path to the directory containing restart files from the test simulation
    test_r_prefix  = ''
    ## File path to the directory containin the history files from the test simulation
    test_h_prefix = ''
    ## File path to the directory containin the history files from the base simulation
    base_h_prefix = ''
    ## Name of the evaluation being performed, this is non-optional
    eval_id = ''
    ## Name for plot labeling of the test case
    test_name = 'test'
    ## Name for plot labeling of the base case
    base_name = 'base'

    try:
        opts, args = getopt.getopt(argv, 'h',["help","plotmode","regressmode",     \
                                              "restartmode","eval-id=",            \
                                              "test-rest-pref=","base-rest-pref=", \
                                              "test-hist-pref=","base-hist-pref=", \
                                              "test-name=","base-name="])

    except getopt.GetoptError as err:
        print('Argument error, see usage')
        usage()
        sys.exit(2)
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif o in ("--plotmode"):
            plotmode = True
        elif o in ("--regressmode"):
            regressionmode = True
        elif o in ("--restartmode"):
            restartmode = True
        elif o in ("--eval-id"):
            eval_id = a
        elif o in ("--test-rest-pref"):
            test_r_prefix = a
        elif o in ("--base-rest-pref"):
            base_r_prefix = a
        elif o in ("--test-hist-pref"):
            test_h_prefix = a
        elif o in ("--base-hist-pref"):
            base_h_prefix = a
        elif o in ("--test-name"):
            test_name = a
        elif o in ("--base-name"):
            base_name = a
        else:
            assert False, "unhandled option"

    if(plotmode):
        print('Plotting is ON')
    else:
        print('Plotting is OFF')

    if(restartmode):
        if(test_r_prefix==''):
            print('You specified restart mode, you must also specify')
            print('the directory prefix to those files')
            usage()
            sys.exit(2)
        
    if(regressionmode):
        print('Regression Testing is ON')
        if(base_r_prefix==''):
            print('In a regression style comparison, you must specify a')
            print('path to baseline restarts. See usage:')
            usage()
            sys.exit(2)
        if(restartmode):
            if(base_r_prefix==''):
                print('When regression and restart modes are active,')
                print('you must provide a path to baseline restart files')
                usage()
                sys.exit(2)
    else:
        print('Regression Testing is OFF')

    if(test_h_prefix==''):
        print('A path to history files is required input, see usage:')
        usage()
        sys.exit(2)

    if(eval_id==''):
        print('You must provide a name/id for this evaluation, see --eval-id:')
        usage()
        sys.exit(2)

    # Remove Whitespace in eval_id string
    eval_id.replace(" ","")
        

    return (plotmode, regressionmode, restartmode, eval_id, test_r_prefix, \
                base_r_prefix, test_h_prefix, base_h_prefix, test_name, base_name)




# ========================================================================================
# ========================================================================================
#                                        Main
# ========================================================================================
# ========================================================================================

def main(argv):

    # Interpret the arguments to the script
    plotmode, regressionmode, restartmode, eval_id, test_r_prefix, \
        base_r_prefix, test_h_prefix, base_h_prefix, \
        test_name, base_name = interp_args(argv)
    
    
    ## list of restart file names for the test-case 
    test_restart_list = []
    
    ## list of restart file names of the base-case
    base_restart_list = []
    
    ## list of history file names for the test-case 
    test_h0_list = [] 
    
    ## list of history file names for the base-case
    base_h0_list = []
    
    ## number of test-case restart files
    n_r_files = int(-9)
    
    ## number of base-case restart files
    n_rb_files = int(-9)
    
    ## number simulation types being used for restart analysis
    # (1 includes only a test, 2 also includes base simulations)
    n_rtypes = 0
    
    ## number of simulation types being used for history analysis 
    # (1 is test, 2 includes baseline)
    n_htypes = 1
    
    ## list of date objects for restarts, date objects contain the file timing information
    # note that only one list (test and base) is ultimately needed, they must be the same
    # list or the script aborts in error
    restart_datelist = []
    
    ## initialize empty list of date objects in the base simulation restart files
    restart_datelist_b = []

    ## list of date objects in the test simulation output, this is non-optional
    test_h0_list = getnclist(test_h_prefix,'h0')
    n_h_files = test_h0_list.__len__()
    if(n_h_files <= 0):
        print('acre could not find any h0 history files, this is non-optional')
        sys.exit(2)

    if(regressionmode):
         base_h0_list = getnclist(base_h_prefix,'h0')
         n_hb_files = base_h0_list.__len__()
         n_htypes = 2
         if(n_h_files != n_hb_files):
             print('Your test and base have different numbers of h0 files')
             sys.exit(2)

    if(restartmode):
        ## list of restart file names for the test-case 
        test_restart_list = getnclist(test_r_prefix,'restart')
        ## number of test-case restart files
        n_r_files = test_restart_list.__len__()
        n_rtypes  = 1
        if(n_r_files <= 0):
            print('acre could not find any restart files, and you optioned for restart analysis')
            sys.exit(2)
        
        if(regressionmode):
            base_restart_list = getnclist(base_r_prefix,'restart')
            n_rb_files = base_restart_list.__len__()
            n_rtypes = 2

            if(n_r_files != n_rb_files):
                print('Your test and base have different numbers of restart files')
                sys.exit(2)


    # ========================================================================================
    # Load up the sites of interest
    # ========================================================================================

    ## a list of all known sites that should be compared to model output grids
    #  this is a list of class sitetype

    sites_known = load_sites(xmlfile,sitetype)
    sites_avail = filter_rest_hist_sites(test_h0_list[0],'history', \
                                              sites_known,geo_thresh)

    # Check to see if the same sites that are available in the restart
    # are also available in the history, if not, abort because that is weird

    if(restartmode):
        
        rsites_avail = filter_rest_hist_sites(test_restart_list[0],'restart', \
                                         sites_known,geo_thresh)
        # Check to see if the valid sites in the history match that of the restart
        # If it passes this step, only one site list is necessary
        if(rsites_avail.__len__() != sites_avail.__len__()):
            print('It seems odd that the hsitory and restart files do not contain')
            print('the same number of sites')
            sys.exit(2)
        else:
            for i in range(rsites_avail.__len__()):
                if(sites_avail[i].name != rsites_avail[i].name):
                    print('site incompatibility in history and restart files')
                    sys.exit(2)
                else:
                    # Copy over the restart grid index
                    sites_avail[i].igr = rsites_avail[i].igr
        rsites_avail = None

    if(len(sites_avail) == 0):
        print('')
        print('No sites of interest were located in your output files.')
        print('Exiting.')
        sys.exit(2)


    # ========================================================================================
    # Load up the output-file time-stamping
    # this is done before the processing of diagnostic variables because it
    # is necessary for array initialization and evaluating what types of output
    # diagnostics can be processed
    # ========================================================================================

    print('Processing File Timing information\n')
    
    if(restartmode):
        for file in test_restart_list:
            # Append the datelist
            restart_datelist = load_restart_dates(file,restart_datelist)

        if(regressionmode):
            for file in base_restart_list:
                restart_datelist_b = load_restart_dates(file,restart_datelist_b)
                if(restart_datelist != restart_datelist_b):
                    print('Baseline restart file-times inconsistent with test')
                    sys.exit(2)
                else:
                    restart_datelist_b = None
 

    # ========================================================================================
    # Evaluate the first history file for dimension information and save it as the hist_dims 
    # class.  Restart analysis is more straightforward, this step is not needed.
    # Openning all 
    # ========================================================================================

    hdims = hutils.hist_dims(test_h0_list[0])
    hdims.timing([test_h0_list[0],test_h0_list[1],test_h0_list[-1]])

    hvarlist = hutils.define_histvars(xmlfile,hdims,test_h0_list[0],n_htypes,test_name,base_name)


    # Initialize the plot file
    # ========================================================================================
    plotfile_name = eval_id+"_plots.pdf"
    if(plotmode):
        pdf = PdfPages(plotfile_name)

    # Initialize the summary output table
    summary_table_name = eval_id+"_summary_table.txt"
    summary_table = open(summary_table_name ,"w")
    tutils.table_header(summary_table)
   

    # ========================================================================================
    # Initialize the history (and if optioned, history) variable class, 
    # then loop through the files and diagnose the critical variables
    # ========================================================================================


    for site in sites_avail:

        start_time = time.time()
        print('Processing Diagnostics for: '+site.name)

        print('')
        print('Loading test case h0 files at site: '+site.name)
        
        
        scr = hutils.scratch_space(test_h0_list[0])
        
        for file in test_h0_list:
            hutils.load_history(file,site,hvarlist,0,scr,hdims)

            
        if(regressionmode):
            for file in base_h0_list:
                hutils.load_history(file,site,hvarlist,1,scr,hdims)

        for hvar in hvarlist:
            hvar.normalize_diagnostics()
            

        # Initialize the restart variables.  This is a list of class hist_vars
        # found in the acre_history_utils module (hutils)

        if(restartmode):
            # Initialize (or re-initialize) the rvars class
            rvar = rutils.rvars(n_rtypes,n_r_files,test_name,base_name)
    
            # Loop through the restart files and populate time-series diagnostics
            print('')
            print('Loading test case restart files at site: '+site.name)
            for file in test_restart_list:
                rvar.load_restart(file,site.igr,0)

            if(regressionmode):
                print('')
                print('Loading base case restart files at site: '+site.name)
                for file in base_restart_list:
                    rvar.load_restart(file,site.igr,1)


        print('Site {} Diagnosed in {} seconds'.format(site.name,(time.time()-start_time)))
        if(plotmode):

            if(hdims.hperiod<=(24*32)):
                putils.multipanel_histplot(site,hvarlist,"MMV",n_htypes,pdf)
            else:
                print('Omitting plots of monthly means, history frequency is too coarse')

            if(hdims.hperiod<12):
                putils.multipanel_histplot(site,hvarlist,"DMV",n_htypes,pdf)
            else:
                print('Omitting plots of diurnal means, history frequency is too coarse')

            if(hdims.nyears>1):
                
                putils.multipanel_histplot(site,hvarlist,"AMV",n_htypes,pdf)
            else:
                print('Omitting plots of annual trends, Less than two years of data provided')

            if(restartmode):
                ## Make some restart analysis plots
                putils.quadpanel_restplots(site,rvar,restart_datelist,n_rtypes,pdf)
        
            print('Close the current figures to advance')
            plt.show()

        else:
            print('Plotting was turned off, no plots to show')
   


        tutils.site_header(summary_table,site)

        for hvar in hvarlist:
            tutils.site_var_write_line(summary_table,hvar,0)
            if(regressionmode):
                tutils.site_var_write_line(summary_table,hvar,1)


    summary_table.close()

    print('ACRE complete.') 
    if(plotmode):
        pdf.close()
        print('Two files have been generated:')
        print(summary_table_name)
        print(plotfile_name)
    else:
        print('One file has been generated:')
        print(summary_table_name)
    exit(0)




# =======================================================================================
# This is the actual call to main
   
if __name__ == "__main__":
    main(sys.argv)








