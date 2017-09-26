import matplotlib.pyplot as plt
import numpy as np
#import code  # For development: code.interact(local=locals())
import sys
import matplotlib.dates as mdates


## Define the plotstruct class structure
#  there are many ways to create multiplanel plots
# multipanel plots are used a lot here because we need to minimize
# all of the things that pop up on the user's screen while maximizing
# all of the cool stuff that is generated.
# @param n_hpan the number of panels in the horizontal
# @param n_vpan the number of panels in the vertical

class plotstruct:

    def __init__(self,n_hpan,n_vpan):
        
        ## the panel order index, on the first call
        # for a new position this will be incremented
        self.ipan = 0

        # 4 panel plot
        if(n_hpan==2 & n_vpan==2):
            ## number of panels in the horizontal
            self.n_hpan = n_hpan
            ## number of panels in the vertical
            self.n_vpan = n_vpan
            ## simulation type labeling for plots
#            self.str_types = ["test","base"]
            ## simulation type plot specifiers (line type, line color)
            self.rline_types = ["k-","b--"]
            self.moline_types = ["ko-","b^--"]
            ## depth of the header where the plot main title is
            self.head_depth = 0.04
            ## vertical offset for the tiling
            self.voff = 0.09
            ## horizontal offset for the tiling
            self.hoff = 0.04
            ## horizontal spacing between panels
            self.mh   = 0.08
            ## vertical spacing between panels
            self.mv   = 0.07
            ## horizontal position (left edge) of current panel
            self.hpos = -9
            ## vertical position (bottom edge) of current panel
            self.vpos = -9
            ## index of the horizontal panel
            self.ivpan = -9
            ## index of the vertical panel
            self.ihpan = -9
            ## vertical width of the current panel
            self.dv = (1.0-self.head_depth-self.voff-(self.n_vpan*self.mv))/self.n_vpan
            ## horizontal width of the current panel
            self.dh = (1.0-2*self.hoff-(self.n_hpan*self.mh))/self.n_hpan

        else:
            print("no support yet for {}x{} paneled plots".format([n_hpan,n_vpan]))
            sys.exit(2)



    ## updatepos updates the horizontal and vertical positions in the figure
    # based on the current panel index
    # of a new figure axes, given the previous order index of the figure, it also 
    # updates the order index.  Note that beginning indices are 1 (not 0)
    # @param ipan the previous order index of the panel
    # @return ipan the current order index
    def updatepos(self):
        
        self.ipan +=1

        if(self.n_hpan==2 & self.n_vpan==2):
            
            ## local variable denoting the vertical index (integer)
            self.ivpan = int(np.floor((self.ipan-1.0)/self.n_hpan))+1
            if(self.ivpan>self.n_vpan):
                print("you have more axes than the bounding dimensions of your figure")
                print("ipan: {}, ivpan: {}".format(self.ipan,self.ivpan))
                sys.exit(2)
                
            ## local variable denoting the horizontal index (integer)
            self.ihpan = self.ipan-(self.ivpan-1)*self.n_hpan
            self.ivpan = self.n_vpan-self.ivpan+1

            self.hpos = self.hoff + (self.ihpan-1)*(self.dh+self.mh) + self.mh
            self.vpos = self.voff + (self.ivpan-1)*(self.dv+self.mv)



#========================================================================================
## plot handle for multiplanel restart comparison.  These plots are all singular 
# time-series data.  The x-axis is always a date in years.
def quadpanel_restplots(site,rvar,restart_datelist,n_rtypes,pdf):
    
    fig0=plt.figure(figsize=[9,8])
    fig0.patch.set_facecolor('linen')
    fig0.suptitle("Site: {}".format(site.name),fontsize=14,horizontalalignment='center')

    ## number of date-ticks
    nyears = restart_datelist[-1].year-restart_datelist[0].year+1
    if(nyears>20):
        nticks = 6
    elif(nyears>6):
        nticks = 4
    else:
        nticks = 3
        
    ## year tick object that defines x-spacing
    yearticks = mdates.YearLocator(max([3,int(nyears/nticks)]), month=1, day=1 )
    
    ## initialize panel functions for a 2x2 multipanel plot
    pdim0 = plotstruct(2,2)
    
    # Structural biomass
    
    # move the panel to the first position
    pdim0.updatepos()
    a = plt.axes([pdim0.hpos,pdim0.vpos,pdim0.dh,pdim0.dv])
    
    for i in range(n_rtypes):
        plt.plot_date(restart_datelist,rvar.bdead_gc[i,:], \
                      "{}".format(pdim0.rline_types[i]), \
                      label="{}".format(rvar.test_name_str[i]))
    plt.ylabel('[Mg/ha]')
    if(pdim0.ivpan==1): 
        plt.xlabel('Year')
    else:
        a.xaxis.set_ticklabels([])
    plt.title("bdead, all live trees")
    plt.grid(True)
    plt.gca().xaxis.set_major_locator(yearticks)

    # Basal Area
    
    # move the panel to the next position
    pdim0.updatepos()
    a = plt.axes([pdim0.hpos,pdim0.vpos,pdim0.dh,pdim0.dv])
        
    for i in range(n_rtypes):
        plt.plot_date(restart_datelist,rvar.ba_gc[i,:], \
                      "{}".format(pdim0.rline_types[i]), \
                      label=rvar.test_name_str[i])
    plt.ylabel('[m2/ha]')
    if(pdim0.ivpan==1): 
        plt.xlabel('Year')
    else:
        a.xaxis.set_ticklabels([])
    plt.title('Basal Area of Live Trees')
    plt.grid(True)
    plt.gca().xaxis.set_major_locator(yearticks)
        

    # Recruit Structural Biomass (dbh<5cm)
        
    # move the panel to the next position
    pdim0.updatepos()
    a = plt.axes([pdim0.hpos,pdim0.vpos,pdim0.dh,pdim0.dv])

    for i in range(n_rtypes):
        plt.plot_date(restart_datelist,rvar.bdead_dbh5_gc[i,:], \
                      "{}".format(pdim0.rline_types[i]), \
                      label=rvar.test_name_str[i])
    if(pdim0.ivpan==1): 
        plt.xlabel('Year')
    else:
        a.xaxis.set_ticklabels([])
    plt.ylabel('[Mg/ha]')
    plt.title('bdead (dbh<5cm)')
    plt.grid(True)
    plt.gca().xaxis.set_major_locator(yearticks)

    # Recruit Structural Biomass (dbh<2cm)

    # move the panel to the next position
    pdim0.updatepos()
    a = plt.axes([pdim0.hpos,pdim0.vpos,pdim0.dh,pdim0.dv])
    
    for i in range(n_rtypes):
        plt.plot_date(restart_datelist,rvar.bdead_dbh2_gc[i,:], \
                      "{}".format(pdim0.rline_types[i]), \
                      label=rvar.test_name_str[i])
    if(pdim0.ivpan==1): 
        plt.xlabel('Year')
    else:
        a.xaxis.set_ticklabels([])
    plt.ylabel('[Mg/ha]')
    plt.title('bdead (dbh<2cm)')
    plt.grid(True)
    plt.gca().xaxis.set_major_locator(yearticks)
    
    # Add a legen to the last axis
    plt.legend(loc='upper left')

    # Add a new axis for the top label
    pdf.savefig(fig0)

    fig1=plt.figure(figsize=[9,8])
    fig1.patch.set_facecolor('linen')
    fig1.suptitle("Site: {}".format(site.name),fontsize=14,horizontalalignment='center')

    ## initialize panel functions for a 2x2 multipanel plot
    pdim1 = plotstruct(2,2)
    
    
    # Cohorts Per Patch
    
    # move the panel to the first position
    pdim1.updatepos()
    a = plt.axes([pdim1.hpos,pdim1.vpos,pdim1.dh,pdim1.dv])
    
    for i in range(n_rtypes):
        plt.plot_date(restart_datelist,rvar.nc_per_pa_gc[i,:], \
                      "{}".format(pdim1.rline_types[i]), \
                      label="{}".format(rvar.test_name_str[i]))
    plt.ylabel('[N]')
    if(pdim1.ivpan==1): 
        plt.xlabel('Year')
    else:
        a.xaxis.set_ticklabels([])
    plt.title("Mean Cohorts Per Patch")
    plt.grid(True)
    plt.gca().xaxis.set_major_locator(yearticks)

    # move the panel to the second position
    pdim1.updatepos()
    a = plt.axes([pdim1.hpos,pdim1.vpos,pdim1.dh,pdim1.dv])
    
    for i in range(n_rtypes):
        plt.plot_date(restart_datelist,rvar.npa_gc[i,:], \
                      "{}".format(pdim1.rline_types[i]), \
                      label="{}".format(rvar.test_name_str[i]))
    plt.ylabel('[N]')
    if(pdim1.ivpan==1): 
        plt.xlabel('Year')
    else:
        a.xaxis.set_ticklabels([])
    plt.title("Number of Patches")
    plt.grid(True)
    plt.gca().xaxis.set_major_locator(yearticks)

    

    # Show the plot
    #return()
    pdf.savefig(fig1)



def multipanel_histplot(site,hvarlist,atype,n_htypes):
    
    count=50
    for hvar in hvarlist:

        if(atype=='MMV'): 
            acheck = hvar.mmv
            unitstr = hvar.mmv_unit
            y_array = hvar.mmv_ar
            x_array = hvar.mmv_x
            marksize = 8
            xlimits = [0.5,12.5]
            x_label  = 'Month'
            title_string = 'Monthly Means: '+site.name
        elif(atype=='DMV'): 
            acheck = hvar.dmv
            unitstr = hvar.dmv_unit
            y_array = hvar.dmv_ar
            x_array = hvar.dmv_x
            marksize = 6
            xlimits = [0,23]
            x_label = 'Hour'
            title_string = 'Diurnal Means: '+site.name
        elif(atype=='AMV'): 
            acheck = hvar.amv
            unitstr = hvar.amv_unit
            y_array = hvar.amv_ar
            x_array = hvar.amv_x     
            marksize = 0
            xlimits = [hvar.amv_x[0,0],hvar.amv_x[-1,0]]
            x_label = 'Year'
            title_string = 'Annual Means: '+site.name
        else:
            print('You specified an unknown plot type for history analysis')
            sys.exit(2)

        if(acheck):
            
            count += 1
            if(count>4):

                figh = plt.figure(figsize=[9,8])
                figh.patch.set_facecolor('linen')
                figh.suptitle("{}".format(title_string), \
                          fontsize=14,horizontalalignment='center')
                # initialize panel functions for a 2x2 multipanel plot
                pdim0 = plotstruct(2,2)
                count = 1

            pdim0.updatepos()
            a = plt.axes([pdim0.hpos,pdim0.vpos,pdim0.dh,pdim0.dv])
            for i in range(n_htypes):
                plt.plot(x_array[:,i],y_array[:,i],
                         "{}".format(pdim0.moline_types[i]),markersize=marksize, \
                             label="{}".format(hvar.test_name_str[i]))


            plt.ylabel(unitstr)
            if(pdim0.ivpan==1): 
                plt.xlabel(x_label)
                # Add a legen to only one panel
                plt.legend(loc='upper left')
            else:
                a.xaxis.set_ticklabels([])
                
            a.set_xlim(xlimits)
            plt.title(hvar.name)
            plt.grid(True)
        

    # Show the plot
    return()
