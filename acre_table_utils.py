import numpy as np
#import code  # For development: code.interact(local=locals())
import sys

max_var_name_len = 24
max_site_name_len = 112
max_var_comp_name_len = 16

def table_header(summary_table):
    summary_table.write("Notes:")
    summary_table.write("MEAN: Mean of yearly means (if available)")
    summary_table.write("      Mean of monthly means (if next available)")
    summary_table.write("      Mean of diurnal means (if next available)")
    summary_table.write("MEAN will accept sparse data (ie ignores nans)")
    summary_table.write("AMIN: Minimum of yearly means (ignores nans)")
    summary_talbe.write("AMAX: Maximum of yearly means (ignores nans)")
    summary_table.write("MMIN: Minimum of monthly means (nans disqualify)")
    summary_talbe.write("MMAX: Maximum of monthly means (nans disqualify)")
    summary_table.write("DMIN: Minimum of diurnal means (nans disqualify)")
    summary_table.write("DMAX: Maximum of diurnal means (nans disqualify)")
    summary_table.write("================================================")


def site_header(summary_table,site):

    # Generate tabular output for this site
    # =============================================================
    summary_table.write("{:16}".format(site.name))
    summary_table.write("\n")
    summary_table.write("{:16}".format(""))
    summary_table.write("{:>24}".format("MEAN"))
    summary_table.write("{:>24}".format("AMIN"))
    summary_table.write("{:>24}".format("AMAX"))
    summary_table.write("{:>24}".format("MMIN"))
    summary_table.write("{:>24}".format("MMAX"))
    summary_table.write("{:>24}".format("DMIN"))
    summary_table.write("{:>24}".format("DMAX"))
    summary_table.write("\n")


def site_var_write_line(summary_table,hvar,htype):


    if(htype==0):
        name_tag = "(t)"
    else:
        name_tag = "(b)"

    summary_table.write("{:16}".format(hvar.name+' '+name_tag))

    # Mean (Try annual first, then monthly, then daily)
    if(hvar.amv):
        num_str = "{:1.3e}".format(np.nanmean(hvar.amv_ar[:,htype]))
        tab_str = "{:>24}".format(num_str+" "+hvar.amv_unit.strip())
        summary_table.write(tab_str)
    elif(hvar.mmv):
        num_str = "{:1.3e}".format(np.nanmean(hvar.mmv_ar[:,htype]))
        tab_str = "{:>24}".format(num_str+" "+hvar.mmv_unit.strip())
        summary_table.write(tab_str)
    elif(hvar.mmv):
        num_str = "{:1.3e}".format(np.nanmean(hvar.dmv_ar[:,htype]))
        tab_str = "{:>24}".format(num_str+" "+hvar.dmv_unit.strip())
        summary_table.write(tab_str)
    else:
        print('No history time-scales selected?')
        print('Exiting.')
        exit(0)

    # Annual Max/Min
    if(hvar.amv):
        num_str = "{:1.3e}".format(np.nanmin(hvar.amv_ar[:,htype]))
        tab_str = "{:>24}".format(num_str+" "+hvar.amv_unit.strip())
        summary_table.write(tab_str)

        num_str = "{:1.3e}".format(np.nanmax(hvar.amv_ar[:,htype]))
        tab_str = "{:>24}".format(num_str+" "+hvar.amv_unit.strip())
        summary_table.write(tab_str)
    else:
        summary_table.write("{:>24}".format("NA"))
        summary_table.write("{:>24}".format("NA"))

    # Diurnal Max/Min
    if(hvar.dmv):
        num_str = "{:1.3e}".format(np.min(hvar.dmv_ar[:,htype]))
        tab_str = "{:>24}".format(num_str+" "+hvar.dmv_unit.strip())
        summary_table.write(tab_str)

        num_str = "{:1.3e}".format(np.max(hvar.dmv_ar[:,htype]))
        tab_str = "{:>24}".format(num_str+" "+hvar.dmv_unit.strip())
        summary_table.write(tab_str)
    else:
        summary_table.write("{:>24}".format("NA"))
        summary_table.write("{:>24}".format("NA"))

    # Monthly Max/Min
    if(hvar.mmv):
        num_str = "{:1.3e}".format(np.min(hvar.mmv_ar[:,htype]))
        tab_str = "{:>24}".format(num_str+" "+hvar.mmv_unit.strip())
        summary_table.write(tab_str)

        num_str = "{:1.3e}".format(np.max(hvar.mmv_ar[:,htype]))
        tab_str = "{:>24}".format(num_str+" "+hvar.mmv_unit.strip())
        summary_table.write(tab_str)
    else:
        summary_table.write("{:>24}".format("NA"))
        summary_table.write("{:>24}".format("NA"))

    summary_table.write("\n")





def str_pad(raw_str,str_len):

    npads = str_len - len(raw_str)
    if(npads>0):
        for i in range(0,npads):
            raw_str += ' '
    else:
        print('A site or variable has a name that is too long:')
        print('String Name:',+raw_str)
        print('Max Lenght:',+(str_len-1))

    return(raw_str)
