
class bmark_vars:

    def __init__(self,symbol,name,dimclass,unit):

        self.symbol = symbol
        self.name   = name
        self.scv    = False    # Size-class variable
        self.sv     = False    # Scalar variable
        self.unit   = unit
        
        self.hfile_id = hfile_id
        

        # time dimension
        if('time' in dimnames):
            self.timedim = [id for id in range(dimnames.__len__())\
                            if dimnames[id] in 'time'][0]
        else:
            print('History Variable: '+name+'has no time dimension?')

        self.dimclass = dimclass

        # This is a size-dimensioned variable
        if(self.dimclass == 'lnd_scpf'):

            # Initialize size-class array
            self.scv = True
            self.scv_mod_ar = np.zeros((hdims.nsize,n_htypes))
            self.scv_mod_n  = np.zeros((hdims.nsize,n_htypes))
            self.scv_obs_ar = np.zeros((hdims.nsize,1))
            
            self.scv_x      = np.zeros((hdims.nsize,1))
            self.scv_x_unit = 'DBH [cm]'

            # Initialize scalar array
            self.sv        = True
            self.sv_mod_ar = np.zeros((n_htypes,1))
            self.sv_mod_n  = np.zeros((n_htypes,1))
            self.sv_obs_ar = np.zeros((1,1))


        if(self.dimclass == 'lnd'):

            # Initialize scalar array
            self.sv        = True
            self.sv_mod_ar = np.zeros((n_htypes,1))
            self.sv_mod_n  = np.zeros((n_htypes,1))
            self.sv_obs_ar = np.zeros((1,1))
