#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Data read/write.

ex.
>>> dataset1.load_VLE(['acetone','water'])
>>> dataset2.load_VLE(['benzene','cyclohexane'])
"""
#%%CSV Handling functions
# [moved to csvDict.py]
   
#%% Load data classes
class ImportData:
    """Load data of coupounds specified with Data.load_pure_data()"""
    def __init__(self):
        #self.name = name
        self.c = []    # creates a new empty list for components 
    
    def load_pure_data(self, Comps):
        from csvDict import load_csv_as_dict
        """
        Returns the pure component VLE data for specified components.
        """
        for i,comp in enumerate(Comps):
            ldstr           = 'Data/Pure_Component/{}.csv'.format(comp)
            try:
                Data = load_csv_as_dict(ldstr)
                Data['name'] = [comp,]
                self.c.append(Data)

            except IOError:  # Raise error if not found
                raise IOError('Data for specified component '  \
                               +'"{}" not found'.format(comp))
                
    def load_VLE(self, Comps):
        """
        Returns multicomponent VLE data from specified components.
        """
        from csvDict import load_csv_as_dict
        # Find file name path for specified components
        filename = '' 
        for i in range(len(Comps)-1):
            filename += '{}_'.format(Comps[i])

        filename +='{}'.format(Comps[len(Comps)-1])
        
        ldstr = 'Data/Binary_VLE/'+filename+'.csv'

        try: # Load data from file path
            Data  = load_csv_as_dict(ldstr)
            self.VLE = Data
        except IOError: # Raise error if not found
            raise IOError('Phase data for system "{}" not found'.format(filename))
     
     
    def test_internal(self): # TEST; DELETE
         self.test_int() # TEST; DELETE
         
    def test_int(self): # TEST; DELETE
         print 'Test Succesful!'

#Comps = ['acetone','water']
#data = import_data()
#data.load_pure_data(Comps)
#data.load_VLE(Comps)
#data.c[1]
#data.VLE
#%%
#def import_data():
#    data = import_data()
#    data.load_pure_data(Compounds) # Call ex. component 1: data.c[0]
#    data.load_VLE(Compounds)       # Call ex. component 1: data.c[0]

 
#%%
if __name__ == '__main__':
    data.test_internal() # TEST; DELETE
    data = ImportData()
    data.load_pure_data(Compounds) # Call ex. component 1: data.c[0]
    if len(Compounds) > 1:
        data.load_VLE(Compounds)       # Call ex. component 1: data.c[0]
    
    
    
    #save_dict_as_csv(data.c[0],'test.csv')
              
    #Order = ['T (K)', 'P (Pa)', 'T_c (K)', 'P_c (Pa)', 'V_c (m3 mol-1)', 'Z_c',
              # 'R (m3 Pa K−1 mol−1)' ,'w' ,'a_c', 'b_c',
    #         'm (Adachi-Lu)', 'm (Soave)','virialT','virialB']

    #save_dict_as_csv(data.c[0],'test2.csv', Order)
