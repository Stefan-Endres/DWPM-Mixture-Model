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

    def load(self, directory, Comps):
        from csvDict import load_csv_as_dict
        # Find file name path for specified components
        filename = '_'.join(Comps)

        # FIXME: Use os.path.join
        ldstr = directory + filename + '.csv'

        # FIXME: Use os.path.exist
        try: # Load data from file path
            Data  = load_csv_as_dict(ldstr)
            self.VLE = Data
        except IOError: # Raise error if not found
            raise IOError('Phase data for system "{}" not found'.format(filename))


    def load_VLE(self, Comps):
        """
        Returns multicomponent VLE data from specified BINARY components.
        """
        self.load('Data/Binary_VLE/', Comps)

    def load_E(self, Comps):
        """
        Returns multicomponent equilibrium data from specified components.
        """
        self.load('Data/nComp_E/', Comps)


    def test_internal(self): # TEST; DELETE
         self.test_int() # TEST; DELETE

    def test_int(self): # TEST; DELETE
         print 'Test Succesful!'

if __name__ == '__main__':
    data.test_internal()
    data = ImportData()
    data.load_pure_data(Compounds) # Call ex. component 1: data.c[0]
    if len(Compounds) > 1:
        data.load_VLE(Compounds)       # Call ex. component 1: data.c[0]
