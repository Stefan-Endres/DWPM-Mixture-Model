#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to save and load dictionary containers with one dimensional arrays
to and from .csv files, the first row of the .csv file is the keys used in the
dictionary.
"""
__author__ = 'Stefan Endres'
__email__ = 'Stefan.C.Endres@gmail.com'
__status__ = 'Development'

#%%CSV Handling functions      
def save_dict_as_csv(Dict, Path, Order=None):
    """
    Saves a dictionary object with keys as column headers, only the first 
    container value will be saved (example an array or list will be saved,
    but a second list inside the first dictionary list will not.).
       
    Parameters
    ----------
    Dict : Dictionary container of lists containing int, floats and strings.
           Input data to be saved.
           
    Path : String.
           Path and/or filename to be saved.
           
    Order : List containing strings, optional.
            The from left to right.
    
    Examples
    --------
    >>> Save_Dict_as_csv(Store, 'Results/Store.csv')
    >>> Save_Dict_as_csv(Store, 'Results/Store.csv', Order=['x1','x2','y1'])

    Dependencies
    ------------
    csv, numpy
    
    !TO DO: TEST NON-ORDERED NON-LIST SAVING (TypeError exception if not Order)
    """
    
    import csv  
    import numpy
    with open(Path, 'wb') as outfile:
        writer = csv.writer(outfile)
        DictVals = 0.0                        
        for key, value in Dict.items():         # Find largest dictionary value
            if numpy.size(value) > DictVals:
                DictVals = numpy.size(value)

        rows = [[]]
        for i in range(DictVals+1):             # Value Rows
            rows.append([])  
            
        if Order is None: # If no order is specified, save in any order
            for key, value in Dict.items():
                rows[0].append(key)    # Save headers
                for j in range(0,DictVals): # Save row of current header value
                    try:
                        if not value[j] == '':
                            rows[j+1].append(value[j])
                        else:
                            rows[j+1].append('\'\'')
                    except NameError: # If entry is a string not a float
                        rows[j+1].append('\''+value[j]+'\'')
                    except IndexError: # If no entry found append empty value
                        rows[j+1].append('\'\'')
                        
                    except TypeError:# Issue with dictionary entry
                        if type(value) is float or type(value) is str \
                        or type(value) is int: # If only 1 entry found
                            if j == 0: # Append val to first entry
                                rows[j+1].append(value)
                            else: # Fill column entries with empty values
                                rows[j+1].append('\'\'')
                        else: # Raise error if unknown object in dictionary
                            raise IOError('Unknown dictionary format, a ' 
                                           'dictionary entry contains an '
                                           'object that is not a list of float'
                                           ' str or int')
        else: # If heading order is specified save in correct order
            for n in Order:
                rows[0].append(n)    # Save headers
                for j in range(0,DictVals):# Save row of current header value
                    try:  
                        if not Dict[n][j] == '':
                            rows[j+1].append(Dict[n][j])
                        else: # If empty entry print empty value
                            rows[j+1].append('\'\'')
                    except NameError:  # If entry is a string not a float
                        rows[j+1].append('\''+Dict[n][j]+'\'')
                    except IndexError:# If no entry found append empty value
                        rows[j+1].append('\'\'')
                    except TypeError:# Issue with dictionary entry
                        if type(Dict[n]) is float or type(Dict[n]) is str \
                        or type(Dict[n]) is int: # If only 1 entry found
                            if j == 0: # Append val to first entry
                                rows[j+1].append(Dict[n])
                            else: # Fill column entries with empty values
                                rows[j+1].append('\'\'')
                        else: # Raise error if unknown object in dictionary
                            raise IOError('Unknown dictionary format, a ' 
                                           'dictionary entry contains an '
                                           'object that is not a list of float'
                                           ' str or int')
            
        for i in range(DictVals+1): # Write all rows
            writer.writerow(rows[i])

    return

#%%
def dict_conv_numstrings(Dict, numberset = set('0123456789e-'),
                               letterset = set('abcdfghijklmnopqrstuvwxyz')):    
    """
    Convert all strings in dictionary containing only numbers to floats.
    
    Parameters
    ----------
    Dict : Dictionary container to be converted.
    
    numberset : set, optional.
                Defines numbers to be converted. Any character defined in this
                set will be considered a number to be converted.
    
    letterset : set, optional. 
                Defines non-number strings, if any character contained in this 
                set is contained in the string it will not be converted to a 
                number. If special non-number characters are used they should 
                be added to this set to avoid a raising a value error:
                            ValueError: invalid literal for float()
                
                Note: By default 'e' is not included in this set so exponential
                numbers such as 1.3806488e-23 will be converted to numbers.
    
    Examples
    --------
    >>> A = {'x': ['123','abc','12c', '12a', '3.1e-05']}
    >>> A
    {'x': ['123', 'abc', '12c', '12a', '3.1e-05']}
    >>> dict_conv_numstrings(A)
    {'x': [123.0, 'abc', '12c', '12a', 3.1e-05]}
   
    >>> A = {'x': ['123','abc','12c', '12ä']}
    >>> dict_conv_numstrings(A, 
                             letterset = set('äabcdefghijklmnopqrstuvwxyz')) 
    {'x': [123.0, 'abc', '12c', '12\xc3\xa4']}
    """

    for key, value in Dict.items():
        for i,n in enumerate(Dict[key]):
            try:
                if any((c in numberset) for c in n):
                    if not any((c in letterset) for c in n):
                        Dict[key][i] = float(Dict[key][i]) 
        
            except TypeError:
                pass      
            
    return Dict
            
#%%            
def load_csv_as_dict(Path, start_row=1, skip_rows=None, conv_numstrings=True, \
                     empty_entries=True):
    """
    Returns a dictionary object with keys as column headers.
    
    Parameters
    ----------
    Path : String.
           Path and/or filename to be loaded.
           
    start_row : Integer, optional.
                Skips to this row number in the csv file to be read. 
                Headings in this row will be used as keys
    
    skip_rows : List with integer entries, optional
                This will skip any row specified in the input list.
               
    conv_numstrings : Boolean, optional
                      Convert all strings in dictionary containing only 
                      numbers to floats.     
                      
    empty_entries : Boolean, optional
                    
                    
    Examples
    --------
    >>> A = load_csv_as_dict('Data/Store.csv')
    >>> A
    {' description y': ['y', 11.0, 12.0, 13.0], ' description z': ['z', 100.0],
    'description x': ['x', 1.0, 'abc', 3.1e-05]}
    
    >>> A = load_csv_as_dict('Data/Store.csv', start_row=2)
    >>> A
    {'y': [11.0, 12.0, 13.0], 'x': [1.0, 'abc', 3.1e-05], 'z': [100.0]}
    
    >>> A = load_csv_as_dict('Data/Store.csv', start_row=2, skip_rows=[4])
    >>> A
    {'y': [11.0, 13.0], 'x': [1.0, 3.1e-05], 'z': [100.0]}
    
    >>> A = load_csv_as_dict('Data/Store.csv', start_row=2, skip_rows=[4],
                             conv_numstrings=False)
    >>> A
    {'y': ['11', '13'], 'x': ['1', '3.1e-05'], 'z': ['100']}
    
    >>> A = load_csv_as_dict('Data/Store.csv', start_row=2, skip_rows=[4],
                             empty_entries=False)
    >>> A
    {'y': [11.0, 13.0], 'x': [1.0, 3.1e-05], 'z': [100.0, '']}
    
    Dependencies
    ------------
    csv, self.dict_conv_numstrings
    """
    import csv
    csv_read = csv.DictReader(open(Path))
    i = 1 
    # Shift skip_rows back 1 to correspond with correct csv_read entry indices
    if skip_rows is not None:
        skip_rows = [skip_rows[k] - 1 for k in range(len(skip_rows))]
        

    Dict = {} # Empty dictionary to return
    for row in csv_read:
        if start_row == 1: # Build dictionary without skipping rows
            if skip_rows is not None:
                if i in skip_rows: # Skip rows specified in input list
                    i += 1    
                    continue
                
            for key, value in row.items():
                Dict.setdefault(key, []).append(value)
                
            i += 1  
            continue
        
        if i  < start_row - 1: # Skip row
            i += 1  
            continue
        
        if i == start_row - 1: # Redefine row keys
            keys = []
            for column, value in row.iteritems():
                keys.append(value)
            i += 1  
            continue 
        
        if i > start_row - 1: # Build new dictionary from column entries
            if skip_rows is not None:
                if i in skip_rows: # Skip rows specified in input list
                    i += 1    
                    continue
                    
            for key, (column, value) in zip(keys, row.iteritems()):
                Dict.setdefault(key, []).append(value)
            
            i += 1      
         
    # Convert number strings to floats
    if conv_numstrings: dict_conv_numstrings(Dict)
    
    # Remove empty entires
    if empty_entries:
        for key, value in Dict.items(): # Filter out empty '' values
            if not value.__class__ == float:
                Dict[key] = list(filter(lambda a: a != '', value))
    
    return Dict

#%% 
def readvarfile(filename):
    import csv
    with open(filename) as f:
        reader = csv.reader(f)
        # The first line contains headings
        reader.next()
        # Now extract the data
        names, valuestrings, units, descriptions = zip(*list(reader))
    # Convert the string representation of numbers to proper floats
    values = [float(s) for s in valuestrings]
    return names, values, units, descriptions

# Read parameter file
# parnames, parvalues, parunits, pardescriptions = readvarfile('parameters.csv')
    