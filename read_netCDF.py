"""
COAPS Simplified Daily Scatterometer Swath netCDF Read Tool

A tool to load data from a COAPS Simplified Daily Scatterometer Swath netCDF
file into memory. See the read function's documentations for more details.

Necessary Modules:
    Python version 2.7.15 or later
    netCDF4 version 1.2.9
    
Usage:
    Command line:
        python read_netCDF.py file_name
        or
        python3 read_netCDF.py file_name
    
    Imported:
        from read_netCDF import read
        data = read(file_name)
        
Authors:
    Arturo Valery, Jocelyn Elya
    
Contact:
    mdc@coaps.fsu.edu
"""


from __future__ import print_function
from netCDF4 import Dataset
import sys

__author__ = "Arturo Valery, Jocelyn Elya"
__credits__ = ["Arturo Valery", "Jocelyn Elya"]
__license__ = "MIT"
__maintainer__ = "Jocelyn Elya"
__email__ = "mdc@coaps.fsu.edu"

def read(file_name):
    """
    Read data from the COAPS Simplified Daily Scatterometer Swath netCDF file
    into a dictionary where the keys are variable names and the values are the
    variable's data.    
    
    Argument:
    file_name - string value containing the relative or absolute path of the
        netcdf file
    
    Returns:
    variable_dict - dictionary where the keys are variable names and the values
        are masked numpy arrays of the variable's data as read from the file
    
    Variables in the COAPS Simplified Daily Scatterometer Swath netCDF files:
        time       
                   A one dimensional LONG array containing the mean time of the 
                   observations. Given as the seconds since 1990-01-01 00:00:00.
                   Each value in the array corresponds to one row of
                   observations: all times in across track (perpendicular to the
                   satellite track) are equal. The index corresponds to the
                   along track time and position.
             
        lat        
                   A two dimensional FLOAT array containing the latitude of the 
                   wind vectors. The first index corresponds to the along track
                   position, and the second index corresponds to the cross
                   track position.
             
        lon         
                   A two dimensional FLOAT array containing the longitude of the 
                   wind vectors. The first index corresponds to the along track
                   position, and the second index corresponds to the cross
                   track position.
             
        eastward_wind        
                   Zonal equivalent neutral wind vector component. Positive
                   values indicate an eastward movement. The first index
                   corresponds to the along track position, and the second index
                   corresponds to the cross track position.
              
        northward_wind        
                   Meridional equivalent neutral wind vector component. Positive
                   values indicate northward movement. The first index
                   corresponds to the along track position and the second index
                   corresponds to the cross track position.
              
        wind_speed       
                   A two dimensional FLOAT array containing the equivalent
                   neutral wind speed. The first index corresponds to the track
                   position, and the second index corresponds to the cross
                   track position.
                   
        wind_to_direction        
                   A two dimensional FLOAT array containing the direction of the
                   wind vector (oceonographic convention). The first index
                   corresponds to the along track position, and the second index
                   corresponds to the cross track position.

        simplified_wvc_quality_flag        
                   A two dimensional INTEGER array with two values, 0 or 1.
                   This flag is a summary of the product developer's non-rain
                   flag, with 0 indicating good data, and 1 indicating a missing
                   flag value, highly suspect data, or data of poor quality. See
                   the metadata for this variable to see flags that are combined
                   to produce this summary flag.
              
        rain_impact (QSCAT and RapidSCAT only)
                   A two dimensional FLOAT array that attempts to indicate the
                   impact of rain on a zero to 100 scale. At this time, the
                   interpretation of these numbers is poorly described for Ku
                   band, and not provided for C band. We anticipate this
                   changing in the future; however, at this time JPL recommends
                   most data users to avoid using the rain_impact variable. The 
                   only exception to this rule is for users who are trying to
                   do very specific studies that are highly adversely impacted
                   by rain.
        
        distance_from_coast (QSCAT and RapidSCAT only)
                   A two dimensional FLOAT array containing the distance from
                   the nearest coastline in kilometers. The first index
                   corresponds to the along track position and the second index
                   corresponds to the cross track position.
        
        ice_prob (ASCAT A and ASCAT B only)
                  A two dimensional FLOAT array containing the ice probability.
                  The first index corresponds to the along track position and
                  the second index corresponds to the cross track position.          
    """
    
    ncdf_file = Dataset(file_name, 'r')
    
    source = str(ncdf_file.source)
    
    if ("QuikSCAT" in source) or ("Rapidscat" in source):
        grid_spacing = str(ncdf_file.cross_track_resolution)
    else:
        grid_spacing = str(ncdf_file.pixel_size_on_horizontal)


    ncdf_attributes = ncdf_file.ncattrs() # List of attributes
    ncdf_vars = ncdf_file.variables # List of variables
    
    # Will hold the variable name as the key and its data as the value
    variable_dict = dict() 
    
    print("\n")
    print("Source:", source)
    print("Grid Spacing:", grid_spacing)
    print("Variables:")
    for var in ncdf_vars:
        print("\t", var)
        print("\t", "-"*len(var))
        if var not in variable_dict:
            # Loads the NetCDF array into a numpy array
            variable_dict[var]= ncdf_file.variables[var][:]  
            print(variable_dict[var])
            print("\n")
    ncdf_file.close()
    return variable_dict



if __name__ == '__main__':
    args = sys.argv
    if len(args) != 2:
        print("ERROR")
        print("USAGE:")
        print("python read_netCDF.py file_name")
        print("or ")
        print("python3 read_netCDF.py file_name")
        sys.exit()
    file_name = args[1]
    data = read(file_name)

