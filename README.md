Contained in this repository are two tools for reading and plotting COAPS Simplified Daily Scatterometer Swath files, *map_swath_data.py* and *read_netCDF.py*.

---
  
# map_swath_data.py - COAPS Simplified Daily Scatterometer Swath Plotting Tool
This is a tool to plot wind surface speed and direction over the oceans using the COAPS Simplified Daily Scatterometer Swath files. 

## Prerequisites and Versions
This tool uses the following software and versions:
```
python -  Version 3.6.6
cartopy - Version 0.16.0
numpy - Version 1.14.5
matplotlib - Version 2.2.3
```

## Set Up
You must first download the two files "map_swath_data.py" and "read_netCDF.py", or clone this git repository.

To clone this repository on a remote server, enter the following command in your terminal:
```
git clone https://github.com/marine-data-center/coaps_daily_swath.git
```

The script requires a dataset file (details in the usage section to follow), which can be found and downloaded at the following site:
```
https://mdc.coaps.fsu.edu/swath
```

Files downloaded from this site are compressed with extension .tgz and can be decompressed via the following command:
```
tar -xvf filename.tgz
```

## Usage
The most basic command required to generate a plot is:
```
$  python3 map_swath_data.py path/to/filename.nc --[vectors|barbs|flags|speed|param param_name]
```

**Required Arguments**

- Filename or path to filename, ending in .nc
- One and only one of the following: --vectors, --barbs, --flags, --speed, --param *param_name*
    - where *param_name* is the (STRING) name of the parameter to plot as it appears in the netCDF file.
        - Note: you cannot plot the one-dimensional time variable. 
    - --vectors creates a plot with vectors
    - --barbs creates a plot with wind barbs
    - --flags creates a plot showing flagged (or bad) data points, based on the COAPS quality control flag
    - --speed creates a scattterplot of wind speed

**Optional Arguments**

- --lon_0 LON_0
    - type: INT
    - Center longitude of map, in range of -180 to 180 (degrees East)
    - 0 by default
- --vmin VMIN
    - type: INT
    - The minimum colorbar value. Data values outside of this will be given special out-of-range value. 
    - Note: this value is assumed to be in m/s for plots involving vectors and speed, and in knots for plots with barbs
    - Default: 0 for plots with barbs, vectors, or speed; Otherwise, it is set to the minimum of the data
    - vmin cannot be used when plotting quality flags
- --vmax VMAX
    - type: INT
    - The maximum colorbar value. Data values outside of this will be given special out-of-range value.
    - Note: this value is assumed to be in m/s for plots involving vectors and speed, and in knots for plots with barbs
    - Default: the maximum value of the data
    - vmax cannot be used when plotting quality flags
- --n N
    - type: INT
    - Plots every nth point in the file. Modify this value to see more or less dense data points on your plot
    - This value scales based on the bounding box as well as the how large the cross_track dimension is of the file
        - See the 'print_details' flag to display such value
- Bounding box
    - In the form --north INT, --east INT, --west INT, --south INT
    - All four values are required at once, and cannot be used in conjunction with a lon_0 argument.
    - On plots crossing the dateline, enter longitude values still according to [-180,180]
- --save_as "SAVE_AS"
    - type: STRING
    - Specify the image name or full path + image name to be saved.
    - Accepts either .png or .jpg extensions
    - Default name is [inputfile]_[vectors|barbs|speed|flags].jpg
    - Default location is current directory.
- --title "TITLE"
    - type: STRING
    - Title of the plot.
- --barb_size BARB_SIZE
    - type: FLOAT
    - The size of barb icons on the plot.
    - Default icon size on global plots is 2.25.
    - For zoomed images, the icons are scaled based on bounding box area.
        - To increase/decrease the size, we recommend using the 'print_details' flag in order to see the current size, and incrementing/decrementing by 0.25-0.75
- --vector_size VECTOR_SIZE
    - type: INT
    - The scaling value of vector icons.
    - Default value on global plots is 1000, and 600 for zoomed
    - For reference, try incrementing/decrementing by 200. 
    - Note: the smaller the value given, the larger the vector icon.
- --print_details
    - a flag for printing details of the plot being created, printed on the command line
    - Since some values (such as n, and the size of icons) scale based on the bounding box, they are not always fixed.
        - This flag is useful to see what those values are for the plot you are creating if you want to increase/decrease these values

**Examples**

To get a global map of wind barbs using the file ASCATA_25_2018256.nc, titled "ASCATA 25km Wind Barbs Sept 13, 2018", with larger sized icons:
```
python3 map_swath_data.py ASCATA_25_2018256.nc --barbs --title "ASCATA 25km Wind Barbs Sept 13, 2018" --barb_size 3.0
```

To get a global map of wind vectors, with a central longitude of 160 W, and saving the image as "ASCATA_vector.jpg":
```
python3 map_swath_data.py ASCATA_25_2018256.nc --vectors --lon_0 -160 --save_as "ASCATA_vector.jpg"
```

To plot wind speed over central america, plotting every 15th point:
```
python3 map_swath_data.py ASCATA_25_2018256.nc --speed --n 15 --north 30 --south 10 --west -100 --east -60
```

To create a global map of wind barbs, coloring only wind speeds between 10 and 25 knots:
```
python3 map_swath_data.py ASCATA_25_2018256.nc --barbs --vmax 25 --vmin 10
```

To plot the quality control flags over the entire globe, and printing the details of the plot to the screen:
```
python3 map_swath_data.py QSCAT_12_2009002.nc --flags --print_details
```

**Getting Help**

Execute the following command to get explanation on arguments:
```
$  python3 map_swath_data.py --help
```

## Authorship
This code was written by Jocelyn Elya, Daniel Kane, and Bethany Sanders
 
---
# read_netCDF.py - COAPS Simplified Daily Scatterometer Swath netCDF File Reader
This is python read code which reads variables into memory from the COAPS Simplified Daily Scatterometer Swath netCDF files.

This code is called upon by the above plotting tool to read the data from a given netCDF file. 

## Prerequisites and Versions
```
python - version 2.7.15
netCDF4 - version 1.2.9
```

## Set Up
Download the file "read_netCDF.py" or clone this repository.
To clone this repository on a remote server, enter the following command in your terminal:
```
git clone https://github.com/marine-data-center/coaps_daily_swath.git
```

The script requires a COAPS Simplified Daily Scatterometer Swath file, which can be found and downloaded at the following site:
```
https://mdc.coaps.fsu.edu/swath
```

Files downloaded from this site are compressed with extension .tgz and can be decompressed via the following command:
```
tar -xvf filename.tgz
```

## Usage
Command line:
```
python read_netCDF.py file_name
```
or
```
python3 read_netCDF.py file_name
```    
Imported:
```
from read_netCDF import read
data = read(file_name)
```

**Argument**

- file_name = string value containing the relative or absolute path of the netCDF file

**Returns**

- variable_dict = dictionary where the keys are variable names and the value are masked numpy arrays of the variable's data as read from the file

## Variables in the netCDF files
**time**
A one dimensional LONG array containing the mean time of the 
observations. Given as the seconds since 1990-01-01 00:00:00.
Each value in the array corresponds to one row of observations:
all times in across track (perpendicular to the satellite track)
are equal. The index corresponds to the along track time and position.

**lat**
A two dimensional FLOAT array containing the latitude of the 
wind vectors. The first index corresponds to the along track
position, and the second index corresponds to the cross
track position.
         
**lon**
A two dimensional FLOAT array containing the longitude of the 
wind vectors. The first index corresponds to the along track
position, and the second index corresponds to the cross
track position.
         
**eastward_wind**
Zonal equivalent neutral wind vector component. Positive values
indicate an eastward movement. The first index corresponds to
the along track position, and the second index corresponds to
the cross track position.

**northward_wind**
Meridional equivalent neutral wind vector component. Positive
values indicate northward movement. The first  index corresponds
to the along track position and the second index corresponds to 
the cross track position
          
**wind_speed**
A two dimensional FLOAT array containing the equivalent neutral 
wind speed. The first index corresponds to the track position, 
and the second index corresponds to the cross track position.
          
**wind_to_direction**
A two dimensional FLOAT array containing the direction of the
wind vector (oceonographic convention). The first index corresponds to the along track 
position, and the second index corresponds to the cross track position.
          
**simplified_wvc_quality_flag**
A two dimensional INTEGER array with two values, 0 or 1.
This flag is a summary of the product developer's non-rain
flag, with 0 indicating good data, and 1 indicating a missing
flag value, highly suspect data, or data of poor quality. See
the metadata for this variable to see flags that are combined
to produce this summary flag.
          
**rain_impact (QSCAT and RapidSCAT only)**
A two dimensional FLOAT array that attempts to indicate the impact of
rain on a zero to 100 scale. At this time, the interpretation of these
numbers is poorly described for Ku band, and not provided for C band. 
We anticipate this changing in the future; however, at this time JPL 
recommends most data users to avoid using the rain_impact variable. The 
only exception to this rule is for users who are trying to do very 
specific studies that are highly adversely impacted by rain.

**distance_from_coast (QSCAT and RapidSCAT only)**
A two dimensional FLOAT array containing the distance from
the nearest coastline in kilometers. The first index
corresponds to the along track position and the second index
corresponds to the cross track position.
        
**ice_prob (ASCAT A and ASCAT B only)**
A two dimensional FLOAT array containing the ice probability.
The first index corresponds to the along track position and
the second index corresponds to the cross track position. 

## Authorship 
This code was written by Arturo Valery and Jocelyn Elya.

---

# Support
We would like to acknowledge support from NASA for the Ocean Vector Winds Science Team.
Please direct concerns, questions, and suggestions to the following email: mdc@coaps.fsu.edu
