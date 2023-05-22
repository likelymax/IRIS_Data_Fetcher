# IRIS Seismic Data Fetcher

This is a code repository to download seismic waves from [IRIS](https://www.iris.edu/hq/). It contains several modules and scripts which perform different tasks from fetching seismic data to converting it into a specific format. Here's a brief overview of the files:

## Code Overview

* **FetchData** and **FetchEvent**: These modules are downloaded directly from IRIS to obtain seismic data and event details, respectively.  
* **find_event.c**: This C program locates the period during which a potential seismic event might occur. The output from this module is the boundaries of this event period.  
* **find_event_com**: This bash script uses 'find_event' to identify the potential event period and FetchEvent to extract the precise event details.  
* **read_sac.py**: A Python script that retrieves the station information and the duration of the seismic wave you wish to download, based on the event details.  
* **irisFetch_YW**: This bash script downloads data from IRIS. The seismic wave data is downloaded in miniSEED time series format.  
* **mseed2sac**: A utility program obtained from [__iris-edu__](https://github.com/iris-edu/mseed2sac/releases) that converts miniSEED time series data into Seismic Analysis Code (SAC) format.  

Please note, the FetchData module has been modified to limit the size of each seismic wave to a maximum of 5e+6 bits.

## Getting Started

To start, you will first need to download the seismic data from IRIS. The required station information should be saved first. I name it as 'name.txt'. The 'read_sac.py' script will use this as input to generate a file containing the station information, including Net, station, stlo, and stla.

Next, use the output information (which includes the time periods of the events you're interested in) to fetch the seismic wave data from IRIS with FetchData. Following this, use 'find_event_com' to extract the event details and integrate this information into the sac file.

Before using 'find_event', you need to compile it into an executable file. This can be done using the command _make find_event_. The input file for 'find_event' should be a file detailing the time periods of the events you want to fetch. The example of input file for _find_event_ is 'event.dat'

To use 'find_event', use the following syntax:

./find_event -Eevent.dat -T$time -D$period/$delta -Otime.txt

Here:

* __-E__ designates the event information file. 'event.dat' is an example of this file that you want the code to process.
* __-T__ specifies the event date/time (year/month/day/hour/minute/sec).
* __-D__ denotes the duration/samplerate.
* __-O__ indicates the output file for the time period. The output format is compatible with FetchEvent.

