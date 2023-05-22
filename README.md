# IRIS Seismic Data Fetcher

This is a code repository to download seismic waves from IRIS. It contains several modules and scripts which perform different tasks from fetching seismic data to converting it into a specific format. Here's a brief overview of the files:

## Code Overview

* **FetchData** and **FetchEvent**: These modules are downloaded directly from IRIS to obtain seismic data and event details, respectively.  
* **find_event.c**: This C program locates the period during which a potential seismic event might occur. The output from this module is the boundaries of this event period.  
* **find_event_com**: This bash script uses 'find_event' to identify the potential event period and FetchEvent to extract the precise event details.  
* **read_sac.py**: A Python script that retrieves the station information and the duration of the seismic wave you wish to download, based on the event details.  
* **irisFetch_YW**: This bash script downloads data from IRIS. The seismic wave data is downloaded in miniSEED time series format.  
* **mseed2sac**: A utility program obtained from iris-edu that converts miniSEED time series data into Seismic Analysis Code (SAC) format.  
