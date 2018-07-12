The historical climate dataset for the whole of Australia here is created in September 2009 for CWYET project.

Data format: 
• There are 321,457 .csv files in the folder. 
• The files are named by an unique identifier which consists of the coordinates of its centroid and has 9 digits: XXXXXYYYY where XXXXX is longitude*100, YYYY is latitude*(-100).
• All the files contain the historical daily climate data from 1889-01-01 to 2009-08-31 for each 0.05 x 0.05 degree cell across Australia. Each of these 321,457 files contains: date, rainfall(mm), areal potential evapotranspiration (APET)(mm), minimum temperature(degree C), maximum temperature(degree C), incoming solar radiation(MJ/day) and relative humidity.

Sample lines from one file:
Date,Rain_mm,APET_mm,Tmin_C,Tmax_C,Rad_MJ,RH
1889-01-01,0,6.7075,20,28.5,28,0.6185
1889-01-02,0,6.9779,20.5,29,29,0.632
1889-01-03,0,6.9164,19.5,29.5,29,0.6061
1889-01-04,0,7.0069,21,29,29,0.6236
1889-01-05,0,6.8207,21,29,28,0.6236
1889-01-06,0.2,6.9728,21,28.5,29,0.6339
...

Data Sources: 
• The rainfall, minimum temperature, maximum temperature and incoming solar radiation are extracted from SILO (see SILO versions for details).
• The APET is estimated from the above climate data (Tmax, Tmin, Rad and VP) using Morton's ET algorithm.
• The relative humidity (actual vapour pressure divided by saturation vapour pressure) is a by product of the above calculations.

Used Silo Version:
1989-01-01 - 2003-12-31:	s0504 			(SEACI extended 2009-05-09:	s0504)
2004-01-01 - 2004-12-31: 	s0602				(SEACI extended 2009-05-09: s0504)
2005-01-01 - 2005-12-31:	s0708				(SEACI extended 2009-05-09: s0607)
2006-01-01 - 2006-12-31:	s0708				(SEACI extended 2009-05-09: s0708)
2007-01-01 - 2008-12-31:	s0905 			(SEACI extended 2009-05-09: s0905)
(folders 2007 and 2008 are copied from \\file-wron\TimeSeries\Climate\silo.eoc.csiro.au\www\data\silo2\rain\flt\ on 19 May 2009.)
2009-01-01 - 2009-08-31:	s0909			(SEACI extended to 2009-05-09: s0905 for from 2009-01-01 to 2009-05-09)

Note:
• SEACI extended to 2009-05-09 is another dataset which only include the cells within MDB, NSW and Victoria. The data period for this data set is from 1889-01-01 to 2009-05-09.
• 2007 and 2008 rainfall data are extracted from \\file-wron\TimeSeries\Climate\silo.eoc.csiro.au\www\data\silo2\rain\flt\ on 19 May 2009. The rest of the data are extracted from \\file-wron\TimeSeries\Climate\usilo.
• 2009 data are latest from NR&M, \\Data2-wron\cy_working\Jai\2009, we call it version s0909, details as following:
				NCOLS   841                                       
				NROWS   681                                       
				XLLCENTER  112.00                                 
				YLLCENTER  -44.00                                 
				CELLSIZE  0.05                                    
				NODATA_VALUE  -99.9                               
				Geodetic:WGS84     0                              
				Elem:rai 20090101                                 
				Type:Daily  1                                     
				Units:mm  0.1                                     
				Created 200909020546                              
				BExtract 20050526                                 
				Algorithm:krige.1  0                              
				Interpolation 2324904                             
				Copyright:StateOfQueensland(NR&M)  1992 

Contact:
This dataset is created as part of CWYET project. Please contact project leader Jai Vaze (jai.vaze@csiro.au) or Jin Teng (jin.teng@csiro.au) or Bill Wang (bill.wang@csiro.au) should you have any enquires about the data.

Reference:
Please remember to properly acknowledge SILO and CWYET project: 

Jeffrey, S. J., Carter, J. O., Moodie, K. B. and Beswick, A. R., 2001. Using spatial interpolation to construct a comprehensive archive of Australian climate data. Environmental Modelling & Software 16 (4), 309-330.

Vaze, J., Chiew, F. H. S., Perraud, JM., Viney, N., Post, D. A., Teng, J., Wang, B., Lerat, J. and Goswami, M., 2011. Rainfall-runoff modelling across southeast Australia: datasets, models and results. Australian Journal of Water Resources, 14 (2), 101-116.

Vaze, J., Perraud, J-M., Teng, J., Chiew, F. H. S., Wang, B. and Yang, Z., 2011. Catchment Water Yield Estimation Tools (CWYET). 34th IAHR World Congress, 27th June to 1st July, Brisbane, Australia.

