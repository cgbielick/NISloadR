# NISloadR

## Summary

An R package to load, clean, and process the National Inpatient Sample (NIS). 

Works with ASCII files for at least years 1993-2022. May work with future years prior to any package updates if formatting stays the same. Saves ASCII as fst, bypassing proprietary Stata, SAS, and SPSS programs. Will still require the ".Do" Stata script files provided by HCUP as an objective reference. 

Future iterations: to include parallelization for lower memory systems, cleaning tools, plug-and-play ICD identification, and machine learning for forecasting and clustering. 
- Name a vector of years of data to include which doesn't need to be repeated

## Functional: 
For all four functions, specify three parameters: 1) "do_path" (include file name for the .Do Stata file), 2) "asc_path" (include NIS ASC file name), and 3) "output_direct" (specify output directory).
- load_nis_core()
- load_nis_severity()
- load nis_hospital()
- load_nis_grps()


## In progress: 
flag_nis() 
- Must provide the numeric year of data (one at a time, will accommodate a vector of years in the future). 
- If 2011 or earlier, default will be to require a path to the new DISCWT's. Can be overridden by argument old_weights = 1 (default 0) to specify the intentional use of the old weights. 
  - For 1993-2011: https://hcup-us.ahrq.gov/db/nation/nis/trendwghts.jsp
  - Must point to a directory containing the FileSpecifications_NIS_###_etc.TXT and the corresponding ASC file.
  - Eventually this part will be automated with permission of AHRQ. 
- If 2015, the diagnosis codes are in the grps file, not the core. Requires special attention due to the split from ICD-9-CM in q1-q3 to ICD-10-CM in q4.

Columns to keep:
- DX (or I10_DX) columns: limit to a finite number of ICD columns, default is all. 0 for none. 
- PR columns: limit to a finite number of PR columns, default is all. 0 for none.
- PRDAY columns: limit to a finite number of PRDAY columns, default is all. 0 for none (should be 0 if PR columns is 0).
- All other columns: give a single character vector of column names to filter. Can include a single large vector of columns which MIGHT be present if processing multiple years of data. 

- Adults = TRUE or (default) FALSE. Keeps both >=18 and NA (where AGE_NEONATE is also NA). These load in by default. The filter will occur AFTER loading in the data with a specified vector of columns, only dropping AGE or AGE_NEONATE later if either is not in that vector.
- Cost: default NULL, if provided a path to the HCUP charges to cost ratio file it will convert charges to estimated cost. Only available from 2001 to present


Codes for which to search: will search over any remaining 
- Any column: default is 1. Will create a column flagging a 0 or 1 if any match occurred. 
- Sum column: default is 1. Will create a column which sums how many matches there were per discharge. Also applies to a list of vectors.
- ICD-9-CM: provide a vector (or list of vectors) of strings to match in any of the DX columns. 
- ICD-10-CM: provide a vector (or list of vectors) of strings to match in any of the I10_DX columns.
- ICD-9-PCS: provide a vector (or list of vectors) of strings to match in any of the PR columns.
- ICD-10-PCS: provide a vector (or list of vectors) of strings to match in any of the I10_PR columns. 



Future: 
flag_nis_ccs()
- utilizes the dx_pr_grps files, DXCCSn for 2015q3 and earlier; I10 for 2018 and later (2016-2017 the data was not present)
- can search for CCS PR as well

Missing data processing: imputation and filtering based on percentages 

Survey object creation 
- May account for lonely PSU

Summary statistics: 
- Utilize https://hcup-us.ahrq.gov/db/nation/nis/nissummstats.jsp as a benchmark 
