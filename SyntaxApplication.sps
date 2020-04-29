* Encoding: UTF-8.
* Encoding: .

*Open the datasets such that the datasets are in the following order:
*DataSet1 = main dataset
*DataSet2 = nclass dataset
*DataSet3 = sex dataset

*Make a unique ID for every observation (teacher id, year and classnumber) in main dataset

DATASET ACTIVATE DataSet1.
COMPUTE New_ID=ID_T * 1000000 + jaar * 10000 + maand*100 + klasnr.
EXECUTE.

*Make a unique ID for every observation (teacher id, year and classnumber) in dataset with classnr

DATASET ACTIVATE DataSet2.
COMPUTE New_ID=ID_T * 1000000 + jaar * 10000 + maand*100 + klasnr.
EXECUTE.

*Because two cases have a missing teacher number and thus ID_T and New_ID cannot be calculated,
we cannot merge the dataset yet. Therefore we will filter out the cases without ID in both datasets.

*Inspect the frequency of missing cases.
DATASET ACTIVATE DataSet1.
FREQUENCIES VARIABLES=New_ID
  /FORMAT=DFREQ
  /ORDER=ANALYSIS.

*Inspect the missing cases and filter them out.

DATASET ACTIVATE DataSet1.
USE ALL.
COMPUTE filter_$=(New_ID  > 0).
VARIABLE LABELS filter_$ 'New_ID  > 0 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

DATASET ACTIVATE DataSet2.
USE ALL.
COMPUTE filter_$=(New_ID  > 0).
VARIABLE LABELS filter_$ 'New_ID  > 0 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

*Delete the cases without teacher ID from both datasets.

DATASET ACTIVATE DataSet1.
FILTER OFF.
USE ALL.
SELECT IF (New_ID  > 0).
EXECUTE.

DATASET ACTIVATE DataSet2.
FILTER OFF.
USE ALL.
SELECT IF (New_ID  > 0).
EXECUTE.

*Merge the two datasets (i.e. main dataset and nclass dataset)

DATASET ACTIVATE DataSet1.
STAR JOIN
  /SELECT t0.schnr, t0.docnr, t0.jaar, t0.maand, t0.klasnr, t0.ID_T, t0.schtype, t0.leerjr, t0.erv, 
    t0.AG, t0.COM, t1.nclass
  /FROM * AS t0
  /JOIN 'DataSet2' AS t1
    ON t0.New_ID=t1.New_ID
  /OUTFILE FILE=*.

*Change the format so errors are visible

formats erv (f6.2).

*Delete irregular years of teachers, so teachers that do not have one of 1, 2 or 3 for experience (erv)
*Do a check first.

DATASET ACTIVATE DataSet1.
USE ALL.
COMPUTE filter_$=(erv = 1  | erv = 2  |  erv = 3).
VARIABLE LABELS filter_$ 'erv = 1  | erv = 2  |  erv = 3 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

*Delete them from the file.

FILTER OFF.
USE ALL.
SELECT IF (erv = 1  | erv = 2  |  erv = 3).
EXECUTE.

*Give a 1 for every first teacher

sort cases by ID_T.
execute.
Match files 
 /file=*
 /by ID_T
/first=N_T.
freq N_T.

*Give a 1 for every first measurement per year per teacher

sort cases by ID_T erv.
execute.
Match files 
 /file=*
 /by ID_T erv
/first=N_Y.
freq N_Y.

*Specify the number of measurements per teacher

AGGREGATE
  /OUTFILE=* MODE=ADDVARIABLES
  /BREAK=ID_T
  /N_YF=SUM(N_Y).

*Put a filter to select every teacher once and then see the frequencies for the years.

USE ALL.
COMPUTE filter_$=(N_T = 1).
VARIABLE LABELS filter_$ 'N_T = 1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

FREQUENCIES N_YF.

USE ALL.

*The following changes are made based on the manual checking of the original data in the folders,
*since that was suggested by Mieke Brekelmans

*Delete duplicate cases for teacher 1004

FILTER OFF.
USE ALL.
SELECT IF (New_ID ~= 1004820502.00).
EXECUTE.

FILTER OFF.
USE ALL.
SELECT IF (New_ID ~= 1004830503.00).
EXECUTE.

FILTER OFF.
USE ALL.
SELECT IF (New_ID ~= 1004830504.00).
EXECUTE.

FILTER OFF.
USE ALL.
SELECT IF (New_ID ~= 1004840503.00).
EXECUTE.

FILTER OFF.
USE ALL.
SELECT IF (New_ID ~= 1004840504.00).
EXECUTE.

*Delete duplicate cases for year 1 of teacher 5035

FILTER OFF.
USE ALL.
SELECT IF (New_ID ~= 5035850606.00).
EXECUTE.

FILTER OFF.
USE ALL.
SELECT IF (New_ID ~= 5035850607.00).
EXECUTE.

FILTER OFF.
USE ALL.
SELECT IF (New_ID ~= 5035850608.00).
EXECUTE.

FILTER OFF.
USE ALL.
SELECT IF (New_ID ~= 5035850609.00).
EXECUTE.

FILTER OFF.
USE ALL.
SELECT IF (New_ID ~= 5035850610.00).
EXECUTE.

*Calculate a mean agency and mean communion score over all classes of the teacher.

AGGREGATE
  /OUTFILE=* MODE=ADDVARIABLES
  /BREAK=erv ID_T
  /AG_mean=MEAN(AG) 
  /COM_mean=MEAN(COM).

*Calculate the number of students per year per teacher

AGGREGATE
  /OUTFILE=* MODE=ADDVARIABLES
  /BREAK=ID_T erv
  /nclass_sum=SUM(nclass).

*Only select the teachers that have data over the three years

DATASET ACTIVATE DataSet1.
USE ALL.
COMPUTE filter_$=(N_YF = 3).
VARIABLE LABELS filter_$ 'N_YF = 3 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

*Copy the dataset with the teachers that have data over the three years to calculate the number of
*students per teacher.

DATASET ACTIVATE DataSet1.

DATASET COPY  AgencyCommunion.
DATASET ACTIVATE  AgencyCommunion.
FILTER OFF.
USE ALL.
SELECT IF (N_YF = 3).
EXECUTE.
DATASET ACTIVATE  DataSet1.

*Check the frequency for the number of teachers and for the number of students per teacher.
DATASET ACTIVATE  AgencyCommunion.
FREQUENCIES VARIABLES=N_T
  /ORDER=ANALYSIS.

*Find the number of students that scored a teacher per year of experience.

DATASET ACTIVATE AgencyCommunion.
USE ALL.
COMPUTE filter_$=(erv = 1).
VARIABLE LABELS filter_$ 'erv = 1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

DESCRIPTIVES VARIABLES=nclass_sum
  /STATISTICS=MEAN STDDEV MIN MAX.

DATASET ACTIVATE AgencyCommunion.
USE ALL.
COMPUTE filter_$=(erv = 2).
VARIABLE LABELS filter_$ 'erv = 2 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

DESCRIPTIVES VARIABLES=nclass_sum
  /STATISTICS=MEAN STDDEV MIN MAX.

DATASET ACTIVATE AgencyCommunion.
USE ALL.
COMPUTE filter_$=(erv = 3).
VARIABLE LABELS filter_$ 'erv = 3 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

DESCRIPTIVES VARIABLES=nclass_sum
  /STATISTICS=MEAN STDDEV MIN MAX.

*Go back to the original datase to continue the formatting

DATASET ACTIVATE Dataset1.

*Select only the first line of measurement for every teacher since we now have the aggregated measure.

FILTER OFF.
USE ALL.
SELECT IF  (NOT(N_Y=0)).
EXECUTE.

*Delete unnecessary variables for main analysis.

DELETE VARIABLES schnr docnr jaar maand klasnr schtype leerjr ag com N_T N_Y filter_$ nclass nclass_sum New_ID.

*Make a wide format

CASESTOVARS
  /ID=ID_T
  /GROUPBY=VARIABLE.

*Make sure the right year/data is used in the rigjht variable when the teacher does not have all 3 years
*Agency

IF(erv.2 = 3) AG_mean.3 = AG_mean.2.
EXECUTE.

IF(erv.1 = 2) AG_mean.2 = AG_mean.1.
EXECUTE.

IF(erv.1 = 3) AG_mean.3 = AG_mean.1.
EXECUTE.

*Delete the wrong values at timepoint 1 (because it belongs to 2 and 3).

IF(erv.1 = 2) AG_mean.1 = 99999.
EXECUTE.

IF(erv.1 = 3) AG_mean.1 = 99999.
EXECUTE.

DO IF (erv.1 ~= 1).
RECODE AG_mean.1 (99999=SYSMIS).
END IF.
EXECUTE.


*Make sure the right year/data is used when the teacher doesn't have all 3 years
*Communion

IF(erv.2 = 3) COM_mean.3 = COM_mean.2.
EXECUTE.

IF(erv.1 = 2) COM_mean.2 = COM_mean.1.
EXECUTE.

IF(erv.1 = 3) COM_mean.3 = COM_mean.1.
EXECUTE.

*Delete the wrong values at timepoint 1 (because it belongs to 2 and 3).

IF(erv.1 = 2) COM_mean.1 = 99999.
EXECUTE.

IF(erv.1 = 3) COM_mean.1 = 99999.
EXECUTE.

DO IF (erv.1 ~= 1).
RECODE COM_mean.1 (99999=SYSMIS).
END IF.
EXECUTE.

*Now change the erv variables so that they are correct

IF(erv.2 = 3) erv.3 = 3.
EXECUTE.

IF(erv.1 = 2) erv.2 = 2.
EXECUTE.

IF(erv.1 = 3) erv.3 = 3.
EXECUTE.

IF(erv.1 = 2) erv.1 = 99999.
EXECUTE.

IF(erv.1 = 3) erv.1 = 99999.
EXECUTE.

IF(erv.2 = 3) erv.2 = 99999.
EXECUTE.

*Now clean the agency and communion based on the right ervaring
*Cases with experience that is not 1, 2 or 3 are cleaned (e.g. erv = 1.50).

IF(erv.1 ~= 1) AG_mean.1 = 99999.
EXECUTE.
IF(erv.1 ~= 1) COM_mean.1 = 99999.
EXECUTE.
IF(erv.1 ~= 1) erv.1 = 99999.
EXECUTE.

IF(erv.2 ~= 2) AG_mean.2 = 99999.
EXECUTE.
IF(erv.2 ~= 2) COM_mean.2 = 99999.
EXECUTE.
IF(erv.2 ~= 2) erv.2 = 99999.
EXECUTE.

IF(erv.3 ~= 3) AG_mean.3 = 99999.
EXECUTE.
IF(erv.3 ~= 3) COM_mean.3 = 99999.
EXECUTE.
IF(erv.3 ~= 3) erv.3 = 99999.
EXECUTE.

*Change 99999 into sysmis, because that is easier when transferring to R.

RECODE erv.1 erv.2 erv.3 (99999=SYSMIS).
EXECUTE.

RECODE AG_mean.1 AG_mean.2 AG_mean.3 (99999=SYSMIS).
EXECUTE.

RECODE COM_mean.1 COM_mean.2 COM_mean.3 (99999=SYSMIS).
EXECUTE.

*Find the median of communion at the first year to make the groups.

FREQUENCIES VARIABLES=COM_mean.1
  /FORMAT=NOTABLE
  /STATISTICS=STDDEV MINIMUM MAXIMUM MEAN MEDIAN SKEWNESS SESKEW KURTOSIS SEKURT
  /ORDER=ANALYSIS.

*New group based on the value of the first year for the people that have data over two years

RECODE COM_mean.1 (Lowest thru 0.218770=0) (0.218770 thru Highest=1) INTO COMGroup.
EXECUTE.

FREQUENCIES VARIABLES=COMGroup
  /STATISTICS=MEAN STDDEV MIN MAX.

*Add sex to datafile

DATASET ACTIVATE DataSet1.
STAR JOIN
  /SELECT t0.N_YF, t0.erv.1, t0.erv.2, t0.erv.3,
    t0.AG_mean.1, t0.AG_mean.2, t0.AG_mean.3, t0.COM_mean.1, t0.COM_mean.2, 
    t0.COM_mean.3, t0.COMGroup, t1.sex
  /FROM * AS t0
  /JOIN 'DataSet3' AS t1
    ON t0.ID_T=t1.ID_T
  /OUTFILE FILE=*.

*Only select the teachers that have measurements at all three timepoints 

USE ALL.
COMPUTE filter_$=(erv.1 = 1 & erv.2 =2  & erv.3 = 3).
VARIABLE LABELS filter_$ 'erv.1 = 1 & erv.2 =2  & erv.3 = 3 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

FILTER OFF.
USE ALL.
SELECT IF  (erv.1 = 1 & erv.2 =2  & erv.3 = 3).
EXECUTE.

*Calculate the number of males and the number of females.

FREQUENCIES VARIABLES=sex
  /ORDER=ANALYSIS.

*Calculate the number of teachers in the lower communion group and the higher communion group.

FREQUENCIES VARIABLES=COMGroup
  /ORDER=ANALYSIS.


