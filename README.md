Zach's "final" processing code for the ARFI lesion : histopath analysis.

regions27_test.m is the driver script.

Google Spreadsheet Master: https://docs.google.com/spreadsheets/d/1YqP2tTq7cqs_wg_ffRaeL_xsTADUAczxyTLCfjle3NI/edit#gid=1159472451

Notes from Zach:

1. loads 49X55 Histopathology worksheet converts to standardized 49X28X4 array
   where:
 * 3dim 1: binary occupation by tumor [built from dim 2]
 * 3dim 2: tumor volume
 * 3dim 3: gleason grade
 * 3dim 4: posterior/anterior location [built from dim 2]

2. repeats similar action for ARFI worksheet, instead of 4 entries in third dim
   like Histopathology workseet, only 2 entries: binary tumor occupation +
   indices of suspicion

3. same comparison technique used to determine if hit or miss as old code

4.  filter technique: checks if a 1 is contained in array after one lesion has
    been compared against all ARFI ROIs, if 1 is found then ARFI hit the
    lesion, otherwise ARFI missed the lesion [this 1 or 0 is mapped to a second
    49 lesion X 1 binary matrix]. The array is then cleared for next lesion. 

5. The 49X1 binary matrix is built up into a 5 column datasheet using the 49X28X4 array
 * column 1: binary tumor hit data
 * column 2: patient index
 * column 3: volume
 * column 4: gleason
 * column 5: posterior/anterior

6. this spreadsheet is then filtered into various sub-spreadsheets which are
   used to calculate sensitivities. These sensitivities are printed, and a bar
   graph is constructed.

7. future effort: repeat for PPV and IOS calculations.
