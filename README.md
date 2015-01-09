Prostate Sensitivity Analysis
=============================

Functions, modules and classes to perform the prostate lesion imaging sensitivity analysis.  This was originally based on Zach's code, but it has been completely re-written to do considerabley more sub-analyses.  The spreadsheet and original csv files that Zach utilized were replaced by individual files in each patient directory.

```write_histology.py``` was used to write Patient*/Histology/HistologyLesions.txt files.  These files have been augmented from their spreadsheet-exported form to include more than just PCA lesion (now atrophy and bph are included).

```write_ARFI_index_IOS.py``` was used to write Patient*/ARFI_Index_Lesion_IOS.txt files.  These files will now differ from the spreadsheets due to location corrections, etc.

```write_MRI_index.py``` was used to write Patient*/MRI_Images/MRI_Index_Region.txt files.

```arfi_histology_analysis.py``` is the main analysis driver script that generated the output in ```arfi_histology_analysis.md```.  It calls the ```LesionAnalysis``` class, which can be used for interactive patient analysis:

```
from lesion_analysis import LesionAnalysis

P79 = LesionAnalysis(79)

print P79
```
