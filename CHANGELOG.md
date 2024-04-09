## 1.5.1 - 2024-04-08
 - Changed default filterSupplemental to tru for sorting task
 - Changed how the above task handled flag (it appears that previous version 
   did not have the right implementation for filterSupplemental flag check)
## 1.5.0 - 2024-03-08
 - updated imported bwaMem (to bwamem2 2.2.1) and star (to 2.3.0)
 - updated samtools from 1.9 to 1.14
 - flag for removing supplemental alignments (ON by default to prevent failures when using updated STAR)
## 1.4.0 - 2022-06-15
 - assembly-specific modules are specified inside the workflow
## 1.3.0 - 2022-03-08
 - multi-lane STAR support for WT/MR data, some re-design of the workflow
## 1.2.1 - 2022-02-08
 - Added support for TS, EX in the workflow
## 1.2 - 2022-01-18
 - Adding ability to run on WT data (using imported STAR)
## 1.1.1 - 2021-05-31
 - Migrating to vidarr
## 1.1 - 2020-10-21
 - Added an ability to remove tags on a list, bwaMem is now a sub-workflow, .json report with totals for different tags is also provisioned 
## 1.0 - 2020-05-07
 - Initial Release
