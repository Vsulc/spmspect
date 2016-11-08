# spmspect
SPM analysis of SPECT data based on STATISCOM and ISAS

# requirements
Matlab 2013a+, student version also working
spm 12 and some spm 8 files
recommended
mricron 
mricro
64 bit system with 8 GB RAM - macOS, Windows 7-10, Linux

# Setting up
This script uses Ictal, Interictal, MRI scan to coregister files to the patient's MRI space and then statistical analysis is done in Matlab and SPM 12
initial setup require copying apriori and templates folder from the SPM 8 installation 
and extracting ch2.nii.gz to ch2.nii from the mricron directory to the spm 12 template directory
you have to change paths at the beginnig of all files and and provide path to isas_denominator.nii which can be created from files provided by on the Yale SPECT page, the script will be uploaded later.

# Workflow
Convert data from dicom to nifit (e.g. in dcm2nii, select spm 8 format).
Start *preprocessing files*, select ictal, interictal and MRI data. (MRI is rigidly coregisted to T1 template to fix gross head movement)
SPM will run segmentation which takes 5-10 minutes, brain mask is the created and SPECT data are coregistered to MRI (r\*prefix) and masked brainmask created from segmented files (m\*prefix) with it. Data are then globally normalized to mean 50  (c\*prefix).
*If you don't have a MRI available* select the nomri preprocessing, ictal and interictal data will be coregistered together and  mean will be coregistered to spm spect template, default brainmask will be then applied.
*Run the analysis*
Input spect images that are already coregistered masked and count normalized (cmr\* files) and rMRI file (if you don't have MRI available, don't select anything, data will be coregistered to mean ictal and interictal image.
Analysis run for 1-5 minutes and outputs the file \_Pts\_ which is the result of analysis coregistered back to patients MRI file.
You can view it in your favourite nifti viewer.



# Sources\:
1. Sulc V, Stykel S, Hanson DP, et al. Statistical SPECT processing in MRI-negative epilepsy surgery. Neurology. 2014;82(11):932-939
2. Scheinost D, Blumenfeld H, Papademetris X. An Improved Unbiased Method for Diffspect Quantification in Epilepsy. Proc IEEE Int Symp Biomed Imaging. 2009;2009:927-930.
3. McNally KA, Paige AL, Varghese G, et al. Localizing value of ictal-interictal SPECT analyzed by SPM (ISAS). Epilepsia. 2005;46(9):1450-1464. 
4. Kazemi NJ, Worrell GA, Stead SM, et al. Ictal SPECT statistical parametric mapping in temporal lobe epilepsy surgery. Neurology. 2010;74(1):70-76.


Copyright (c) 2016, Vlastimil Sulc, Department of Neurology, 2nd Faculty of Medicine, Charles University and Motol University Hospital, Prague, Czech Republic
vlsulc+github(at)gmail.com
