function hawk_mri_spect_preprocessing_nnm

% If you find this script useful consider citing the original articles:
%
% 1. Sulc V, Stykel S, Hanson DP, et al. Statistical SPECT processing in MRI-negative epilepsy surgery. Neurology. 2014;82(11):932-939
% 2. Scheinost D, Blumenfeld H, Papademetris X. An Improved Unbiased Method for Diffspect Quantification in Epilepsy. Proc IEEE Int Symp Biomed Imaging. 2009;2009:927-930.
% 3. McNally KA, Paige AL, Varghese G, et al. Localizing value of ictal-interictal SPECT analyzed by SPM (ISAS). Epilepsia. 2005;46(9):1450-1464. 
% 4. Kazemi NJ, Worrell GA, Stead SM, et al. Ictal SPECT statistical parametric mapping in temporal lobe epilepsy surgery. Neurology. 2010;74(1):70-76.
%
% Licensed under The BSD 3-Clause License 
% Copyright (c) 2016, Vlastimil Sulc, Department of Neurology, 2nd Faculty
% of Medicine, Charles University and Motol University Hospital, Prague,
% Czech Republic
% vlsulc@gmail.com

addpath('C:\Users\INM\Documents\IctalSPECT\spm12')
5 addpath('E:\Neuro\spm12')
% clear classes

spm_jobman('initcfg');
DELETE_Temp_Files = 1;
[mri.name, mri.pathname, ~] = uigetfile('*MRI*.nii', 'Pick MRI file');
[ictal.name, ictal.path, ~] = uigetfile('*Ictal*.nii', 'Select Ictal file');
[interictal.name, interictal.path, ~] = uigetfile('*Interictal*.nii', 'Select Interictal file');
ictal.fullfile = fullfile(ictal.path,ictal.name);
interictal.fullfile = fullfile(interictal.path,interictal.name);
cd(mri.pathname)

matlabbatch{1}.spm.spatial.preproc.channel.vols = {fullfile(mri.pathname,mri.name)};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,1')}; 
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,2')};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,3')};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,4')};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,5')};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,6')};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];

spm_jobman('run',matlabbatch);
clear matlabbatch



% create binary brain mask

matlabbatch{1}.spm.util.imcalc.input = {
                                       fullfile(mri.pathname,['c1' mri.name])
                                       fullfile(mri.pathname,['c2' mri.name])
                                       fullfile(mri.pathname,['c3' mri.name])
                                        };
matlabbatch{1}.spm.util.imcalc.output = 'brainmask_0';
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = 'i1+i2+i3';
%matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

matlabbatch{2}.spm.spatial.smooth.data(1) = {fullfile(mri.pathname,'brainmask_0.nii')};
matlabbatch{2}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch{2}.spm.spatial.smooth.dtype = 0;
matlabbatch{2}.spm.spatial.smooth.im = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = 's';

matlabbatch{3}.spm.util.imcalc.input(1) = {fullfile(mri.pathname,'sbrainmask_0.nii')};
matlabbatch{3}.spm.util.imcalc.output = 'brainmask_1';
matlabbatch{3}.spm.util.imcalc.outdir = {''};
matlabbatch{3}.spm.util.imcalc.expression = 'i1>0.5';
matlabbatch{3}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{3}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{3}.spm.util.imcalc.options.mask = 0;
matlabbatch{3}.spm.util.imcalc.options.interp = 1;
matlabbatch{3}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',matlabbatch);
clear matlabbatch

%%  rigidly COREG to MRI template
disp('Coregistering MRI to template (rigid coregistration)...')
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {fullfile(spm('Dir'),'canonical','ch2.nii,1')};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {fullfile(mri.pathname,mri.name)};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {
                                                fullfile(mri.pathname,['c1' mri.name])
                                                fullfile(mri.pathname,['c2' mri.name])
                                                fullfile(mri.pathname,['c3' mri.name])
                                                    };
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
spm_jobman('run',matlabbatch);
clear matlabbatch


%% create binary brain mask
disp('Creating binary mask...')

matlabbatch{1}.spm.util.imcalc.input = {
                                       fullfile(mri.pathname,['rc1' mri.name])
                                       fullfile(mri.pathname,['rc2' mri.name])
                                       fullfile(mri.pathname,['rc3' mri.name])
                                        };
matlabbatch{1}.spm.util.imcalc.output = 'brainmask_0';
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = 'i1+i2+i3';
%matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

matlabbatch{2}.spm.spatial.smooth.data(1) = {fullfile(mri.pathname,'brainmask_0.nii')};
matlabbatch{2}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch{2}.spm.spatial.smooth.dtype = 0;
matlabbatch{2}.spm.spatial.smooth.im = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = 's';

matlabbatch{3}.spm.util.imcalc.input(1) = {fullfile(mri.pathname,'sbrainmask_0.nii')};
matlabbatch{3}.spm.util.imcalc.output = 'brainmask_1';
matlabbatch{3}.spm.util.imcalc.outdir = {''};
matlabbatch{3}.spm.util.imcalc.expression = 'i1>0.5';
matlabbatch{3}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{3}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{3}.spm.util.imcalc.options.mask = 0;
matlabbatch{3}.spm.util.imcalc.options.interp = 1;
matlabbatch{3}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',matlabbatch);
clear matlabbatch

%% coreg PET and SPECT to rigidly normalized MRI and mask it 
disp('Coregistering data to rigidly normalized MRI and masking it...')
%matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {fullfile(mri.pathname,['' mri.name])}; %fix for bad mri?

ictal.rfullfile = fullfile(ictal.path,['r' ictal.name]);
interictal.rfullfile = fullfile(interictal.path,['r' interictal.name]);
ictal.mrfullfile = fullfile(ictal.path,['mr' ictal.name]);
interictal.mrfullfile = fullfile(interictal.path,['mr' interictal.name]);

matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {fullfile(mri.pathname,['r' mri.name])};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {ictal.fullfile};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
spm_jobman('run',matlabbatch);
clear matlabbatch
 

matlabbatch{1}.spm.util.imcalc.input = {
                                       ictal.rfullfile
                                       fullfile(mri.pathname,'brainmask_1.nii')
                                        };
matlabbatch{1}.spm.util.imcalc.output = ictal.mrfullfile;
matlabbatch{1}.spm.util.imcalc.outdir = {ictal.path};
matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2';
%matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',matlabbatch);
clear matlabbatch


matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {fullfile(mri.pathname,['r' mri.name])};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {interictal.fullfile};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
spm_jobman('run',matlabbatch);
clear matlabbatch
 

matlabbatch{1}.spm.util.imcalc.input = {
                                       interictal.rfullfile
                                       fullfile(mri.pathname,'brainmask_1.nii')
                                        };
matlabbatch{1}.spm.util.imcalc.output = interictal.mrfullfile;
matlabbatch{1}.spm.util.imcalc.outdir = {interictal.path};
matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2';
%matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',matlabbatch);
clear matlabbatch

%% Global normalization
disp('Normalizing SPECT count ...')



%ictal
ictal.hdr = spm_vol(ictal.mrfullfile);
ictal.vol = spm_read_vols(ictal.hdr);
ictal.nvol = ictal.vol;

ictal.nvol=ictal.vol(ictal.vol>0);

rThr = mean(ictal.nvol(:));

ictal.rcvol = ictal.vol ./ rThr .* 50;

ictal.rchdr = ictal.hdr;
ictal.rchdr.fname = fullfile(ictal.path,['cmr' ictal.name]);
ictal.rchdr.private.dat.fname = ictal.rchdr.fname;
ictal.rchdr.private.dat.dtype = 'FLOAT32-LE'; %'INT16-LE', 'FLOAT32-LE'
ictal.rchdr.dt = [16 0]; % 4=16-bit integer; 16=32-bit real datatype

spm_write_vol(ictal.rchdr,ictal.rcvol);

clear rThr

%interictal
interictal.hdr = spm_vol(interictal.mrfullfile);
interictal.vol = spm_read_vols(interictal.hdr);
interictal.nvol = interictal.vol;
iinterictalctal.nvol=interictal.vol(ictal.vol>0);
rThr = mean(interictal.nvol(:));
interictal.rcvol = interictal.vol ./ rThr .* 50;

interictal.rchdr = interictal.hdr;
interictal.rchdr.fname = fullfile(interictal.path,['cmr' interictal.name]);
interictal.rchdr.private.dat.fname = interictal.rchdr.fname;
interictal.rchdr.private.dat.dtype = 'FLOAT32-LE'; %'INT16-LE', 'FLOAT32-LE'
interictal.rchdr.dt = [16 0]; % 4=16-bit integer; 16=32-bit real datatype

spm_write_vol(interictal.rchdr,interictal.rcvol);

clear rThr

%% Cleaning up

if DELETE_Temp_Files == 1;
    delete(fullfile(mri.pathname,['c*' mri.name]));
    delete(fullfile(mri.pathname,'*brainmask_0.nii'));   
    

else
end


end







































