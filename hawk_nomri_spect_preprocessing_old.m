function hawk_nomri_spect_preprocessing_v01
% When promted select spect images to be coregistered masked and count
% normalized
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

addpath('/Users/vlasta/Neuro/spm12')
spm_jobman('initcfg');
[ictal.name, ictal.path, ~] = uigetfile('*Ictal*.nii', 'Select Ictal file');
[interictal.name, interictal.path, ~] = uigetfile('*Baseline*.nii', 'Select Baseline file');
ictal.fullfile = fullfile(ictal.path,ictal.name);
interictal.fullfile = fullfile(interictal.path,interictal.name);

meanspect.oldfullfile=fullfile(ictal.path,['mean' ictal.name]);

meanspect.fullfile=fullfile(interictal.path,'meanspect.nii');
meanspect.path=interictal.path;
cd(interictal.path)

%% realign estimate and write and get mean image
realign{1}.spm.spatial.realign.estwrite.data = {{ictal.fullfile; interictal.fullfile}}';
realign{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
realign{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
realign{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 8;
realign{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
realign{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
realign{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
realign{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
realign{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
realign{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
realign{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
realign{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
realign{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
spm_jobman('run',realign);
clear jobs;

copyfile(meanspect.oldfullfile,meanspect.fullfile);
delete(meanspect.oldfullfile);

specttemplate.fullfile = fullfile(spm('Dir'),'templates','SPECT.nii');
specttemplate.cpfile = fullfile(ictal.path,'SPECT.nii') ;
brainmask.fullfile = fullfile(spm('Dir'),'apriori','brainmask.nii');
brainmask.cpfile = fullfile(ictal.path,'brainmask.nii');
copyfile(brainmask.fullfile,brainmask.cpfile);
copyfile(specttemplate.fullfile,specttemplate.cpfile);

matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {meanspect.fullfile};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {specttemplate.cpfile};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {
                                                  brainmask.cpfile
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

%% mask it

matlabbatch{1}.spm.util.imcalc.input = {
                                       fullfile(ictal.path,['r' ictal.name])
                                       fullfile(ictal.path,'rbrainmask.nii,1')
                                        };
matlabbatch{1}.spm.util.imcalc.output = fullfile(ictal.path,['mr' ictal.name]);
matlabbatch{1}.spm.util.imcalc.outdir = {ictal.path};
matlabbatch{1}.spm.util.imcalc.expression = 'i1.*(i2>0.5)';
%matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;


matlabbatch{2}.spm.util.imcalc.input = {
                                       fullfile(interictal.path,['r' interictal.name])
                                       fullfile(ictal.path,'rbrainmask.nii,1')
                                        };
matlabbatch{2}.spm.util.imcalc.output = fullfile(interictal.path,['mr' interictal.name]);
matlabbatch{2}.spm.util.imcalc.outdir = {interictal.path};
matlabbatch{2}.spm.util.imcalc.expression = 'i1.*(i2>0.5)';
matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{2}.spm.util.imcalc.options.mask = 0;
matlabbatch{2}.spm.util.imcalc.options.interp = 1;
matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',matlabbatch);
clear matlabbatch


%% Global normalization
disp('Normalizing SPECT count ...')


%ictal
ictal.hdr = spm_vol(fullfile(ictal.path,['mr' ictal.name]));
ictal.vol = spm_read_vols(ictal.hdr);
ictal.nanvol = ictal.vol;
ictal.nanvol(~ictal.nanvol)=NaN;
rThr = nanmean(ictal.nanvol(:));
ictal.rcvol = ictal.vol ./ rThr .* 50;

ictal.rchdr = ictal.hdr;
ictal.rchdr.fname = fullfile(ictal.path,['cmr' ictal.name]);
ictal.rchdr.private.dat.fname = ictal.rchdr.fname;
ictal.rchdr.private.dat.dtype = 'FLOAT32-LE'; %'INT16-LE', 'FLOAT32-LE'
ictal.rchdr.dt = [16 0]; % 4=16-bit integer; 16=32-bit real datatype

spm_write_vol(ictal.rchdr,ictal.rcvol);

clear rThr

%interictal
interictal.hdr = spm_vol(fullfile(interictal.path,['mr' interictal.name]));
interictal.vol = spm_read_vols(interictal.hdr);
interictal.nanvol = interictal.vol;
interictal.nanvol(~interictal.nanvol)=NaN;
rThr = nanmean(interictal.nanvol(:));
interictal.rcvol = interictal.vol ./ rThr .* 50;

interictal.rchdr = interictal.hdr;
interictal.rchdr.fname = fullfile(interictal.path,['cmr' interictal.name]);
interictal.rchdr.private.dat.fname = interictal.rchdr.fname;
interictal.rchdr.private.dat.dtype = 'FLOAT32-LE'; %'INT16-LE', 'FLOAT32-LE'
interictal.rchdr.dt = [16 0]; % 4=16-bit integer; 16=32-bit real datatype

spm_write_vol(interictal.rchdr,interictal.rcvol);

clear rThr

end


















