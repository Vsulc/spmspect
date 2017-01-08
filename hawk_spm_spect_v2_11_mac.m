function hawk_spm_spect_v2_11_mac
% Input spect images that are already coregistered masked and count
% normalized
% apriori and templates folders have to be copied from spm8
% isas_denominator.nii has to be present
% /Users/vlasta/Neuro/_togo/yale-controls by default
% if you don't have MRI just cancel when prompted for one, results are
% worse though
%
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

%% setting up
addpath('E:\Neuro\spm12')
controls_dir = 'E:\Neuro\_togo\yale-controls';
groupControlName = 'Yale';  %YALE HMPAO
FWHM_size = [16 16 16];
FWHM.dir = ['FWHM' num2str(FWHM_size(1))];
denominator.fullfile = fullfile(controls_dir,FWHM.dir,'isas_denominator.nii'); % path to the file

spm_jobman('initcfg');
         
DELETE_Temp_Files = 0;


  [ictal.file, work_dir, ictal.exists] = uigetfile('cmr*Iktal*.nii', 'Pick Ictal file');
  [interictal.file, ~, interictal.exists] = uigetfile('cmr*Baseline*.nii', 'Pick Baseline file');
  [mriT1.file, ~, mriT1.exists] = uigetfile('*rcoT1*.nii', 'Pick MRI file');
  
  brainmask_mri = fullfile(work_dir,'brainmask.nii');
  
if ictal.exists==0 || interictal.exists==0
    return
end

%     Ictal =  iktalPool(groupSize);
%     Interictal = baselinePool(groupSize);
%     Mri = mriPool(ceil(groupSize/2));

    [~, ictal.name , ictal.ext]  = fileparts(ictal.file); 
    [~, interictal.name , interictal.ext]  = fileparts(interictal.file);

    
 if mriT1.exists==0   
    meanimage.file = fullfile(work_dir,['meanspect' ictal.ext]);
    [~, meanimage.name , meanimage.ext]  = fileparts(meanimage.file);
 elseif mriT1.exists==2  || mriT1.exists==1
    [~, mriT1.name , mriT1.ext]  = fileparts(mriT1.file);
 end

    


%% clean it a bit
datum = fix(clock);
%temp_dir_name = ['SPM_SPECT_', groupControlName,'_FWHM',groupFWHMName];
temp_dir_name = sprintf('SPM_SPECT_%s_%d-%02d-%02d-%d-%d',groupControlName,datum(1:5));
temp_dir = [work_dir,filesep,temp_dir_name];
mkdir(temp_dir)
clear datum

ictal.cpfile = fullfile(temp_dir,[ictal.name ictal.ext]);
copyfile(ictal.file,ictal.cpfile)
if strcmpi(ictal.ext,'.img,1')
    copyfile([ictal.file(1:end-5) 'hdr'],[ictal.cpfile(1:end-5) 'hdr']);
end


interictal.cpfile = fullfile(temp_dir,[interictal.name interictal.ext]);
copyfile(interictal.file,interictal.cpfile)
if strcmpi(interictal.ext,'.img,1')
    copyfile([interictal.file(1:end-5) 'hdr'],[interictal.cpfile(1:end-5) 'hdr']);
end


if mriT1.exists==0
        meanimage.cpfile = fullfile(temp_dir,['meanspect' ictal.ext]);
        copyfile(meanimage.file,meanimage.cpfile)
            if strcmpi(meanimage.ext,'.img,1')
                copyfile([meanimage.file(1:end-5) 'hdr'],[meanimage.cpfile(1:end-5) 'hdr']);
            end
elseif mriT1.exists==2  || mriT1.exists==1
            mriT1.cpfile = fullfile(temp_dir,[mriT1.name mriT1.ext]);
            copyfile(mriT1.file,mriT1.cpfile)
            if strcmpi(mriT1.ext,'.img,1')
                copyfile([mriT1.file(1:end-5) 'hdr'],[mriT1.cpfile(1:end-5) 'hdr']);
            end
end
%% realign (r)

%% Real analysis starts here
ictal.hdr = spm_vol(ictal.cpfile);
ictal.vol = spm_read_vols(ictal.hdr);
interictal.hdr = spm_vol(interictal.cpfile);
interictal.vol = spm_read_vols(interictal.hdr);

%% subtract Ictal and Interictal image
disp('Creating difference image ...')

DiffVolume.hdr = spm_vol(ictal.cpfile);
DiffVolume.hdr.fname = fullfile(temp_dir,['DiffVolume' '.nii']);
DiffVolume.fullfile = fullfile(temp_dir,['DiffVolume' '.nii']);
DiffVolume.hdr.private.dat.fname = DiffVolume.hdr.fname;
DiffVolume.hdr.private.dat.dtype = 'FLOAT32-LE';
DiffVolume.hdr.dt = [16 0];
DiffVolume.hdr.descrip = 'Count normalized and subtracted';

DiffVolume.vol = ictal.vol - interictal.vol;
spm_write_vol(DiffVolume.hdr,DiffVolume.vol) ;



%% Normalize to standard space (w) old norm
disp('Normalizing to standard space ...')

if mriT1.exists==0
normalize2std{1}.spm.tools.oldnorm.estwrite.subj.source{1} = meanimage.cpfile;
normalize2std{1}.spm.tools.oldnorm.estwrite.subj.resample = {meanimage.cpfile;DiffVolume.fullfile;interictal.cpfile;ictal.cpfile};  
normalize2std{1}.spm.tools.oldnorm.estwrite.eoptions.template = {fullfile(spm('Dir'),'templates','SPECT.nii')};
disp('MRI does not exist, normalizing only using SPECT data ...')
elseif  mriT1.exists==2 || mriT1.exists==1
normalize2std{1}.spm.tools.oldnorm.estwrite.subj.source{1} = mriT1.cpfile;
normalize2std{1}.spm.tools.oldnorm.estwrite.subj.resample = {mriT1.cpfile;DiffVolume.fullfile;interictal.cpfile;ictal.cpfile};  
normalize2std{1}.spm.tools.oldnorm.estwrite.eoptions.template = {fullfile(spm('Dir'),'templates','T2.nii')}; %T1
disp('MRI normalization ...')
end
normalize2std{1}.spm.tools.oldnorm.estwrite.subj.wtsrc = '';
normalize2std{1}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
normalize2std{1}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
normalize2std{1}.spm.tools.oldnorm.estwrite.eoptions.smoref = 0;
normalize2std{1}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
normalize2std{1}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
normalize2std{1}.spm.tools.oldnorm.estwrite.eoptions.nits = 16;
normalize2std{1}.spm.tools.oldnorm.estwrite.eoptions.reg = 1;
normalize2std{1}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
normalize2std{1}.spm.tools.oldnorm.estwrite.roptions.bb = [Inf Inf Inf; Inf Inf Inf]; %% ? need to change [[-90,-126,-72];[90,90,108]]
normalize2std{1}.spm.tools.oldnorm.estwrite.roptions.vox = [2 2 2];
normalize2std{1}.spm.tools.oldnorm.estwrite.roptions.interp = 1;
normalize2std{1}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
normalize2std{1}.spm.tools.oldnorm.estwrite.roptions.prefix = 'w';


spm_jobman('run',normalize2std);
clear normalize2std;

%% Smooth with 16 mm FWHM %% siscom smooth with 12 mm
disp('Smoothing ...')

wDiffVolume = fullfile(temp_dir,['wDiffVolume' '.nii']);
wcmeanimage = fullfile(temp_dir,['wmeanspect'  '.nii']);
interictal.wcrmname = fullfile(temp_dir,['w' interictal.name '.nii']);
ictal.wcrmname = fullfile(temp_dir,['w' ictal.name '.nii']);

jobs{1}.spm.spatial.smooth.data = {wDiffVolume};
jobs{1}.spm.spatial.smooth.fwhm = FWHM_size;
jobs{1}.spm.spatial.smooth.dtype = 0;
jobs{1}.spm.spatial.smooth.im = 0;
jobs{1}.spm.spatial.smooth.prefix = 's';

jobs{2}.spm.spatial.smooth.data = {interictal.wcrmname};
jobs{2}.spm.spatial.smooth.fwhm = FWHM_size;
jobs{2}.spm.spatial.smooth.dtype = 0;
jobs{2}.spm.spatial.smooth.im = 0;
jobs{2}.spm.spatial.smooth.prefix = 's';

jobs{3}.spm.spatial.smooth.data = {ictal.wcrmname};
jobs{3}.spm.spatial.smooth.fwhm = FWHM_size;
jobs{3}.spm.spatial.smooth.dtype = 0;
jobs{3}.spm.spatial.smooth.im = 0;
jobs{3}.spm.spatial.smooth.prefix = 's';

if mriT1.exists==0
        jobs{4}.spm.spatial.smooth.data = {wcmeanimage};
        jobs{4}.spm.spatial.smooth.fwhm = FWHM_size;
        jobs{4}.spm.spatial.smooth.dtype = 0;
        jobs{4}.spm.spatial.smooth.im = 0;
        jobs{4}.spm.spatial.smooth.prefix = 's';
    elseif mriT1.exists==1  
end

spm_jobman('run',jobs);
clear jobs;


%% Mask so it matches with the mask used for control group analysis
disp('Masking ...')

swcmeanimage = fullfile(temp_dir,['swmeanspect' '.nii']);
swDiffVolume = fullfile(temp_dir,['swDiffVolume' '.nii']);
interictal.swcrmfile = fullfile(temp_dir,['sw' interictal.name '.nii']);
ictal.swcrmfile = fullfile(temp_dir,['sw' ictal.name '.nii']);
brainmask.path = fullfile(spm('Dir'),'apriori','brainmask.nii');


jobs{1}.spm.util.imcalc.input = {swDiffVolume;brainmask.path};
jobs{1}.spm.util.imcalc.output = 'mswDiffVolume.nii';
jobs{1}.spm.util.imcalc.outdir = {temp_dir}; %was ''
jobs{1}.spm.util.imcalc.expression = 'i1.*(i2>.5)';
jobs{1}.spm.util.imcalc.options.dmtx = 0;
jobs{1}.spm.util.imcalc.options.mask = 0;
jobs{1}.spm.util.imcalc.options.interp = 1;
jobs{1}.spm.util.imcalc.options.dtype = 4;



jobs{2}.spm.util.imcalc.input = {interictal.swcrmfile;brainmask.path};
jobs{2}.spm.util.imcalc.output = ['msw' interictal.name '.nii'];
jobs{2}.spm.util.imcalc.outdir = {temp_dir}; %was ''
jobs{2}.spm.util.imcalc.expression = 'i1.*(i2>.5)';
jobs{2}.spm.util.imcalc.options.dmtx = 0;
jobs{2}.spm.util.imcalc.options.mask = 0;
jobs{2}.spm.util.imcalc.options.interp = 1;
jobs{2}.spm.util.imcalc.options.dtype = 4;


jobs{3}.spm.util.imcalc.input = {ictal.swcrmfile;brainmask.path};
jobs{3}.spm.util.imcalc.output = ['msw' ictal.name '.nii'];
jobs{3}.spm.util.imcalc.outdir = {temp_dir}; %was ''
jobs{3}.spm.util.imcalc.expression = 'i1.*(i2>.5)';
jobs{3}.spm.util.imcalc.options.dmtx = 0;
jobs{3}.spm.util.imcalc.options.mask = 0;
jobs{3}.spm.util.imcalc.options.interp = 1;
jobs{3}.spm.util.imcalc.options.dtype = 4;

if mriT1.exists==0
        jobs{4}.spm.util.imcalc.input = {swcmeanimage;brainmask.path};
        jobs{4}.spm.util.imcalc.output = ['mswmeanspect''.nii'];
        jobs{4}.spm.util.imcalc.outdir = {temp_dir}; %was ''
        jobs{4}.spm.util.imcalc.expression = 'i1.*(i2>.5)';
        jobs{4}.spm.util.imcalc.options.dmtx = 0;
        jobs{4}.spm.util.imcalc.options.mask = 0;
        jobs{4}.spm.util.imcalc.options.interp = 1;
        jobs{4}.spm.util.imcalc.options.dtype = 4;
elseif mriT1.exists==1  
end

spm_jobman('run',jobs);
clear jobs;

%% STATS 

%% Stats ISAS
disp('Running ISAS ...')
denominator.hdr = spm_vol(denominator.fullfile);
mswDiff.hdr = spm_vol(fullfile(temp_dir,['mswDiffVolume' '.nii']));
mswDiff.vol = spm_read_vols(mswDiff.hdr);
%denominator = spm_read_vols(spm_vol(fullfile(controls_dir,'controls-diffvolumes-masked','denominator.nii'))); % this is not the right one
%denominator = spm_read_vols(spm_vol(fullfile(controls_dir,'controls-diffvolumes-masked','isas_denominator.nii')));
denominator.vol = spm_read_vols(denominator.hdr);
RawISAS.vol = mswDiff.vol./denominator.vol;
RawISAS.hdr = mswDiff.hdr;
RawISAS.hdr.fname = fullfile(temp_dir,'SPMSPECT_raw_results.nii');

% T.vol =
% diff./(squeeze((sqrt(pi)/2)*sum(abs(alldiff_vol))/size(alldiff_vol,1))/sqrt(size(alldiff_vol,1))
RawISAS.hdr.dt=[16 0];
RawISAS.vol(isnan(RawISAS.vol))=0;




%% set threshold % seems OK
% if COMPESATE_SD_range==1 % to make it look more like a standard SISCOM images
    spm_write_vol(RawISAS.hdr,RawISAS.vol);
    RawISAS.stdev = RawISAS.vol(:);
    RawISAS.stdev(RawISAS.stdev == 0) = NaN;
    RawISAS.vol = RawISAS.vol./nanstd(RawISAS.stdev(:));
%     isasHyper.msk = (RawISAS.vol>1.5);
%     isasHypo.msk = (-RawISAS.vol>1.5);
    RawISAS.hdr.fname = fullfile(temp_dir,'SPMSPECT_SDcorr.nii');
    spm_write_vol(RawISAS.hdr,RawISAS.vol);
% else
%     isasHyper.msk = (RawISAS.vol>2);
%     isasHypo.msk = (-RawISAS.vol>2);
%     spm_write_vol(RawISAS.hdr,RawISAS.vol);
% end

% apply cluster threshold of 125 (5x5x5) voxels with 2x2x2 mm voxels ( = 1
% cm3 )
% isasHyper.hdr = mswDiff.hdr;
% isasHyper.tmsk = bwareaopen(isasHyper.msk,125,26);
% isasHyper.vol = RawISAS.vol .* isasHyper.tmsk;
% isasHyper.hdr.fname = fullfile(temp_dir,'ISAS_hyperperfusion.nii');
% spm_write_vol(isasHyper.hdr,isasHyper.vol);
% 
% isasHypo.hdr = mswDiff.hdr;
% isasHypo.tmsk = bwareaopen(isasHypo.msk,125,26);
% isasHypo.vol = -RawISAS.vol .* isasHypo.tmsk;
% isasHypo.hdr.fname = fullfile(temp_dir,'ISAS_hypoperfusion.nii');
% spm_write_vol(isasHypo.hdr,isasHypo.vol);
% time.isas = toc;

%% AND now in MRI SPACE

siscom.rawhdr = spm_vol(fullfile(temp_dir,'DiffVolume.nii'));
siscom.rawvol = spm_read_vols(siscom.rawhdr);
siscom.std  = std(siscom.rawvol(siscom.rawvol~=0));
siscom.hdr = siscom.rawhdr;

siscom.hyper = (siscom.rawvol/siscom.std).*(siscom.rawvol>(2*siscom.std));
siscom.hdr.fname = fullfile(temp_dir,'siscom_hyper_2dot0.nii');
spm_write_vol(siscom.hdr,siscom.hyper);

siscom.hypo = (-siscom.rawvol/siscom.std).*(siscom.rawvol<(-2*siscom.std));
siscom.hdr.fname = fullfile(temp_dir,'siscom_hypo_2dot0.nii');
spm_write_vol(siscom.hdr,siscom.hypo);

siscom.hyper = (siscom.rawvol/siscom.std).*(siscom.rawvol>(1.5*siscom.std));
siscom.hdr.fname = fullfile(temp_dir,'siscom_hyper_1dot5.nii');
spm_write_vol(siscom.hdr,siscom.hyper);

siscom.hypo = (-siscom.rawvol/siscom.std).*(siscom.rawvol<(-1.5*siscom.std));
siscom.hdr.fname = fullfile(temp_dir,'siscom_hypo_1dot5.nii');
spm_write_vol(siscom.hdr,siscom.hypo);

%% Backnormalize
% from template space to meanspect (ictal) cmean_sn.mat
% from ictal spect space to MRI
if mriT1.exists==0
backnorm{1}.spm.util.defs.comp{1}.sn2def.matname = {fullfile(temp_dir,'meanspect_sn.mat')};
elseif mriT1.exists==2 || mriT1.exists==1
backnorm{1}.spm.util.defs.comp{1}.sn2def.matname = {fullfile(temp_dir,[mriT1.name '_sn.mat'])};    
end
backnorm{1}.spm.util.defs.comp{1}.sn2def.vox = [NaN NaN NaN];
backnorm{1}.spm.util.defs.comp{1}.sn2def.bb = [NaN NaN NaN; NaN NaN NaN];
backnorm{1}.spm.util.defs.out{1}.push.fnames = {
													fullfile(temp_dir,'SPMSPECT_SDcorr.nii')
                                                    fullfile(temp_dir,['msw' ictal.name '.nii'])
                                                    fullfile(temp_dir,['msw' interictal.name '.nii'])
                                                    fullfile(temp_dir,'wT1.nii')
                                                   };
backnorm{1}.spm.util.defs.out{1}.push.weight = {''};
backnorm{1}.spm.util.defs.out{1}.push.savedir.saveusr = {temp_dir};
if mriT1.exists==0
backnorm{1}.spm.util.defs.out{1}.push.fov.file = {fullfile(temp_dir,'meanspect.nii')};
elseif mriT1.exists==2 || mriT1.exists==1
backnorm{1}.spm.util.defs.out{1}.push.fov.file = {fullfile(temp_dir,[mriT1.name mriT1.ext])};    
end
backnorm{1}.spm.util.defs.out{1}.push.preserve = 0;
backnorm{1}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
backnorm{1}.spm.util.defs.out{1}.push.prefix = '_Pts_';

spm_jobman('run',backnorm);
clear backnorm


%% DELETE temp files

if DELETE_Temp_Files == 1;
%    delete([temp_dir filesep 'c*.nii'])
    delete([temp_dir filesep 'sw*.nii'])
    delete([temp_dir filesep 'w*.nii'])
%    delete([temp_dir filesep 'wc*.nii'])
%    delete([temp_dir filesep 'wm*.nii'])
    delete([temp_dir filesep 'ms*.nii'])
    delete([temp_dir filesep '*mean*'])
    delete([temp_dir filesep 'nsiscom_*.nii'])
%    delete([temp_dir filesep 'Statiscom_*.nii'])
%    delete([temp_dir filesep 'ISAS_*.nii'])
    delete([temp_dir filesep '*r' ictal.name '.nii'])
    delete([temp_dir filesep '*r' interictal.name '.nii'])
%    delete([temp_dir filesep 'wDiffVolume.nii'])
    delete([temp_dir filesep 'rp*.txt'])
else
end

% close all
clear ALL
disp('=-=-=-=-=-==-=-=-=-= READY FOR THE NEXT ONE =-=-=-=-=-==-=-=-=-=')
end












