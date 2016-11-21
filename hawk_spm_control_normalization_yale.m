function hawk_spm_control_normalization_yale
%% Preprocessing of SPECT data for spmspect
% http://spect.yale.edu/downloads.html
% templates and apriori has to be coppied from SPM 8 dir
% to work with spm8 and 12 spect data has to be preprocessed using the
% hawk_spm_fix_yale_controls script
% 7 and 10 are not well segmented

addpath('/Users/vlasta/Neuro/spm12')
spm_jobman('initcfg');

FWHM.size = [16 16 16];
controls.dir = '/Users/vlasta/Neuro/controls_preprocessing/Raw_Healthy_Normals';
controls.ext = '.nii';
controls.prefix = 'f';

ictalpool =   {
                    fullfile(controls.dir,[controls.prefix 'HN001_D1' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN002_D1' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN003_D1' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN004_D1' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN005_D1' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN006_D1' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN007_D1' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN008_D1' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN009_D1' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN010_D1' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN011_D1' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN012_D1' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN013_D1' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN014_D1' controls.ext])  
                  };


interictalpool =   {
                    fullfile(controls.dir,[controls.prefix 'HN001_D2' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN002_D2' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN003_D2' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN004_D2' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN005_D2' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN006_D2' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN007_D2' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN008_D2' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN009_D2' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN010_D2' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN011_D2' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN012_D2' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN013_D2' controls.ext])  
                    fullfile(controls.dir,[controls.prefix 'HN014_D2' controls.ext])  
                  };




for r = 1:size(interictalpool,1)
Ictal = ictalpool{r};
Interictal = interictalpool{r};

[pth, name , ext]  = fileparts(Ictal);
[~, intername , interext]  = fileparts(Interictal);


% if exist('/Users/vlasta/Neuro/spm8/apriori/brainmask.nii','file');
%     BrainMask = spm_vol('/Users/vlasta/Neuro/spm8/apriori/brainmask.nii','file');
% else
%     BrainMask = spm_vol(spm_get(1,'IMAGE', 'Select Brain Mask'));
% end

%% coregister (r)

tic
% jobs{1}.spm.spatial.realign.estwrite.data{1}{1} = Ictal;
% jobs{1}.spm.spatial.realign.estwrite.data{1}{2} = Interictal;
jobs{1}.spm.spatial.realign.estwrite.data = {{Ictal; Interictal}}';
jobs{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
jobs{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
jobs{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5; %% 7 in ISAS
jobs{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
jobs{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
jobs{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
jobs{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
jobs{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];

jobs{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
jobs{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
jobs{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
jobs{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

spm_jobman('run',jobs);
clear jobs;
time.coregister_e(r) = toc;

rIctal = fullfile(pth,['r' name ext]); %fullfile(pth,['r' name ext]);
rInterictal = fullfile(pth,['r' intername interext]);
meanimage = fullfile(pth,['mean' name ext]);

%% Segment 
% spm 5 style segment to get brain masks
% (spm 8 segmentation takes forever - 52 minutes)
tic

fastseg{1}.spm.tools.oldseg.data = {meanimage};
fastseg{1}.spm.tools.oldseg.output.GM = [0 0 1];
fastseg{1}.spm.tools.oldseg.output.WM = [0 0 1];
fastseg{1}.spm.tools.oldseg.output.CSF = [0 0 1];
fastseg{1}.spm.tools.oldseg.output.biascor = 1;
fastseg{1}.spm.tools.oldseg.output.cleanup = 2;
fastseg{1}.spm.tools.oldseg.opts.tpm = {
                                               fullfile(spm('Dir'),'toolbox','OldSeg','grey.nii')
                                               fullfile(spm('Dir'),'toolbox','OldSeg','white.nii')
                                               fullfile(spm('Dir'),'toolbox','OldSeg','csf.nii')
                                               };
fastseg{1}.spm.tools.oldseg.opts.ngaus = [2
                                                 2
                                                 2
                                                 4];
fastseg{1}.spm.tools.oldseg.opts.regtype = 'mni';
fastseg{1}.spm.tools.oldseg.opts.warpreg = 1;
fastseg{1}.spm.tools.oldseg.opts.warpco = 25;
fastseg{1}.spm.tools.oldseg.opts.biasreg = 0.01;
fastseg{1}.spm.tools.oldseg.opts.biasfwhm = 60;
fastseg{1}.spm.tools.oldseg.opts.samp = 3;
fastseg{1}.spm.tools.oldseg.opts.msk = {''};

spm_jobman('run',fastseg);

% clear fastseg
% 
% matlabbatch{1}.spm.spatial.preproc.channel.vols = {meanimage};
% matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
% matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
% matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
% matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,1')}; 
% matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
% matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
% matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
% matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,2')};
% matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
% matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
% matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
% matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,3')};
% matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
% matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
% matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
% matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,4')};
% matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
% matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
% matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
% matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,5')};
% matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
% matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
% matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
% matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,6')};
% matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
% matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
% matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
% matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
% matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
% matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
% matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
% matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
% matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
% matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
% 
% spm_jobman('run',matlabbatch);
% clear matlabbatch

time.fastseg(r) = toc;

%% Brainmasking
tic
% 
% c1mean.hdr = spm_vol(fullfile(pth,['c1mean' name ext]));
% c1mean.vol = spm_read_vols(c1mean.hdr);
% c2mean.hdr = spm_vol(fullfile(pth,['c2mean' name ext]));
% c2mean.vol = spm_read_vols(c2mean.hdr);
% c3mean.hdr = spm_vol(fullfile(pth,['c3mean' name ext]));
% c3mean.vol = spm_read_vols(c3mean.hdr);
% brainmask.vol = (c1mean.vol + c2mean.vol + c3mean.vol)>0.5;
% brainmask.vol(brainmask.vol>0)=1; %new
% brainmask.hdr = c1mean.hdr;
% brainmask.hdr.fname = fullfile(pth,['brainmask_0' name '.nii']);
% spm_write_vol(brainmask.hdr,brainmask.vol);

matlabbatch{1}.spm.util.imcalc.input = {
                                       fullfile(pth,['c1mean' name ext])
                                       fullfile(pth,['c2mean' name ext])
                                       fullfile(pth,['c3mean' name ext])
                                        };
matlabbatch{1}.spm.util.imcalc.output = ['brainmask_0' name '.nii'];
matlabbatch{1}.spm.util.imcalc.outdir = {pth};
matlabbatch{1}.spm.util.imcalc.expression = 'i1+i2+i3';
%matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

matlabbatch{2}.spm.spatial.smooth.data(1) = {fullfile(pth,['brainmask_0' name '.nii'])};
matlabbatch{2}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch{2}.spm.spatial.smooth.dtype = 0;
matlabbatch{2}.spm.spatial.smooth.im = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = 's';

matlabbatch{3}.spm.util.imcalc.input(1) = {fullfile(pth,['sbrainmask_0' name '.nii'])};
matlabbatch{3}.spm.util.imcalc.output = ['brainmask_1' name '.nii'];
matlabbatch{3}.spm.util.imcalc.outdir = {pth};  
matlabbatch{3}.spm.util.imcalc.expression = 'i1>0.5';
matlabbatch{3}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{3}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{3}.spm.util.imcalc.options.mask = 0;
matlabbatch{3}.spm.util.imcalc.options.interp = 1;
matlabbatch{3}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',matlabbatch);
clear matlabbatch

brainmask.hdr = spm_vol(fullfile(pth,['brainmask_1' name '.nii']));
brainmask.vol = spm_read_vols(brainmask.hdr);

mIctal.hdr = spm_vol(rIctal);
mIctal.vol = spm_read_vols(mIctal.hdr);
mIctal.vol = mIctal.vol  .* brainmask.vol;
mIctal.vol(mIctal.vol == 0) = NaN; %new
mIctal.hdr.fname = fullfile(pth,['mr' name '.nii']); %new
spm_write_vol(mIctal.hdr,mIctal.vol);


mInterictal.hdr = spm_vol(rInterictal);
mInterictal.vol = spm_read_vols(mInterictal.hdr);
mInterictal.vol = mInterictal.vol  .* brainmask.vol;
mInterictal.vol(mInterictal.vol == 0) = NaN;%new
mInterictal.hdr.fname = fullfile(pth,['mr' intername '.nii']);%new
spm_write_vol(mInterictal.hdr,mInterictal.vol);

mmeanimage.hdr = spm_vol(meanimage);
mmeanimage.vol = spm_read_vols(mmeanimage.hdr);
mmeanimage.vol = mmeanimage.vol  .* brainmask.vol;
mmeanimage.vol(mmeanimage.vol == 0) = NaN;%new
mmeanimage.hdr.fname = fullfile(pth,['mmean' name '.nii']);%new
spm_write_vol(mmeanimage.hdr,mmeanimage.vol);

time.brainmasking(r) = toc;

%% better brainmask with smoothing


%% normalize SPECT count (c)

disp('Normalizing SPECT count ...')
tic

rThr = nanmean(mIctal.vol(:));
rcIctal.vol = mIctal.vol ./ rThr .* 50;

rcIctal.hdr = mIctal.hdr;
rcIctal.hdr.fname = fullfile(pth,['cmr' name '.nii']);
rcIctal.hdr.private.dat.fname = rcIctal.hdr.fname;
rcIctal.hdr.private.dat.dtype = 'FLOAT32-LE';
rcIctal.hdr.dt = [16 0];

spm_write_vol(rcIctal.hdr,rcIctal.vol);

clear rThr

rThr = nanmean(mInterictal.vol(:));
rcInterictal.vol = mInterictal.vol ./ rThr .* 50;

rcInterictal.hdr = mInterictal.hdr;
rcInterictal.hdr.fname = fullfile(pth,['cmr' intername  '.nii']);
rcInterictal.hdr.private.dat.fname = rcInterictal.hdr.fname;
rcInterictal.hdr.private.dat.dtype = 'FLOAT32-LE';
rcInterictal.hdr.dt = [16 0];

spm_write_vol(rcInterictal.hdr,rcInterictal.vol);

clear rThr

rThr = nanmean(mmeanimage.vol(:));
rcmeanimage.vol = mmeanimage.vol ./ rThr .* 50;

rcmeanimage.hdr = mmeanimage.hdr;
rcmeanimage.hdr.fname = fullfile(pth,['cmmean' name '.nii']);
rcmeanimage.hdr.private.dat.fname = rcmeanimage.hdr.fname;
rcmeanimage.hdr.private.dat.dtype = 'FLOAT32-LE';
rcmeanimage.hdr.dt = [16 0];

spm_write_vol(rcmeanimage.hdr,rcmeanimage.vol);

clear rThr

time.cnormalize(r) = toc;

%% subtract Ictal and Interictal image
disp('Creating difference image ...')
tic

DiffVolume.hdr = spm_vol(rIctal);
DiffVolume.hdr.fname = fullfile(pth,['DiffVolume' '_' name '-' intername '.nii']);
DiffVolume.hdr.private.dat.fname = DiffVolume.hdr.fname;
DiffVolume.hdr.private.dat.dtype = 'FLOAT32-LE';
DiffVolume.hdr.dt = [16 0];
DiffVolume.hdr.descrip = 'Count normalized and subtracted';

DiffVolume.vol = rcIctal.vol - rcInterictal.vol;
spm_write_vol(DiffVolume.hdr,DiffVolume.vol) ;

time.subtract(r) = toc;

%% Normalize to standard space (w) spm.tools.oldnorm.estwrite
disp('Normalizing clinical data ...')
tic

jobs{1}.spm.tools.oldnorm.estwrite.subj.source{1} = rcmeanimage.hdr.fname;
jobs{1}.spm.tools.oldnorm.estwrite.subj.wtsrc = '';
jobs{1}.spm.tools.oldnorm.estwrite.subj.resample = {rcmeanimage.hdr.fname;DiffVolume.hdr.fname;rcInterictal.hdr.fname;rcIctal.hdr.fname};  

jobs{1}.spm.tools.oldnorm.estwrite.eoptions.template ={fullfile(spm('Dir'),'templates','SPECT.nii')};
jobs{1}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
jobs{1}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
jobs{1}.spm.tools.oldnorm.estwrite.eoptions.smoref = 0;
jobs{1}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
jobs{1}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
jobs{1}.spm.tools.oldnorm.estwrite.eoptions.nits = 16;
jobs{1}.spm.tools.oldnorm.estwrite.eoptions.reg = 1;

jobs{1}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
jobs{1}.spm.tools.oldnorm.estwrite.roptions.bb = [Inf Inf Inf; Inf Inf Inf]; %% ? need to change [[-90,-126,-72];[90,90,108]]
jobs{1}.spm.tools.oldnorm.estwrite.roptions.vox = [2 2 2];
jobs{1}.spm.tools.oldnorm.estwrite.roptions.interp = 1;
jobs{1}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
jobs{1}.spm.tools.oldnorm.estwrite.roptions.prefix = 'w';
spm_jobman('run',jobs);
clear jobs;

time.normalize2std(r) = toc;

%% Smooth with 12 mm FWHM % originally 16
disp('Smoothing ...')
tic

wDiffVolume = fullfile(pth,['wDiffVolume' '_' name '-' intername '.nii']);
wcmeanimage = fullfile(pth,['wcmmean' name '.nii']);
wcrInterictal = fullfile(pth,['wcmr' intername '.nii']);
wcIctal = fullfile(pth,['wcmr' name '.nii']);

jobs{1}.spm.spatial.smooth.data = {wDiffVolume};
jobs{1}.spm.spatial.smooth.fwhm = FWHM.size;
jobs{1}.spm.spatial.smooth.dtype = 0;
jobs{1}.spm.spatial.smooth.im = 0;
jobs{1}.spm.spatial.smooth.prefix = 's';

jobs{2}.spm.spatial.smooth.data = {wcmeanimage};
jobs{2}.spm.spatial.smooth.fwhm = FWHM.size;
jobs{2}.spm.spatial.smooth.dtype = 0;
jobs{2}.spm.spatial.smooth.im = 0;
jobs{2}.spm.spatial.smooth.prefix = 's';

jobs{3}.spm.spatial.smooth.data = {wcrInterictal};
jobs{3}.spm.spatial.smooth.fwhm = FWHM.size;
jobs{3}.spm.spatial.smooth.dtype = 0;
jobs{3}.spm.spatial.smooth.im = 0;
jobs{3}.spm.spatial.smooth.prefix = 's';

jobs{4}.spm.spatial.smooth.data = {wcIctal};
jobs{4}.spm.spatial.smooth.fwhm = FWHM.size;
jobs{4}.spm.spatial.smooth.dtype = 0;
jobs{4}.spm.spatial.smooth.im = 0;
jobs{4}.spm.spatial.smooth.prefix = 's';

spm_jobman('run',jobs);
clear jobs;
time.smooth(r) = toc;

%% Mask
disp('Masking ...')
tic

swcmeanimage = fullfile(pth,['swcmmean' name '.nii']);
swDiffVolume = fullfile(pth,['swDiffVolume' '_' name '-' intername '.nii']);
swcrInterictal = fullfile(pth,['swcmr' intername '.nii']);
swcrIctal = fullfile(pth,['swcmr' name '.nii']);
brainmask.path = fullfile(spm('Dir'),'apriori','brainmask.nii');


jobs{1}.spm.util.imcalc.input = {swDiffVolume;brainmask.path};
jobs{1}.spm.util.imcalc.output = ['mswDiffVolume' '_' name '-' intername '.nii'];
jobs{1}.spm.util.imcalc.outdir = {pth}; %was ''
jobs{1}.spm.util.imcalc.expression = 'i1.*(i2>.5)';
jobs{1}.spm.util.imcalc.options.dmtx = 0;
jobs{1}.spm.util.imcalc.options.mask = 0;
jobs{1}.spm.util.imcalc.options.interp = 1;
jobs{1}.spm.util.imcalc.options.dtype = 4;


jobs{2}.spm.util.imcalc.input = {swcmeanimage;brainmask.path};
jobs{2}.spm.util.imcalc.output = ['mswcmmean' name  '.nii'];
jobs{2}.spm.util.imcalc.outdir = {pth}; %was ''
jobs{2}.spm.util.imcalc.expression = 'i1.*(i2>.5)';
jobs{2}.spm.util.imcalc.options.dmtx = 0;
jobs{2}.spm.util.imcalc.options.mask = 0;
jobs{2}.spm.util.imcalc.options.interp = 1;
jobs{2}.spm.util.imcalc.options.dtype = 4;

jobs{3}.spm.util.imcalc.input = {swcrInterictal;brainmask.path};
jobs{3}.spm.util.imcalc.output = ['mswcmr' intername  '.nii'];
jobs{3}.spm.util.imcalc.outdir = {pth}; %was ''
jobs{3}.spm.util.imcalc.expression = 'i1.*(i2>.5)';
jobs{3}.spm.util.imcalc.options.dmtx = 0;
jobs{3}.spm.util.imcalc.options.mask = 0;
jobs{3}.spm.util.imcalc.options.interp = 1;
jobs{3}.spm.util.imcalc.options.dtype = 4;


jobs{4}.spm.util.imcalc.input = {swcrIctal;brainmask.path};
jobs{4}.spm.util.imcalc.output = ['mswcmr' name  '.nii'];
jobs{4}.spm.util.imcalc.outdir = {pth}; %was ''
jobs{4}.spm.util.imcalc.expression = 'i1.*(i2>.5)';
jobs{4}.spm.util.imcalc.options.dmtx = 0;
jobs{4}.spm.util.imcalc.options.mask = 0;
jobs{4}.spm.util.imcalc.options.interp = 1;
jobs{4}.spm.util.imcalc.options.dtype = 4;

spm_jobman('run',jobs);
clear jobs;
time.mask(r) = toc;
end


tic
controlpool =   {
                    fullfile(controls.dir,['mswDiffVolume' '_' controls.prefix 'HN001_D1' '-' controls.prefix 'HN001_D2' '.nii'])  
                    fullfile(controls.dir,['mswDiffVolume' '_' controls.prefix 'HN002_D1' '-' controls.prefix 'HN002_D2' '.nii'])  
                    fullfile(controls.dir,['mswDiffVolume' '_' controls.prefix 'HN003_D1' '-' controls.prefix 'HN003_D2' '.nii'])  
                    fullfile(controls.dir,['mswDiffVolume' '_' controls.prefix 'HN004_D1' '-' controls.prefix 'HN004_D2' '.nii'])  
                    fullfile(controls.dir,['mswDiffVolume' '_' controls.prefix 'HN005_D1' '-' controls.prefix 'HN005_D2' '.nii'])  
                    fullfile(controls.dir,['mswDiffVolume' '_' controls.prefix 'HN006_D1' '-' controls.prefix 'HN006_D2' '.nii'])  
                    fullfile(controls.dir,['mswDiffVolume' '_' controls.prefix 'HN007_D1' '-' controls.prefix 'HN007_D2' '.nii'])  
                    fullfile(controls.dir,['mswDiffVolume' '_' controls.prefix 'HN008_D1' '-' controls.prefix 'HN008_D2' '.nii'])  
                    fullfile(controls.dir,['mswDiffVolume' '_' controls.prefix 'HN009_D1' '-' controls.prefix 'HN009_D2' '.nii'])  
                    fullfile(controls.dir,['mswDiffVolume' '_' controls.prefix 'HN010_D1' '-' controls.prefix 'HN010_D2' '.nii'])  
                    fullfile(controls.dir,['mswDiffVolume' '_' controls.prefix 'HN011_D1' '-' controls.prefix 'HN011_D2' '.nii'])  
                    fullfile(controls.dir,['mswDiffVolume' '_' controls.prefix 'HN012_D1' '-' controls.prefix 'HN012_D2' '.nii'])  
                    fullfile(controls.dir,['mswDiffVolume' '_' controls.prefix 'HN013_D1' '-' controls.prefix 'HN013_D2' '.nii'])  
                    fullfile(controls.dir,['mswDiffVolume' '_' controls.prefix 'HN014_D1' '-' controls.prefix 'HN014_D2' '.nii'])   
                  };
alldiff.vol = zeros(size(controlpool,1),91,109,91); 
for c = 1:size(controlpool,1)
    
alldiff.vol(c,:,:,:) = spm_read_vols(spm_vol(controlpool{c}));

end

% I guess this is to look what part of brain is most susceptible to
% perfusion changes

avdiff.hdr = spm_vol(controlpool{1});
avdiff.hdr.fname = fullfile(pth,'avdiff.nii');
avdiff.vol = squeeze(sum(alldiff.vol))/size(controlpool,1);
avdiff.hdr.dt = [16 0];
spm_write_vol(avdiff.hdr,avdiff.vol)

% Yeah this is wrong version; T.vol = (squeeze((sqrt(pi)/2)*sum(abs(alldiff.vol))/size(controlpool,1))/sqrt(size(controlpool,1)));
T.vol = (sqrt(pi)/2)*squeeze(sum(abs(alldiff.vol)))/size(controlpool,1)*sqrt(1+1/size(controlpool,1));
% corrected according to the internet  version of the article, not the pdf
T.hdr = spm_vol(controlpool{1});
T.hdr.fname = fullfile(pth,'isas_denominator.nii');
spm_write_vol(T.hdr,T.vol)

Variance.hdr = spm_vol(controlpool{1});
Variance.vol = (sqrt(pi/2))*squeeze(sum(abs(alldiff.vol)))/size(controlpool,1);
Variance.hdr.fname = fullfile(pth,'variance.nii');
spm_write_vol(Variance.hdr,Variance.vol)


time.stats = toc;
end





















































