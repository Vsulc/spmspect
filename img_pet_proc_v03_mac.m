function img_pet_proc_v03_mac

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
%
% SPECT has to be manualy realigned for coregistration to match, roughly AC

addpath('/Users/vlasta/Neuro/spm12')
spm_jobman('initcfg');
pet.smooth = 4; %technically c1 smooth

[pet.name, pet.path, ~] = uigetfile('*PET*.nii', 'Pick PET file');
pet.fullfile = fullfile(pet.path,pet.name);

[mriT1.name, mriT1.path, mriT1.exists] = uigetfile('*MRI*.nii', 'Pick MRI or rMRI file');
mriT1.fullfile = fullfile(mriT1.path,mriT1.name);

[c1mri.name, c1mri.pathname, c1mri.exists] = uigetfile('*c1*.nii', 'Pick c1 file');
c1mri.fullfile = fullfile(c1mri.pathname,c1mri.name);

if c1mri.exists==1
            c2mri.name= strrep(c1mri.name, 'c1', 'c2');
            c3mri.name= strrep(c1mri.name, 'c1', 'c3');
       if   exist(c2mri.name,'file')  && exist(c3mri.name,'file')
           skipseg=1;
       else
           skipseg=0;
       end
end


 if mriT1.exists==0   
        normalize2std{1}.spm.tools.oldnorm.estwrite.subj.source{1} = pet.fullfile;
        normalize2std{1}.spm.tools.oldnorm.estwrite.subj.resample = {pet.fullfile;};  
        normalize2std{1}.spm.tools.oldnorm.estwrite.eoptions.template = {fullfile(spm('Dir'),'templates','PET.nii')};
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
            clear jobs;

 elseif mriT1.exists==1
        disp('Coregistering PET to MRI...')
        matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {mriT1.fullfile};
        matlabbatch{1}.spm.spatial.coreg.estwrite.source = {pet.fullfile};
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
        
            if skipseg==0   
                disp('Segmenting MRI...')
                matlabbatch{1}.spm.spatial.preproc.channel.vols = {mriT1.fullfile};
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
                matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
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
                
                c1mri.name = ['c1' mriT1.name];
                c1mri.fullfile = fullfile(mriT1.path,c1mri.name);

            end
        
        disp('Smoothing grey matter mask from MRI...')
        smoothpet{1}.spm.spatial.smooth.data = {c1mri.fullfile};
        smoothpet{1}.spm.spatial.smooth.fwhm = [pet.smooth pet.smooth pet.smooth];
        smoothpet{1}.spm.spatial.smooth.dtype = 0;
        smoothpet{1}.spm.spatial.smooth.im = 0;
        smoothpet{1}.spm.spatial.smooth.prefix = 's';
        spm_jobman('run',smoothpet);
        
        c1mri.sname = ['s'  c1mri.name];
        c1mri.sfullfile =  fullfile(mriT1.path,c1mri.sname); %% fix 01
        c1mri.hdr = spm_vol(c1mri.sfullfile);
        c1mri.vol = spm_read_vols(c1mri.hdr);
        
        pet.rname = ['r' pet.name];
        disp(pet.rname)
        pet.rfullfile = fullfile(pet.path,pet.rname);

        pet.hdr = spm_vol(pet.rfullfile);
        
        pet.vol = spm_read_vols(pet.hdr);
        pet.rvol = pet.vol.*(c1mri.vol>.2);
        pet.rvol(pet.rvol  == 0) = NaN;
        
%         disp('Saving data...')
%         pet.rhdr = pet.hdr;
%         pet.rhdr.fname = fullfile(pet.path,['sr' pet.name]);
%         pet.rhdr.private.dat.fname = pet.hdr.fname;
%         pet.rhdr.private.dat.dtype = 'FLOAT32-LE';
%         pet.rhdr.dt = [16 0];
% 
%         spm_write_vol(pet.rhdr,pet.rvol);
        
        disp('Saving data mean 50 data ...')
        pet.rThr = nanmean(pet.rvol(:));
        pet.rcvol = pet.rvol ./ pet.rThr .* 50;
        pet.rchdr = pet.hdr;
        pet.rchdr.fname = fullfile(pet.path,['src' pet.name]);
        pet.rchdr.private.dat.fname = pet.hdr.fname;
        pet.rchdr.private.dat.dtype = 'FLOAT32-LE';
        pet.rchdr.dt = [16 0];


        spm_write_vol(pet.rchdr,pet.rcvol);
        
        
        %% global normalization
        
        if skipseg==0 
            matlabbatch{1}.spm.util.imcalc.input = {
                                       fullfile(mriT1.path,['c1' mriT1.name])
                                       fullfile(mriT1.path,['c2' mriT1.name])
                                       fullfile(mriT1.path,['c3' mriT1.name])
                                                };
        elseif skipseg==1
            matlabbatch{1}.spm.util.imcalc.input = {
                                        fullfile(c1mri.pathname,c1mri.name)
                                        fullfile(c1mri.pathname,c2mri.name)
                                        fullfile(c1mri.pathname,c3mri.name)
                                        };
            
        end
        matlabbatch{1}.spm.util.imcalc.output = 'brainmask_0';
        matlabbatch{1}.spm.util.imcalc.outdir = {''}; % fix to ictal.path
        matlabbatch{1}.spm.util.imcalc.expression = 'i1+i2+i3';
        %matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

        matlabbatch{2}.spm.spatial.smooth.data(1) = {fullfile(mriT1.path,'brainmask_0.nii')};
        matlabbatch{2}.spm.spatial.smooth.fwhm = [8 8 8];
        matlabbatch{2}.spm.spatial.smooth.dtype = 0;
        matlabbatch{2}.spm.spatial.smooth.im = 0;
        matlabbatch{2}.spm.spatial.smooth.prefix = 's';

        matlabbatch{3}.spm.util.imcalc.input(1) = {fullfile(mriT1.path,'sbrainmask_0.nii')};
        matlabbatch{3}.spm.util.imcalc.output = 'brainmask';
        matlabbatch{3}.spm.util.imcalc.outdir = {''};  % fix to ictal.path
        matlabbatch{3}.spm.util.imcalc.expression = 'i1>0.5';
        matlabbatch{3}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{3}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{3}.spm.util.imcalc.options.mask = 0;
        matlabbatch{3}.spm.util.imcalc.options.interp = 1;
        matlabbatch{3}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run',matlabbatch);
        clear matlabbatch
        
        
        matlabbatch{1}.spm.util.imcalc.input = {
                                       fullfile(pet.path,['r' pet.name])
                                       fullfile(pet.path,'brainmask.nii,1')
                                        };
        matlabbatch{1}.spm.util.imcalc.output = fullfile(pet.path,['mr' pet.name]);
        matlabbatch{1}.spm.util.imcalc.outdir = {''};
        matlabbatch{1}.spm.util.imcalc.expression = 'i1.*(i2>0.5)';
        %matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

        spm_jobman('run',matlabbatch);
        clear matlabbatch

        %ictal
        pet.hdr = spm_vol(fullfile(pet.path,['mr' pet.name]));
        pet.vol = spm_read_vols(pet.hdr);
        pet.nanvol = pet.vol;
        pet.nanvol(~pet.nanvol)=NaN;
        rThr = nanmean(pet.nanvol(:));
        pet.rcvol = pet.vol ./ rThr .* 50;

        pet.rchdr = pet.hdr;
        pet.rchdr.fname = fullfile(pet.path,['cmr' pet.name]);
        pet.rchdr.private.dat.fname = pet.rchdr.fname;
        pet.rchdr.private.dat.dtype = 'FLOAT32-LE'; %'INT16-LE', 'FLOAT32-LE'
        pet.rchdr.dt = [16 0]; % 4=16-bit integer; 16=32-bit real datatype

        spm_write_vol(pet.rchdr,pet.rcvol);

 end




end

