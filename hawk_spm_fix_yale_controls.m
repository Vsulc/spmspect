function hawk_spm_fix_yale_controls

% addpath('/Users/vlasta/Neuro/spm12')
% spm_jobman('initcfg');

controls_dir = '/Users/vlasta/Neuro/controls_preprocessing/Raw_Healthy_Normals';

% fixmatrix =    [ 	1     0     0     0
%                  	0     1     0     0
%                     0     0    -1     0
%                     0     0     0     1];

ictalpool =   {
                    'HN001_D1.img,1'  
                    'HN002_D1.img,1'  
                    'HN003_D1.img,1'  
                    'HN004_D1.img,1'  
                    'HN005_D1.img,1'  
                    'HN006_D1.img,1'  
                    'HN007_D1.img,1'  
                    'HN008_D1.img,1'  
                    'HN009_D1.img,1'  
                    'HN010_D1.img,1'  
                    'HN011_D1.img,1'  
                    'HN012_D1.img,1'  
                    'HN013_D1.img,1'  
                    'HN014_D1.img,1'  
                    'HN001_D2.img,1'  
                    'HN002_D2.img,1'  
                    'HN003_D2.img,1'  
                    'HN004_D2.img,1'  
                    'HN005_D2.img,1'  
                    'HN006_D2.img,1'  
                    'HN007_D2.img,1'  
                    'HN008_D2.img,1'  
                    'HN009_D2.img,1'  
                    'HN010_D2.img,1'  
                    'HN011_D2.img,1'  
                    'HN012_D2.img,1'  
                    'HN013_D2.img,1'  
                    'HN014_D2.img,1'  
                    };
                
    for f=1:size(ictalpool)
        spect.fullfile = fullfile(controls_dir,char(ictalpool(f)));                
        [spect.path,spect.name,spect.ext] = fileparts(spect.fullfile);
        ictal.hdr = spm_vol(spect.fullfile);
        ictal.vol = spm_read_vols(ictal.hdr);
        ictal.fvol = flipdim(ictal.vol,3);
        ictal.fhdr = ictal.hdr;
        ictal.fhdr.fname = fullfile(controls_dir,['f' spect.name '.nii']);
        ictal.fhdr.private.dat.fname = ictal.fhdr.fname;
        ictal.fhdr.private.dat.dtype = 'FLOAT32-LE'; %'INT16-LE', 'FLOAT32-LE'
        ictal.fhdr.dt = [16 0]; %
        spm_write_vol(ictal.fhdr,ictal.fvol);

    end
end
