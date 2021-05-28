function check_mcmc(need_dram)
cd(fileparts(mfilename('fullpath')));
cd('..'); % TT-IRT root
cd('samplers')

if (nargin<1)
    need_dram = true;
end

if ((exist('dramrun', 'file')==0)&&(need_dram))
    if (exist('dramcode', 'dir')==0)
        if (exist('dramcode.zip', 'file')==0)
            try
                fprintf('DRAM is not found. Downloading...\n');
                opts = weboptions; opts.CertificateFilename=(''); % SUXX
                websave('dramcode.zip', 'http://helios.fmi.fi/~lainema/dram/dramcode.zip', opts);
            catch ME
                error('%s. Please download DRAM from http://helios.fmi.fi/~lainema/dram/', ME.message);
            end
            fprintf('Success!\n');
        end
        try
            mkdir('dramcode')
            movefile('dramcode.zip', 'dramcode')
            cd('dramcode')
            unzip('dramcode.zip');
            cd('utils')
            movefile('*', '..')
            cd('..') % in dramcode
            cd('..') % in matlab
        catch ME
            error('%s. Automatic unzipping failed. Please put ALL files from dramcode.zip into dramcode directory', ME.message);
        end        
    end
    cd('dramcode')
    addpath(pwd)
    cd('..') % samplers
end

if (exist('UWerr', 'file')==0)
    try
        fprintf('UWerr is not found. Downloading...\n');
        opts = weboptions; opts.CertificateFilename=(''); % SUXX
        websave('UWerr.m', 'https://www.physik.hu-berlin.de/de/com/UWerr.m', opts);
    catch ME
        error('%s. Please download UWerr from https://www.physik.hu-berlin.de/de/com/ALPHAsoft', ME.message);
    end
    fprintf('Success!\n');
end

cd('..'); % TT-IRT root
end

