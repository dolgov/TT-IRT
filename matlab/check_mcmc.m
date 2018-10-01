function check_mcmc()
if (exist('dramrun', 'file')==0)
    if (exist('dramcode', 'dir')==0)
        if (exist('dramcode.zip', 'file')==0)
            try
                fprintf('DRAM is not found. Downloading...\n');
                urlwrite('http://helios.fmi.fi/~lainema/dram/dramcode.zip', 'dramcode.zip');
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
    cd('..')
end

if (exist('UWerr', 'file')==0)
    try
        fprintf('UWerr is not found. Downloading...\n');
        urlwrite('https://www.physik.hu-berlin.de/de/com/UWerr.m', 'UWerr.m');
    catch ME
        error('%s. Please download UWerr from https://www.physik.hu-berlin.de/de/com/ALPHAsoft', ME.message);
    end
    fprintf('Success!\n');
end
end

