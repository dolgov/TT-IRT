function check_svn()
% Check and download SVN repository

cd(fileparts(mfilename('fullpath')));
cd('..'); % examples
cd('..'); % TT-IRT root
cd('samplers')

if (exist('SVN_H', 'file')==0)
    if (exist('Stein-variational-samplers-master', 'dir')==0)
        if (exist('Stein-variational-samplers-master.zip', 'file')==0)
            try
                fprintf('Stein-variational-samplers is not found, downloading...\n');
                opts = weboptions; opts.CertificateFilename=(''); % SUXX
                websave('Stein-variational-samplers-master.zip', 'https://github.com/gianlucadetommaso/Stein-variational-samplers/archive/master.zip', opts);
            catch ME
                error('%s. Automatic download failed. Please download Stein-variational-samplers-master.zip from https://github.com/gianlucadetommaso/Stein-variational-samplers', ME.message);
            end
            fprintf('Success!\n');
        end
        try
            unzip('Stein-variational-samplers-master.zip');
        catch ME
            error('%s. Automatic unzipping failed. Please extract Stein-variational-samplers-master.zip here', ME.message);
        end
    end
    cd('Stein-variational-samplers-master');
    cd('MATLAB');
    cd('samplers');
    cd('stein');
    addpath(pwd)
    cd('..'); % samplers
    cd('..'); % MATLAB
    cd('..'); % SVSM
    cd('..'); % samplers
    cd('..'); % TT-IRT root
end
end