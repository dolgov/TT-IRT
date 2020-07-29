% Check/download lattice vector
function check_qmc()
cd(fileparts(mfilename('fullpath')));
cd('..'); % TT-IRT root
cd('samplers')

if (exist('lattice-39102-1024-1048576.3600.txt', 'file')==0)
    try
        fprintf('lattice-39102-1024-1048576.3600.txt is not found, downloading...\n');
        opts = weboptions; opts.CertificateFilename=(''); % SUXX
        websave('lattice-39102-1024-1048576.3600.txt', 'http://web.maths.unsw.edu.au/~fkuo/lattice/lattice-39102-1024-1048576.3600', opts);
    catch ME
        error('%s. Automatic download failed. Please download lattice-39102-1024-1048576.3600.txt into this directory', ME.message);
    end
    fprintf('Success!\n');
end

cd('..'); % TT-IRT root
end
