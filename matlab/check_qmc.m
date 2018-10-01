% Check/download lattice vector
function check_qmc()
if (exist('lattice-39102-1024-1048576.3600.txt', 'file')==0)
    try
        fprintf('lattice-39102-1024-1048576.3600.txt is not found, downloading...\n');
        urlwrite('http://web.maths.unsw.edu.au/~fkuo/lattice/lattice-39102-1024-1048576.3600', 'lattice-39102-1024-1048576.3600.txt');
    catch ME
        error('%s. Automatic download failed. Please download lattice-39102-1024-1048576.3600.txt into this directory', ME.message);
    end
end
end
