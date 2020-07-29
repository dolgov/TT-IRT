% Check that TT-IRT directories are in the path
function check_ttirt()
cd(fileparts(mfilename('fullpath')));
cd('..'); % TT-IRT root

try % Check if we have anything from constructors/ in the path
    amen_cross_s(2*ones(3,1), @(i)sum(i,2), 1e-1, 'verb', 0);
catch
    fprintf('Adding constructors/ to the path\n');
    cd('constructors')
    addpath(pwd)
    cd('..')
end

try % Check if we have anything from samplers/ in the path
    randref('N4', 2,2);
catch
    fprintf('Adding samplers/ to the path\n');
    cd('samplers')
    addpath(pwd)
    cd('..')
end

try % Check if we have anything from utils/ in the path
    tracemultm(rand(3,5), randi(5,3,1));
catch
    fprintf('Adding utils/ to the path\n');
    cd('utils')
    addpath(pwd)
    cd('..')    
end

end
