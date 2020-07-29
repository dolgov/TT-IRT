% Installation script for TT-IRT toolbox

cd(fileparts(mfilename('fullpath'))); % Change to TT-IRT root
cd('utils')
check_ttirt; % add our own paths

FilesToPatch = {'samplers/tt_irt_lin.m', 'samplers/tt_irt_sqr.m', 'samplers/tt_irt_fourier.m'};
fprintf('Checking if tracemult MEX is available on your architecture...');
try
    tracemult(rand(3,5), randi(5,3,1));
    fprintf('yes\n');
    for fname = FilesToPatch
        fprintf('Patching %s to MEX tracemult\n', fname{1});
        txt = fileread(fname{1});
        if (contains(txt, 'tracemultm'))
            txt = strrep(txt, 'tracemultm', 'tracemult');
            f = fopen(fname{1}, 'w');
            fwrite(f, txt);
            fclose(f);
        end
    end
    
catch ME
    fprintf('no\n');
    req = input('Compile tracemult MEX? ([true]/false) ');
    if (isempty(req))
        req = 1;
    end
    fprintf('\tGot %g\n', req);
    if (req)
        ME2 = [];
        try
            cd('utils');
            mex -v -largeArrayDims -lmwblas -lmwlapack -O tracemult.c
            fprintf('\n');
            cd('..');
        catch ME2
        end
        if (isempty(ME2))
            % Compilation succeeded
            for fname = FilesToPatch
                fprintf('Patching %s to MEX tracemult\n', fname{1});
                txt = fileread(fname{1});
                if (contains(txt, 'tracemultm'))
                    txt = strrep(txt, 'tracemultm', 'tracemult');
                    f = fopen(fname{1}, 'w');
                    fwrite(f, txt);
                    fclose(f);
                end
            end
        else
            % Compilation failed
            warning(ME2.message)
            for fname = FilesToPatch
                fprintf('Patching %s to native tracemultm\n', fname{1});
                txt = fileread(fname{1});
                if (~contains(txt, 'tracemultm'))
                    txt = strrep(txt, 'tracemult', 'tracemultm');
                    f = fopen(fname{1}, 'w');
                    fwrite(f, txt);
                    fclose(f);
                end
            end
        end
    else % compilation is not wanted
        for fname = FilesToPatch
            fprintf('Patching %s to native tracemultm\n', fname{1});
            txt = fileread(fname{1});
            if (~contains(txt, 'tracemultm'))
                txt = strrep(txt, 'tracemult', 'tracemultm');
                f = fopen(fname{1}, 'w');
                fwrite(f, txt);
                fclose(f);
            end
        end
    end
end

req = input('Use parfor for different runs? (true/[false]) ');
if (isempty(req))
    req = 0;
end
fprintf('\tGot %g\n', req);
FilesToPatch = {'examples/shock_absorber/test_shock_absorber_tt.m', ...
                'examples/shock_absorber/test_shock_absorber_dram.m', ...
                'examples/diffusion/test_diffusion_tt.m', ...
                'examples/diffusion/test_diffusion_dram.m', ...
                'examples/diffusion/test_diffusion_qmcrat.m', ...
                'examples/predator_prey/test_predator_prey_dirt.m', ...
                'examples/predator_prey/test_predator_prey_dram.m', ...
                'examples/predator_prey/test_predator_prey_svn.m'};
if (req)    
    for fname=FilesToPatch
        fprintf('Patching %s to parfor\n', fname{1});
        txt = fileread(fname{1});
        if (~contains(txt, 'parfor irun'))
            txt = strrep(txt, 'for irun', 'parfor irun');
            f = fopen(fname{1}, 'w');
            fwrite(f, txt);
            fclose(f);
        end
    end
else
    for fname=FilesToPatch
        fprintf('Patching %s to sequential for\n', fname{1});
        txt = fileread(fname{1});
        if (contains(txt, 'parfor irun'))
            txt = strrep(txt, 'parfor irun', 'for irun');
            f = fopen(fname{1}, 'w');
            fwrite(f, txt);
            fclose(f);
        end
    end
end



fprintf('Checking if als-cross MEXes are available on your architecture...');
try
    UAU = rand(3,3,5);
    % Make the left environment well-conditined for testing purposes
    UAU = UAU + permute(UAU,[2,1,3]) + repmat(eye(3),1,1,5);
    project_blockdiag_mex(UAU, rand(5,7,6), rand(3,7,4));
    solve_blockdiag_mex(UAU, rand(5,7,4), rand(3,7,4));
    fprintf('yes\n');
    
catch ME
    fprintf('no\n');
    req = input('Compile MEXes for als-cross? ([true]/false) ');
    if (isempty(req))
        req = 1;
    end
    fprintf('\tGot %g\n', req);
    if (req)
        cd('utils');
        mex -v -largeArrayDims -lmwblas -lmwlapack -O solve_blockdiag_mex.c
        fprintf('\n');
        mex -v -largeArrayDims -lmwblas -lmwlapack -O project_blockdiag_mex.c
        fprintf('\n');
        cd('..');
    end
end

% % =======================================================================
% % MEXed tt_irt_sqr and fourier would be much more difficult to maintain,
% % hence the MEXed irt is not supported anymore. In principle, all
% % operations except tracemult are well vectorised, so once tracemult is
% % MEXed, the code should be fast enough.
% % =======================================================================
% 
% req = input('Compile IRT MEX? (true/false) ');
% if (isempty(req))
%     req = 1;
% end
% fprintf('\tGot %g\n', req);
% fname = 'tt_irt_debias.m';
% if (req)
%     ME = [];
%     try
%         mex -v -largeArrayDims -lmwblas -lmwlapack -O tt_irt_mex.c tt_irt1_int64.c
%         fprintf('\n');
%     catch ME
%     end
%     if (isempty(ME))
%         % Compilation succeeded
%         fprintf('Patching %s to tt_irt_mex\n', fname);
%         txt = fileread(fname);
%         if (contains(txt, 'tt_irt(xsf, f, Z)'))
%             txt = strrep(txt, 'tt_irt(xsf, f, Z)', 'tt_irt_mex(f.n, cell2mat(xsf), f.r, f.core, Z)');
%             f = fopen(fname, 'w');
%             fwrite(f, txt);
%             fclose(f);
%         end
%     else
%         % Compilation failed
%         fprintf('Patching %s to native tt_irt\n', fname);
%         txt = fileread(fname);
%         if (contains(txt, 'tt_irt_mex(f.n, cell2mat(xsf), f.r, f.core, Z)'))
%             txt = strrep(txt, 'tt_irt_mex(f.n, cell2mat(xsf), f.r, f.core, Z)', 'tt_irt(xsf, f, Z)');
%             f = fopen(fname, 'w');
%             fwrite(f, txt);
%             fclose(f);
%         end
%         rethrow(ME)
%     end
% else
%     fprintf('Patching %s to native tt_irt\n', fname);
%     txt = fileread(fname);
%     if (contains(txt, 'tt_irt_mex(f.n, cell2mat(xsf), f.r, f.core, Z)'))
%         txt = strrep(txt, 'tt_irt_mex(f.n, cell2mat(xsf), f.r, f.core, Z)', 'tt_irt(xsf, f, Z)');
%         f = fopen(fname, 'w');
%         fwrite(f, txt);
%         fclose(f);
%     end
% end
% 
