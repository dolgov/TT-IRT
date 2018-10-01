% Installation script for TT-IRT toolbox

req = input('Compile tracemult MEX? (true/false) ');
if (isempty(req))
    req = 1;
end
fprintf('\tGot %g\n', req);
fname = 'tt_irt.m';
if (req)
    ME = [];
    try
        mex -v -largeArrayDims -lmwblas -lmwlapack -O tracemult.c
        fprintf('\n');
    catch ME        
    end
    if (isempty(ME))
        % Compilation succeeded
        fprintf('Patching %s to MEX tracemult\n', fname);
        txt = fileread(fname);
        if (contains(txt, 'tracemultm'))
            txt = strrep(txt, 'tracemultm', 'tracemult');
            f = fopen(fname, 'w');
            fwrite(f, txt);
            fclose(f);
        end
    else
        % Compilation failed
        fprintf('Patching %s to native tracemultm\n', fname);
        txt = fileread(fname);
        if (~contains(txt, 'tracemultm'))
            txt = strrep(txt, 'tracemult', 'tracemultm');
            f = fopen(fname, 'w');
            fwrite(f, txt);
            fclose(f);
        end
        rethrow(ME)
    end
else
    fprintf('Patching %s to native tracemultm\n', fname);
    txt = fileread(fname);
    if (~contains(txt, 'tracemultm'))
        txt = strrep(txt, 'tracemult', 'tracemultm');
        f = fopen(fname, 'w');
        fwrite(f, txt);
        fclose(f);
    end
end




req = input('Compile IRT MEX? (true/false) ');
if (isempty(req))
    req = 1;
end
fprintf('\tGot %g\n', req);
fname = 'tt_irt_debias.m';
if (req)
    ME = [];
    try
        mex -v -largeArrayDims -lmwblas -lmwlapack -O tt_irt_mex.c tt_irt1_int64.c
        fprintf('\n');
    catch ME
    end
    if (isempty(ME))
        % Compilation succeeded
        fprintf('Patching %s to tt_irt_mex\n', fname);
        txt = fileread(fname);
        if (contains(txt, 'tt_irt(xsf, f, Z)'))
            txt = strrep(txt, 'tt_irt(xsf, f, Z)', 'tt_irt_mex(f.n, cell2mat(xsf), f.r, f.core, Z)');
            f = fopen(fname, 'w');
            fwrite(f, txt);
            fclose(f);
        end
    else
        % Compilation failed
        fprintf('Patching %s to native tt_irt\n', fname);
        txt = fileread(fname);
        if (contains(txt, 'tt_irt_mex(f.n, cell2mat(xsf), f.r, f.core, Z)'))
            txt = strrep(txt, 'tt_irt_mex(f.n, cell2mat(xsf), f.r, f.core, Z)', 'tt_irt(xsf, f, Z)');
            f = fopen(fname, 'w');
            fwrite(f, txt);
            fclose(f);
        end
        rethrow(ME)
    end
else
    fprintf('Patching %s to native tt_irt\n', fname);
    txt = fileread(fname);
    if (contains(txt, 'tt_irt_mex(f.n, cell2mat(xsf), f.r, f.core, Z)'))
        txt = strrep(txt, 'tt_irt_mex(f.n, cell2mat(xsf), f.r, f.core, Z)', 'tt_irt(xsf, f, Z)');
        f = fopen(fname, 'w');
        fwrite(f, txt);
        fclose(f);
    end
end



req = input('Use parfor for different runs? (true/false) ');
if (isempty(req))
    req = 1;
end
fprintf('\tGot %g\n', req);
if (req)    
    for fname={'test_shock_absorber_tt.m', 'test_shock_absorber_dram.m', 'test_diffusion_tt.m', 'test_diffusion_dram.m'}
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
    for fname={'test_shock_absorber_tt.m', 'test_shock_absorber_dram.m', 'test_diffusion_tt.m', 'test_diffusion_dram.m'}
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



req = input('Compile MEXes for als-cross? (true/false) ');
if (isempty(req))
    req = 1;
end
fprintf('\tGot %g\n', req);
if (req)
    mex -v -largeArrayDims -lmwblas -lmwlapack -O solve_blockdiag_mex.c
    fprintf('\n');
    mex -v -largeArrayDims -lmwblas -lmwlapack -O project_blockdiag_mex.c
    fprintf('\n');    
end

