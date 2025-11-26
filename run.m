function run(tabelaPath)
% =========================================================================
% Exemple: run('dataset.csv')
% =========================================================================

    data = readtable(tabelaPath);

    nCases = height(data);
    fprintf('Run %d cases...\n', nCases);

    time = tic;

    for i = 1:nCases
        fprintf('\n===== Run examples %d de %d =====\n', i, nCases);
        fprintf('Parameters: nelx=%d, nely=%d, volfrac=%.2f, penal=%.1f, rmin=%.2f, newF=%.1f, problem=%d, geometry=%d\n', ...
            data.nelx(i), data.nely(i), data.volfrac(i), data.penal(i), data.rmin(i), ...
            data.newF(i), data.problem(i), data.geometry(i));

        % Call SIMP
        try
            SIMP(data.nelx(i), ...
                 data.nely(i), ...
                 data.volfrac(i), ...
                 data.penal(i), ...
                 data.rmin(i), ...
                 data.newF(i), ...
                 data.problem(i), ...
                 data.geometry(i));
        catch ME
            warning('Error when running examples %d: %s\n', i, ME.message);
        end

    timeTotal = toc(time);
    fprintf('\nAll cases were executed in %.2f seconds (%.2f minutes).\n', ...
            timeTotal, timeTotal/60);
end