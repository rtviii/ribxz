function msg = alpha(alpha_radius, dist_cutoff)
    % This file computes the surface ratios of a ribosome
    % PARAMETERS:
    % alpha_radius : radius of the alpha shape
    % dist_cutoff  : cutoff distance defined to be close to surface
    % Notes: data can be an input file
    data = readtable("../ribosomes/5AFI.csv");
    
    % Data cleaning -- formatting
    % This process will be different from different csv files
    data(:,1) = []; % delete the first column
    data.Properties.VariableNames = table2array(data(1,:));
    data(1,:) = [];

    % Data cleaning -- removing RNAs
    % Make a copy of the csv file and remove the RNAS
    dataPtns = data;
    dataPtns(:, "a") = [];
    dataPtns(:, "v") = [];
    dataPtns(:, "w") = [];
    dataPtns(:, "x") = [];
    dataPtns(:, "y") = [];
    dataPtns(:, "A") = [];
    dataPtns(:, "B") = [];

    % The rest are the same!
    % Compute alpha shape
    P = subunitCoords(data);
    shp = alphaShape(P, alpha_radius);

    % Compute surface ratios
    colNames = dataPtns.Properties.VariableNames;
    surfaceRatio = zeros(1,length(colNames));

    for j = 1:length(colNames)
        % This dists is the distance distribution per protein, maybe output this!
        dists = distsToSurface(data, colNames(j), shp);
        surfaceRatio(j) = sum(dists<dist_cutoff)/length(dists);
    end

    % Save surface ratios
    tb = table(transpose(colNames), transpose(surfaceRatio));
    tb.Properties.VariableNames = {'name', 'ratio'};
    writetable(tb, 'surface_ratio_5AFI_test.csv');
    msg = true;
end
