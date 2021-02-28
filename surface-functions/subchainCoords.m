function P = subchainCoords(data, subchainName)
    % Extract the coordinates of a subchain of a ribosome
    % PARAMETERS:
    % data : a table storing coordinates as strings "[x, y, z]"
    % subchainName: the column name of the subchain
    % RETURNS:
    % An n by 3 matrix with the 3D coordinates of all atoms in the subchain
    
    subchain = table2cell(rmmissing(data(:,subchainName)));
    n = length(subchain);
    X = zeros(n, 1);
    Y = zeros(n, 1);
    Z = zeros(n, 1);
    
    for i=1:n
        coords = rmmissing(strsplit(subchain{i}(2:end-1)));
        X(i) = str2double(coords{1});
        Y(i) = str2double(coords{2});
        Z(i) = str2double(coords{3});
    end
    
    P = [X, Y, Z];
end