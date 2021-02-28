function P = subunitCoords(data)
    % Extract the 3D coordinates in the cells of the table
    % PARAMETERS:
    % data : a table storing coordinates as strings "[x, y, z]"
    % RETURNS: 
    % An n by 3 matrix with the 3D coordinates of all atoms in the table
    
    Xo = [];
    Yo = [];
    Zo = [];
    
    for j = 1:width(data)
        subchain = rmmissing(table2array(data(:,j)));
        nj = length(subchain);
        X = zeros(1, nj);
        Y = zeros(1, nj);
        Z = zeros(1, nj);
        
        for i = 1:nj
            coords = rmmissing(strsplit(subchain{i}(2:end-1)));
            X(i) = str2double(coords{1});
            Y(i) = str2double(coords{2});
            Z(i) = str2double(coords{3});
        end
    
        Xo = [Xo, X];
        Yo = [Yo, Y];
        Zo = [Zo, Z];
    end

    P = [transpose(Xo), transpose(Yo), transpose(Zo)];
end

