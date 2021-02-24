function dists = distsToSurface(data, subchainName, shape)
    % Compute the distance distribution of a ribosome subchain (Protein)
    % PARAMETERS:
    % data : a table storing coordinates as strings "[x, y, z]"
    % subchainName: the column name of the subchain
    % shape: the alpha shape object of the ribosome
    % RETURNS: 
    % A n by 1 array storing the surface distances of all atoms in the subchain

    P = subchainCoords(data, subchainName);
    [~, bfP] = boundaryFacets(shape);
    
    % Stores the distances
    n = size(P,1);
    dists = zeros(1,n);
    
    for j = 1:n
        dists(j) = min(vecnorm(bfP - P(j,:), 2, 2));
    end
end