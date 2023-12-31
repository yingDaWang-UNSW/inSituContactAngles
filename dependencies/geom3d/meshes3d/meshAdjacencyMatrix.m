function adj = meshAdjacencyMatrix(faces, varargin)
%MESHADJACENCYMATRIX Compute adjacency matrix of a mesh from set of faces
%
%   ADJMAT = meshAdjacencyMatrix(FACES)
%   Returns a sparse NV-by-NV matrix (NV being the maximum vertex index)
%   containing vertex adjacency of the mesh represented by FACES.
%   FACES is either a NF-by-3, a NF-by-4 index array, or a Nf-by-1 cell
%   array.
%
%   Example
%     [v f] = createCube;
%     adj = meshAdjacencyMatrix(f);
%
%   See also
%     meshes3d, triangulateFaces, smoothMesh

% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2013-04-30,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2013 INRA - Cepia Software Platform.

% Ensures faces is a N-by-3 or N-by-4 array
if iscell(faces) || (isnumeric(faces) && size(faces, 2) > 4)
    faces = triangulateFaces(faces);
end

% forces faces to b efloating point array, for sparse function
if ~isfloat(faces)
    faces = double(faces);
end
    
% populate a sparse matrix
if size(faces, 2) == 3
    adj = sparse(...
        [faces(:,1); faces(:,1); faces(:,2); faces(:,2); faces(:,3); faces(:,3)], ...
        [faces(:,3); faces(:,2); faces(:,1); faces(:,3); faces(:,2); faces(:,1)], ...
        1.0);
elseif size(faces, 2) == 4
    adj = sparse(...
        [faces(:,1); faces(:,1); faces(:,2); faces(:,2); faces(:,3); faces(:,3); faces(:,4); faces(:,4)], ...
        [faces(:,4); faces(:,2); faces(:,1); faces(:,3); faces(:,2); faces(:,4); faces(:,3); faces(:,1)], ...
        1.0);
    
elseif size(faces, 2) == 2
    adj = sparse(...
        [faces(:,1); faces(:,2)], ...
        [faces(:,2); faces(:,1)], ...
        1.0);
end
   
% remove double adjacencies
adj = min(adj, 1);
