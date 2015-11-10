function Occ = ConvertDyads2Occ(Dyads, a)
% CONVERTDYADS2OCC  Converts a cell array containing the nucleosome dyads
% distributions into a cell array containing the nucleosome occupancy
% profiles.
%
% Inputs:
%   Dyads - cell array, each cell corresponding the nucleosome dyad counts
%           for a separate chromosome
%   a     - footprint of the particle (typical size for nucleosomes = 147)
%
% Output:
%   Occ   - cell array, each cell corresponding the nucleosome occupancy
%           for a separate chromosome
%
% Example:
% load('Dyads_AGH01-1_120_160.mat', 'Dyads') % load Dyads
% Occ = ConvertDyads2Occ(Dyads, 147)         % compute corresponding Occ
%
% See also DETECT_NDRS_AND_FLANKING_NUCLEOSOMES.
%
% Copyright (C) Razvan V. Chereji

halfParticle = floor(a/2);
noChr = numel(Dyads);
Occ = cell(size(Dyads));
for c = 1:noChr
    Occ{c} = filter(ones(1,a), 1, [Dyads{c}(halfParticle + 1 : end), zeros(1, halfParticle)]);
end