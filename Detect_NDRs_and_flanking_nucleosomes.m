function Detect_NDRs_and_flanking_nucleosomes(DyadsFilename, DatasetLabel)
% DETECT_NDRS_AND_FLANKING_NUCLEOSOMES  Detect the locations of the
% nucleosome depleted regions (NDRs) and the flanking nucleosomes (+1/-1)
% corresponding to each gene promoter.
%
% Inputs:
%   DyadsFilename - filename containing the cell array Dyads, each cell
%                   corresponding to a separate chromosome
%   DatasetLabel  - dataset label to include in the output file 
%
% Output:
%   A file named ['NDR_with_flanking_nucleosomes_',DatasetLabel,'.mat'] is
%   generated, which will contain the coordinates for the NDR center, and
%   +1/-1 nucleosomes
%
% Example:
% Detect_NDRs_and_flanking_nucleosomes('Dyads_AGH01-1_120_160.mat', 'AGH01-1')
%
% Note: A helping function is used (ConvertDyads2Occ), which converts a 
%   cell array containing the nucleosome dyads distributions into a cell
%   array containing the nucleosome occupancy profiles.
%
% See also CONVERTDYADS2OCC.
%
% Copyright (C) Razvan V. Chereji

if nargin < 2
	error('Not enough input arguments supplied!');
end
    
load(DyadsFilename, 'Dyads')
% Extend the dyads symmetrically to 101 bp, i.e. compute the nucleosome
% occupancy, enhancing the linkers by trimming the particle size to 101 bp
Occ = ConvertDyads2Occ(Dyads, 101);

%% Align TSS for all genes 
beforeTSS = 1000;
afterTSS = 1000;

% Chromosome lengths
chrLen = [230218,813184,316620,1531933,...
    576874,270161,1090940,562643,...
    439888,745751,666816,1078177,...
    924431,784333,1091291,948066];

% Load transcript annotations
load('Yeast_Transcripts.mat', 'ORF', 'Chr', 'Watson', 'TSS');

noGenes = numel(ORF);
AlignedOcc = nan(noGenes, 1 + beforeTSS + afterTSS);
for g = 1:noGenes
    if Watson(g) % if gene on the Watson (+) strand
        leftEdge = max([TSS(g) - beforeTSS, 1]);
        rightEdge = min([TSS(g) + afterTSS, chrLen(Chr(g))]);
        AlignedOcc(g, beforeTSS + 1 - (TSS(g) - leftEdge)...
            : beforeTSS + 1 + (rightEdge - TSS(g))) = ...
            Occ{Chr(g)}(leftEdge : rightEdge);
    else % if gene on the Crick (-) strand
        leftEdge = max([TSS(g) - afterTSS, 1]);
        rightEdge = min([TSS(g) + beforeTSS, chrLen(Chr(g))]);
        AlignedOcc(g, beforeTSS + 1 - (rightEdge - TSS(g))...
            : beforeTSS + 1 + (TSS(g) - leftEdge)) = ...
            fliplr(Occ{Chr(g)}(leftEdge : rightEdge)); % flip the coordinates for the genes on the Crick (-) strand
    end
end

RelativePlus1 = nan(noGenes, 1);
RelativeMinus1 = nan(noGenes, 1);

%%
NDRThreshold = 0.2 * mean(cat(2, [Occ{:}])); % Detect only nucleosomes with Occ > threshold (=0.2 x genome-wide average occupancy)
ProblematicGeneIndex = [];

for g = 1:noGenes
    y = AlignedOcc(g, :);
    
    % smooth the profile using a mobing average filter
    ys = smooth(y, 75);
    
    if max(ys) < NDRThreshold
        ProblematicGeneIndex = [ProblematicGeneIndex; g];
    end
    
    PeakThreshold = 0.6 * mean(ys);
    [pks, locs] = findpeaks(ys(1:1+beforeTSS+200), 'minPeakHeight', PeakThreshold, 'minPeakDistance', 100); % search fpr +1 up to TSS+200 (from 150)
    locs = locs - (1 + beforeTSS);
    if numel(locs) >= 2
        peakLocDiff = diff(locs);
        PossibleNDRIdx = find(peakLocDiff > 190);
        if ~isempty(PossibleNDRIdx)
            NDRcandidates = round(mean(locs([PossibleNDRIdx'; PossibleNDRIdx'+1])));
            
            Score = ys(1+beforeTSS+NDRcandidates)' ./ mean([pks(PossibleNDRIdx)'; pks(PossibleNDRIdx + 1)']);
            NDRcandidates(Score > 0.5) = [];
            Score(Score > 0.5) = [];
            
            if ~isempty(NDRcandidates)
                %                 [~, ClosestNDRIdx] = min(abs(NDRcandidates));
                [~, ClosestNDRIdx] = sort(abs(NDRcandidates), 'ascend');
                Score = Score(ClosestNDRIdx);
                %                 check whether we identified 2 NDRs close to TSS
                if sum(abs(NDRcandidates(ClosestNDRIdx)) < 200) > 1
                    ClosestNDRIdx = ClosestNDRIdx(abs(NDRcandidates(ClosestNDRIdx)) < 200);
                    Score = Score(abs(NDRcandidates(ClosestNDRIdx)) < 200);
                    
                    [~, k] = min(Score);
                    ClosestNDRIdx = ClosestNDRIdx(k);
                else
                    ClosestNDRIdx = ClosestNDRIdx(1);
                end
                
                NDRapprox = NDRcandidates(ClosestNDRIdx);
                
                k = find(locs < NDRapprox, 1, 'last');
                RelativeMinus1(g) = locs(k);
                k = find(locs > NDRapprox, 1, 'first');
                RelativePlus1(g) = locs(k);
              
            end
        end
    end
    
    if isnan(RelativeMinus1(g))
        ProblematicGeneIndex = [ProblematicGeneIndex; g];
    end

end

%%
% Construct the absolute Plus1, Minus1, coordinates
Plus1 = nan(size(RelativePlus1));
Plus1(Watson) = TSS(Watson) + RelativePlus1(Watson);
Plus1(~Watson) = TSS(~Watson) - RelativePlus1(~Watson);

Minus1 = nan(size(RelativeMinus1));
Minus1(Watson) = TSS(Watson) + RelativeMinus1(Watson);
Minus1(~Watson) = TSS(~Watson) - RelativeMinus1(~Watson);

NDRcenter = round(mean([Plus1, Minus1], 2));
NDRwidth = abs(Plus1 - Minus1) - 147;

ProblematicGenes = ORF(ProblematicGeneIndex);

clearvars -except ORF Watson Chr TSS RelativePlus1 RelativeMinus1 Plus1 Minus1 NDRcenter NDRwidth ProblematicGenes ProblematicGeneIndex DatasetLabel
save(['NDR_with_flanking_nucleosomes_',DatasetLabel,'.mat'])