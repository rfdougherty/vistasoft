function [grads, maxB,delta,Delta] = dtiGradsBuildCharmed(bvals, nDirs, nReps, maxG, Delta, delta, interact)
%
% [grads, maxB] = dSimGetCharmedGrads(bvals, nDirs, [nReps], [maxG = 50.0], [Delta=optimal value to minimize TE], [delta=Delta])
%
% bvals is a list of b-values that you want (i.e., shells) in units of ms/um^2. Don't forget the b-0!
% nDirs is a list (same size as bvals) for the # of directions for each b-value. (For b=0, this is the # of b=0 images.)
%
% E.g., for 2-shells at b=1500,3000 (1.5,3.0 in ms/um^2 units) with 75 directions each and 10 b=0 images:
%
% grads = dtiGradsBuildCharmed([0,1.5,3.0], [10,75,75]);
%
% for insertion into a GE tensor.dat file:
% fp = fopen('/tmp/tensor.dat','w'); fprintf(fp, '%d\n', size(grads,2)); for i=1:size(grads,2), fprintf(fp, '%0.6f %06f %06f\n', grads(:,i)); end; fclose(fp);
%
% More examples:
%
% A single shell with b=2500, 96 directions, and 10 b=0:
%
% grads = dtiGradsBuildCharmed([0,2.5], [10,96]);
%
% NODI-optimized 2-shell and 3-shell:
%
% grads = dtiGradsBuildCharmed([0,0.7,2.5], [9,30,60]);
%
% grads = dtiGradsBuildCharmed([0,0.9,1.8,2.7], [11,30,45,65]);
%
% A set of shells for 'CHARMED':
%
% bvals = [0 714 1428 2285 3214 4286 5357 6429 7500 8571 10000]./1000;
% nDirs = [1   6    7   11   13   15   17   19   21   29    31];
% nReps = [16  2    2    2    2    2    2    2    2    2     2];
% grads = dtiGradsBuildCharmed(bvals, nDirs, nReps);
%
% HISTORY:
% 2009.05.13 RFD & AM wrote it.% bvals = [0  900 1800 2700]./1000;

if(~exist('nReps','var') || isempty(nReps))
    nReps = 1;
end
if(numel(nReps)<numel(bvals))
    nReps = repmat(nReps(1), 1, numel(bvals));
end
if(~exist('maxG','var') || isempty(maxG))
    maxG = 50.0; % Max gradient in mT/m
end
if(~exist('interact','var') || isempty(interact))
    interact = 0;
end

g = 42576.0; % gyromagnetic ratio for H in kHz/T = (cycles/millisecond)/T
nSlices = 65;

maxB = max(bvals);

% Sort to ensure that we have the lowest bvalue first, and then all
% remaining bvlaues in descending order. This will make the interleaving
% (below) work better at keeping high bvalue scans separated from other
% high b-value scans.
[junk,si] = sort(bvals,'descend');
%si = si([end,1:end-1]);
bvals = bvals(si);
nDirs = nDirs(si);
nReps = nReps(si);

if(~exist('delta','var'))
    % DW gradient duration, in msec
    delta = []; % We'll set it below
end

if(~exist('Delta','var') || isempty(Delta))
    if(~isempty(delta))
        Delta = maxB/((2*pi*g).^2 * (maxG*1e-9).^2 * delta.^2) + delta/3;
    else
        % Find the optimal Delta for the max b-value (i.e. where Delta=delta).
        % This ugly equation is simply the solution for this equation:
        %   Delta = max(bvals)/((2*pi*g).^2 * (maxG*1e-9).^2 * delta.^2) + delta/3;
        % when we set delta=Delta and solve for Delta again.
        Delta = 2*46874999999999994^(1/3)*maxB^(1/3)/(pi^(2/3)*(g^2)^(1/3)*(maxG^2)^(1/3));
    end
end

if(isempty(delta))
    % DW gradient duration, in msec
    delta = Delta;
end

% On our GE Signa, there is ~81msec of overhead for RF, imaging grads,
% read-out, etc.
overhead = 81;
sliceTime = Delta+delta+overhead;
TR = sliceTime/1000*nSlices;
nVols = sum(nReps.*nDirs);
totalTime = TR*nVols/60;
if interact==1
fprintf('Estimated TR for %d slices = %0.2f. Total scan time for %d volumes is %0.1f minutes.\n\n',nSlices,TR,nVols,totalTime);
end;

% Now compute the Gradient amplitudes for all the b-values
G = sqrt(bvals./((2*pi*g).^2 * delta.^2 * (Delta-delta/3))) * 1e9;

ptsDir = fullfile(fileparts(mfilename('fullpath')),'caminoPts')

n = numel(bvals);
if interact==1
figure; axis; grid on; axis equal; hold on; c = 'rgbycmk';
end;
grads = repmat(NaN,3,nVols);
availInds = [1:nVols];
for(ii=1:n)
    if(nDirs(ii)>=3)
        pts = dlmread(fullfile(ptsDir,sprintf('Elec%03d.txt',nDirs(ii))));
        pts = reshape(pts(2:end),[3 nDirs(ii)]);
    elseif(nDirs(ii)==2 || nDirs(ii)<1 || nDirs(ii)>150)
        error('No supported');
    else
        % nDirs must == 1
        pts = [1;0;0];
    end
    curGrads = G(ii)./maxG .* pts;
    curGrads = repmat(curGrads,[1 nReps(ii)]);
    % flip half the directions
    if(bvals(ii)>0)
        flipThese = rand(1,size(curGrads,2))>=0.5;
        curGrads(:,flipThese) = -curGrads(:,flipThese);
    end
    inds = floor([1:numel(availInds)/size(curGrads,2):numel(availInds)]);
    grads(:,availInds(inds)) = curGrads;
    availInds = setdiff(availInds,availInds(inds));
end

return;
