function [ coil_position, coil_dir ] = find_my_beans_points( varargin )
% BSD 2-Clause License

% Copyright (c) 2019, fmri-at
% All rights reserved.

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% parse the input
p = inputParser;
validPoints = @(x) all(size(x) == [3 3]);
addOptional(p, 'points', nan*ones(3), validPoints);
addParameter(p, 'nii', '', @(x) ischar(x));
addParameter(p, 'use_world_space', false, @(x) islogical(x));
addParameter(p, 'create_niftis', true, @(x) islogical(x));

parse(p, varargin{:});
r = p.Results;

% see if there is a voxel to world transformation present from a nifti (use
% SPM's nifti class)
nifti_exists = false;

if ~isempty(r.nii)
    if exist(r.nii, 'file')
        nii = nifti(r.nii);
        mat = nii.mat;
        sz = nii.dat.dim;
        nifti_exists = true;
    end
end
if ~r.use_world_space
    if ~nifti_exists
        error('You have to supply a nifti if you want to use voxel space.')
    end
end

points = r.points;

if any(isnan(points))
    n_points = 0;
    while n_points < 3
        i = n_points + 1;
        
        if ~r.use_world_space
            p_i = -1;
            p_j = -1;
            p_k = -1;
            fprintf('Please supply the zero-based voxel coordinates of one bean:\n')

            while ~isnumeric(p_i) || ~isscalar(p_i) || p_i < 0 || p_i > (sz(1) - 1)
                p_i = input('i: ');
            end

            while ~isnumeric(p_j) || ~isscalar(p_j) || p_j < 0 || p_j > (sz(2) - 1)
                p_j = input('j: ');
            end

            while ~isnumeric(p_j) || ~isscalar(p_j) || p_k < 0 || p_k > (sz(3) - 1)
                p_k = input('k: ');
            end

            p_add = [p_i p_j p_k] + 1;
            p_add_world = mat * [p_add' ; 1];
            valid = input(sprintf('Is the input [%u, %u, %u] (%f, %f, %f) correct? [Y,n]: ', p_i, p_j, p_k, p_add_world(1), p_add_world(2), p_add_world(3)), 's');
            if isempty(valid)
                valid = 'y';
            end

            if ~isempty(regexpi(strtrim(valid), '^y.*', 'match'))
                points(:,i) = p_add_world(1:3);
                n_points = n_points + 1;
            end
        else
            fprintf('Please supply the world coordinates of one bean:\n')
            p_x = nan;
            p_y = nan;
            p_z = nan;
            
            while ~isnumeric(p_x) || ~isscalar(p_x) || isnan(p_x)
                p_x = input('x: ');
            end
            
            while ~isnumeric(p_y) || ~isscalar(p_y) || isnan(p_y)
                p_y = input('y: ');
            end
            
            while ~isnumeric(p_z) || ~isscalar(p_z) || isnan(p_z)
                p_z = input('z: ');
            end
            
            p_add_world = [p_x p_y p_z].';
                        
            valid = input(sprintf('Is the input (%f, %f, %f) correct? [Y,n]: ', p_add_world), 's');
            if isempty(valid)
                valid = 'y';
            end

            if ~isempty(regexpi(strtrim(valid), '^y.*', 'match'))
                points(:,i) = p_add_world;
                n_points = n_points + 1;
            end
        end
    end
end
    
coil_dir = cross(points(:,1) - points(:,2), points(:,1) - points(:,3)); % calculate normal to spheres
coil_dir = coil_dir / sqrt(coil_dir' * coil_dir);
    
coil_position = mean(points, 2); % mean of sphere coordinates lies on stimulation point

if r.create_niftis
    if nifti_exists
        [ii, jj, kk] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
        [path, name] = fileparts(r.nii);
        
        % create nifti with visible path
        nii.dat.fname = fullfile(path, ['tms_path_' name '.nii']);
        nii.dat.dtype  = 'FLOAT32-LE';
        create(nii)
        
        pp = mat * [ ii(:)' ; jj(:)' ; kk(:)' ; ones(1,numel(ii)) ]; % transform every voxel coordinate to nifti coordinates
        
        coil_dir_hom = [coil_dir ; 0];
        coil_position_hom   = [coil_position ; 1];
        
        proj = coil_dir_hom.' * (coil_position_hom(:,ones(1,numel(ii))) - pp);      % calculate distance from line of every voxel
        dd_v = (coil_position_hom(:,ones(1,numel(ii))) - pp) - proj(ones(1,4),:) .* coil_dir_hom(:,ones(1,numel(ii)));
        dd =  dd_v(1,:).^2 + dd_v(2,:).^2 + dd_v(3,:).^2;
    
        nii.dat(:,:,:) = reshape(exp(-dd/4), size(ii)); % use gaussian like function to visualise stimulation path
    
    end
end

