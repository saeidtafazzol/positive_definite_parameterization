clc,clear,close all

%% Choose Parameters

eig_vals = [3 5 2]; % eigenvalues (its length determines dimension of positive-definite matrix)
thetas = [0.5 0.3 0.2]; % angles (there should be (N^2-N)/2 number of angles where N = length(eig_vals))

%% Call Function to Get Positive-Definite Matrix, K

[K,L,Q] = get_positive(eig_vals,thetas);
disp(K)

%% Subfunctions

function [K,L,Q] = get_positive(eig_vals, thetas)
    N = length(eig_vals); % determine dimension of matrix based on number of eigenvalues

    if numel(thetas) ~= (N^2-N)/2 % check to ensure that correct number of angle parameters were chosen
        error('incorrect number of angles')
    end

    thetas_mat = zeros(N-1, N-1); % initialize matrix of angles

    starting_idx = 1; % initialize indexing

    for i = 1:N-1 % loop that constructs matrix of angles
        thetas_mat(i,1:N-i) = thetas(starting_idx:starting_idx + N - i - 1);
        starting_idx = starting_idx + N - i;
    end

    Q = get_orthogonal(thetas_mat); % get matrix of eigenvectors (N-dimensional rotation matrix)

    L = diag(eig_vals); % construct eigenvalue matrix
    
    K = Q * L * Q'; % multiply following eigendecompisition
end

function Q = get_orthogonal(thetas)

    Q = get_v(thetas(1,:)); % construct first eigenvector

    for i = 2:size(thetas,2) % loop to recursively build Q
        B = null(Q'); % get nullspace of Q
        v = B * reshape(get_v(thetas(i,1:end-i+1)), [], 1); % project back into N space
        Q = [Q v]; % append new eigenvector to Q
    end

    B = null(Q'); % get last eigenvector
    Q = [Q B]; % append
end

function v = get_v(thetas) % function to construct eigenvectors (N-dimensional unit vector in spherical coordinates)
    if isempty(thetas)
        v = 1.0;
        return;
    end
    vp = get_v(thetas(2:end));
    v = [cos(thetas(1)); sin(thetas(1)) * vp]; % see Eq. 4 in paper
end
