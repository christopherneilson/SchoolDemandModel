
function [ Q_NODES, Q_WEIGHTS ] = GHQuadInit( nDim_, nNodes_ )


% Get one-dimensional Gauss-Hermite nodes and weights

[vNodes,vWeights]=gauher(nNodes_);
% extract quadrature information for one dimension
%vNodes = x ;
%vWeights = w ;
% calculate three dimensional nodes and weights
% Q_WEIGHTS = kron( vWeights, kron( vWeights, vWeights ) ) ;
Q_WEIGHTS = vWeights ;
for ix = 2 : nDim_
Q_WEIGHTS = kron( vWeights, Q_WEIGHTS ) ;
end
% Make sure that the right-most dimension (ixDim = nDim_) varies
% most quickly and the left-most (ixDim = 1) most slowly
Q_NODES = zeros( nDim_, nNodes_^nDim_ ) ;
for ixDim = 1 : nDim_
Q_NODES( ixDim, : ) = kron( ones( nNodes_^(ixDim - 1), 1 ), ...
kron( vNodes, ones( nNodes_^(nDim_ - ixDim), 1 ) ) ) ;
end
% Correct for Gaussian kernel versus normal density
%Q_WEIGHTS = Q_WEIGHTS; %* ( pi ^ ( nDim_ / 2 ) ) ;
%Q_NODES = Q_NODES; %*sqrt( 2 ) ;