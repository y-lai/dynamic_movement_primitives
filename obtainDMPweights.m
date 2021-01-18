function [weights,f] = obtainDMPweights(dmpstruct,K,position,options)
% OBTAINDMPWEIGHTS Obtain DMP weights from input data
%     weights = obtainDMPweights(DMP,K,pos)
%     DMP is the initialised DMP structure. Must be initialised with the
%     correct data collection frequency and data length.
%     K is the stiffness parameter value in DMP.
%     pos is an n x d array


arguments
    dmpstruct
    K {mustBeNonNan,mustBeFinite}
    position {mustBeNonNan}
    options.D {mustBeNonNan,mustBeFinite} = []
    options.velocity {mustBeNonNan} = []
    options.acceleration {mustBeNonNan} = []
end
% optional name value pairs
if(isempty(options.velocity))
    vel = diff(position)*dmpstruct.datafreq;
    vel = [vel; vel(end,:)];
end
if(isempty(options.acceleration))
    accel = diff(vel)*dmpstruct.datafreq;
    accel = [accel; accel(end,:)];
end
if(isempty(options.D))
    D = K/4;
else
    D = options.D;
end
if (sum(size(position)~=size(vel))+sum(size(position)~=size(accel))>0)
    error(['Invalid size of input velocity and acceleration arrays. ' ...
        'Please make sure the size of the position, velocity, ' ...
        'and acceleration arrays are the same.']);
end
% Calculate forcing function from data
f = (dmpstruct.tauval^2*accel+D*dmpstruct.tauval*vel) ...
    ./(K*repmat(dmpstruct.phase,1,size(position,2))) ...
    -(position(end,:)-position) ...
    ./repmat(dmpstruct.phase,1,size(position,2)) ...
    +position(end,:)-position(1,:);
% Linear Least Squares to obtain DMP weights
weights = (dmpstruct.psi'*f)./(dmpstruct.psi'*dmpstruct.phase+1e-5);
end