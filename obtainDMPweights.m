function [weights,dmpout,f] = obtainDMPweights(dmpstruct, K, position, options)
% OBTAINDMPWEIGHTS Obtain DMP weights from input data.
%     weights = obtainDMPweights(DMP, K, pos)
%     Recommendation: use dmp output to update DMP input.
%     [weights,DMP] = obtainDMPweights(DMP, K, pos)
% 
%     Arguments:
%     DMP is the initialised DMP structure. Must be initialised with the
%     correct data collection frequency and data length.
%     K is the stiffness parameter value in DMP.
%     pos is an n x d array of the position values.
% 
%     Example:
%         dmpout = initialiseDMP(350) % prereq to obtain initialised dmp
%         wgt = obtainDMPweights(dmpout, 80, dataposition)
%         % Example to obtain forcing function values
%         [wgt,dmpout,forcing] = obtainDMPweights(dmpout, 80, dataposition)
% 
%     weights = obtainDMPweights(... , options) specifies other parameters
%     and/or velocity and acceleration data. Other optional parameters are:
% 
%     'D'           - Damping coefficient. (default = K/4)
%     'velocity'    - an n x d array of the velocities. Note that the array
%                     size has to be the same as the position array.
%                     (Default is using integration from position data. 
%                     Last sample repeats)
%     'acceleration'- an n x d array of the accelerations. Note that the 
%                     array size has to be the same as the velocity array.
%                     (Default is using integration from velocity data.
%                     Last sample repeats)
% 
%     Example:
%         % Example for different damping coefficient
%         dmpout = initialiseDMP(350);
%         wgt = obtainDMPweights(dmpout, 80, dataposition, 'D', 15)
%         % Example for original velocity and acceleration data
%         [wgt,~,f] = obtainDMPweights(dmpout, 80, dataposition, 'velocity',
%                     datavelocity,'acceleration', dataacceleration)
% 
% 	For more information on this particular implementation, refer to <a
%     href="https://link.springer.com/chapter/10.1007/11008941_60">this article</a>.

arguments
    dmpstruct
    K                    {mustBeNumeric, mustBeNonNan, mustBeNonnegative}
    position             {mustBeNonNan, mustBeFinite}
    options.D            {mustBeNumeric, mustBeNonNan, mustBeNonnegative} = []
    options.velocity     {mustBeNonNan, mustBeFinite} = []
    options.acceleration {mustBeNonNan, mustBeFinite} = []
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
dmpout = dmpstruct; dmpout.D = D; dmpout.K = K; 
% Calculate forcing function from data
f = (dmpstruct.tauval^2*accel+D*dmpstruct.tauval*vel) ...
    ./(K*repmat(dmpstruct.phase,1,size(position,2))) ...
    -(position(end,:)-position) ...
    ./repmat(dmpstruct.phase,1,size(position,2)) ...
    +position(end,:)-position(1,:);
% Linear Least Squares to obtain DMP weights
weights = (dmpstruct.psi'*f)./(dmpstruct.psi'*dmpstruct.phase+1e-5);
end