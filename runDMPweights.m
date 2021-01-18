function [pos,vel,accel,f] = runDMPweights(dmpstruct, datalength, weights, obs, ref, options)
% RUNDMPWEIGHTS Obtain trajectory from weights, observation and reference
%     pos = runDMPweights(DMP,wgt,startpos,refpos)
%     [pos,vel,accel] = runDMPweights(DMP,datanum,wgt,startpos,refpos)
% 
%     Arguments:
%     DMP is the output DMP structure from obtainDMPweights. Require K and
%     D value from obtaining weights. If new K and D required, use options.
%     datanum is length of trajectory.
%     wgt are the m x 1 array of weights to run DMP trajectory generation.
%     startpos is an n x d array of observed position (minimum 1 sample).
%     refpos is an c x d array of reference position (minimum 1 sample) .
% 
%     Example:
%         dmpout = initialiseDMP(350)
%         [wgt, dmpout] = obtainDMPweights(dmpout, 80, dataposition)
%         [pos,~,~,f] = runDMPweights(dmpout, wgt, dataobs, refdata)
% 
%     pos = runDMPweights(... , options) specifies other data input and
%           parameters. Optional parameters are:
% 
%     'obsvel'      - l x d array of observed velocities if available. Note
%                     that the array size has to be the same as the obs array.
%     'obsaccel'    - l x d array of observed acceleration if available. Note
%                     that the array size has to be the same as the obs array.
%     'phase_alpha' - alpha parameter for DMP. (default = dmpstruct.alpha)
%     'numrbf'      - Number of Radial Basis Function. (default = dmpstruct.numrbf)
%     'width'       - Width of RBFs. (default = dmpstruct.width)
%     'datafreq'    - Collection frequency of the dataset. (default = dmpstruct.datafreq)
%     'K'           - Stiffness coefficient. (default is dmpstruct.K)
%     'D'           - Damping coefficient. (default is dmpstruct.D)
%     'tauval'      - Tau parameter for DMP. (default utilises exp2 fit 
%                     from phrip article)
% 
%     pos = runDMPweights(..., 'obsvel', datavel, 'obsaccel', dataaccel)
% 
%     For more information on DMP, refer to <a 
%     href="https://link.springer.com/chapter/10.1007/11008941_60">this article</a>.
%     For more information on the tau calculation, type 'help initialiseDMP'

arguments
    dmpstruct
    datalength          int64
    weights             {mustBeNumeric, mustBeNonNan, mustBeFinite}
    obs                 {mustBeNumeric, mustBeNonNan, mustBeFinite}
    ref                 {mustBeNumeric, mustBeNonNan, mustBeFinite}
    options.obsvel      {mustBeNumeric, mustBeNonNan, mustBeFinite} = []
    options.obsaccel    {mustBeNumeric, mustBeNonNan, mustBeFinite} = []
    options.phase_alpha {mustBeNumeric, mustBeNonNan, mustBeNonnegative} = []
    options.numrbf      {mustBeNumeric, mustBeNonNan, mustBeNonnegative} = []
    options.width       {mustBeNumeric, mustBeNonNan, mustBeNonnegative} = []
    options.datafreq    {mustBeNumeric, mustBeNonNan, mustBeNonnegative} = []
    options.K           {mustBeNumeric, mustBeNonNan, mustBeNonnegative} = []
    options.D           {mustBeNumeric, mustBeNonNan, mustBeNonnegative} = []
    options.tauval      {mustBeNumeric, mustBeNonNan, mustBeNonnegative} = []
end
% optional name value pairs
if(isempty(options.phase_alpha))
    alpha = dmpstruct.alpha;
else
    alpha = options.phase_alpha;
end
if(isempty(options.numrbf))
    numrbf = dmpstruct.numrbf;
else
    numrbf = options.numrbf;
end
if(isempty(options.width))
    width = dmpstruct.width;
else
    width = options.width;
end
if(isempty(options.datafreq))
    datafreq = dmpstruct.datafreq;
else
    datafreq = options.datafreq;
end
if(isempty(options.K))
    K = dmpstruct.K;
else
    K = options.K;
end
if(isempty(options.D))
    D = K/4;
else
    D = options.D;
end
if(isempty(options.tauval))
    tauval = [];
else
    tauval = options.tauval;
end
dmp = initialiseDMP(datalength, 'phase_alpha', alpha, 'numrbf', numrbf, ...
    'width', width, 'datafreq', datafreq, 'tauval', tauval);
if(isempty(options.obsvel))
    obsvel = diff(obs)*dmp.datafreq;
    obsvel = [obsvel; obsvel(end,:)];
else
    obsvel = options.obsvel;
end
if(isempty(options.obsaccel))
    obsaccel = diff(obsvel)*dmp.datafreq;
    obsaccel = [obsaccel; obsaccel(end,:)];
else
    obsaccel = options.obsaccel;
end
if (sum(size(obs)~=size(obsvel))+sum(size(obs)~=size(obsaccel))>0)
    error(['Invalid size of input velocity and acceleration arrays. ' ...
        'Please make sure the size of the position, velocity, ' ...
        'and acceleration arrays are the same.']);
end
% Obtained forcing function from weights
f = ((dmp.psi*weights).*repmat(dmp.phase,1,size(weights,2))) ...
    ./repmat((sum(dmp.psi,2)+1e-3),1,size(weights,2));
refend = ref(end,:); dt = 1/dmp.datafreq;
% Integration of acceleration to velocity and position
pos = zeros(datalength,size(weights,2)); vel = pos; accel = pos;
pos(1:length(obs),:) = obs; vel(1:length(obs),:) = obsvel;
accel(1:length(obs),:) = obsaccel;
for i=length(obs)+1:datalength
    pos(i,:) = pos(i-1,:) + vel(i-1,:)*dt;
    vel(i,:) = vel(i-1,:) + accel(i-1,:)*dt;
    accel(i,:) = (K*(refend-pos(i,:))+dmp.phase(i)*K*(f(i,:)- ...
        (refend-pos(1,:))))/dmp.tauval^2 - (D*vel(i,:)/dmp.tauval);
end
end