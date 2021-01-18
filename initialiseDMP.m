function dmpout = initialiseDMP(datalength,options)
% INITIALISEDMP Initialise a DMP structure.
%     DMP = initialiseDMP(X) creates a DMP structure.
%     X is a scalar value for the length of the data.
%     
%     Example:
%         dmpout = initialiseDMP(500)
%
%     DMP = initialiseDMP(...,options) specifies other parameters for ...
%     the initialisation of DMP. Other optional parameters are:
%     
%     'phase_alpha' - alpha parameter for DMP. (default = 0.8)
%     'numrbf'      - Number of Radial Basis Function (default = 30)
%     'width'       - Width of RBFs. (default = 0.55)
%     'datafreq'    - Collection frequency of the dataset. (default = 100)
%     'tauval'      - Tau parameter for DMP. (default utilises exp2 fit from
%                     phrip article)
% 
%     Example:
%         % Example for data collected at 125 hz
%         dmpout = initialiseDMP(... , 'datafreq' , 125)
%         % Example for different number of RBF (50) for forcing function
%         dmpout = initialiseDMP(... , 'numrbf' , 50)
%
%     For more information on DMP, refer to <a
%     href="https://link.springer.com/chapter/10.1007/11008941_60">this article</a>.
%     For more information on the tau calculation, refer to <a
%     href="https://www.researchgate.net/publication/348575169_User_Intent_Estimation_during_robot_learning_using_Physical_Human_Robot_Interaction_Primitives">this article</a>.

arguments
    datalength          int64
    options.tauval      {mustBeNumeric, mustBeNonNan, mustBeNonnegative} = []
    options.phase_alpha {mustBeNumeric, mustBeNonNan, mustBeNonnegative} = 0.8
    options.numrbf      {mustBeNumeric, mustBeNonNan, mustBeNonnegative} = 30
    options.width       {mustBeNumeric, mustBeNonNan, mustBeNonnegative} = 0.55
    options.datafreq    {mustBeNumeric, mustBeNonNan, mustBeNonnegative} = 100
end
% assign default values
dmpout.datalength = datalength;
% optional name value pairs
dmpout.numrbf = options.numrbf; dmpout.width = options.width;
dmpout.datafreq = options.datafreq; dmpout.alpha = options.phase_alpha;
if(isempty(options.tauval))
    % exp2 fit for tauval is based on data length. Refer to phrip article.
    % link to arxiv/preprint
    % Values ensure sum of forcing function activation are above 0.5
    dmpout.tauval = 2.22*exp(-0.01292*cast(dmpout.datalength,'double')) ...
        + 0.4684*exp(-0.001595*cast(dmpout.datalength,'double'));
else
    if(length(options.tauval)>1)
        error('Incorrect size for argument. Optional argument tauval must be a single scalar value');
    end
    dmpout.tauval = options.tauval;
end

dmpout.centers = exp(-dmpout.alpha*(linspace(0,1,dmpout.numrbf)'));
X = zeros(dmpout.datalength,1); dX = X; x = 1; dt = 1/dmpout.datafreq;
for i=1:dmpout.datalength
    X(i) = x; dX(i) = -dmpout.alpha*x*dmpout.tauval; x = x + dX(i)*dt;
end
dmpout.phase = X; d = (diff(dmpout.centers)*dmpout.width).^2;
d = 1./[d; d(end)]; dmpout.psi = exp(-0.5*((X*ones(1,length(dmpout.centers)) ...
    -ones(dmpout.datalength,1)*dmpout.centers').^2) ...
    .*(ones(dmpout.datalength,1)*d'));
end