clear
clc
load('data.mat');

% initialise dmp
dmp = initialiseDMP(length(learnPos));
% learn dmp weights
[weights,dmp] = obtainDMPweights(dmp,params.K,learnPos);
% [weights,dmp] = obtainDMPweights(dmp,130,learnPos);
% generate trajectory using DMP weights with different translation
percent = 5; 
tPos = testPos(1:floor(percent/100*length(testPos)),:);
% add translation difference
tPos = tPos + [0.2 -0.12]; refPos = refPos + [-0.15 0];
[pos,vel,accel] = runDMPweights(dmp,length(refPos),weights,tPos,refPos);

%% plot figures
figure; hold on;
plot(learnPos(:,1),learnPos(:,2),'r'); plot(pos(:,1),pos(:,2),'b');
plot(testPos(:,1),testPos(:,2),'k'); plot(refPos(:,1),refPos(:,2),'m');
legend({'Learning','DMP generated','observed','Reference'},'Location','northwest');