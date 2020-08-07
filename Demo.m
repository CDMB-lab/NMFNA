clear all;
load('ME.mat');
load('CN.mat');

% load('PAAD_ME.mat');
% load('PAAD_CN.mat');
% ME=PAAD_ME;
% CN=PAAD_CN;
%% Calculate Pearson coefficient matrix
R11=corr(ME','type','pearson');
R22=corr(CN','type','pearson');
R12=corr(ME',CN','type','pearson');
X1 = mapminmax(R11, 0, 1);
X1(isnan(X1)==1) = 0;
X2 = mapminmax(R22, 0, 1);
X2(isnan(X2)==1) = 0;
X3 = mapminmax(R12, 0, 1);
X3(isnan(X3)==1) = 0;

%% One-step finding communities of PAAD
 % maxiter = 200;
 % options.NeighborMode = 'KNN';
 % options.k = 5;
 % options.WeightMode = 'Binary';  
 % options.t = 1;
 % W1 = constructW(R11',options);  
 % W2 = constructW(R22',options); 
 % r=6834;
 % lamda1 = 0.04;
 % lamda2=lamda1;

 % NMFNA function
 % tic
 % [S11,S22,G1,G2,res,errorx] = NMFNA(X1, X2, X3,W1,W2, r, maxiter,lamda1,lamda2);
 % toc

 % identification of community
 % xita=2;
 % [ communityFinal1,mSim1] = communitySelection( G1, xita );
 % [ communityFinal2,mSim2 ] = communitySelection( G2, xita );

 % mSim1(isnan(mSim1)) = 0;
 % mSim2(isnan(mSim2)) = 0;
 % mSim=mSim1+mSim2;

%% Step-by-step finding communities
 % finding the K
 % [U,S,V]=svd(X1);
 % [m,n]=size(S);
 % SS=ones(n,1); 
 % WQ=S*SS;
 % plot(WQ);
 % A=WQ(1:n);
 % t1 = [1 : size(A,1)];
 % p1 = polyfit(t1,A',10);
 % Y1 = polyval(p1,t1);
 % DY1 = diff(Y1,2);
 % plot(t1,Y1,'b-',t1,A','r-')

 % initilaize parameters 
 % maxiter = 200;
 % options.NeighborMode = 'KNN';
 % options.k = 5;
 % options.WeightMode = 'Binary';  
 % options.t = 1;
 % W1 = constructW(R11',options);  
 % W2 = constructW(R22',options); 

 % selecting lambda
 % lamda1 = 0.01:0.01:0.1;
 % lamda2 =  lamda1;
 % for i=1:size(lamda1,2)
 %  [S11,S22,G1,G2,res,errorx] = NMFNA(X1, X2, X3,W1,W2, r, maxiter,lamda1(i),lamda2(i));
 %   xita=2;
 %  [ communityFinal1,mSim1] = moduleNodesSelection( G1, xita );
 %  [ communityFinal2,mSim2 ] = moduleNodesSelection( G2, xita );
 %  mSim(i,j)=mSim1+mSim2;
 % end
%% initilaize parameters 
maxiter = 200;
options.NeighborMode = 'KNN';
options.k = 5;
options.WeightMode = 'Binary';  
options.t = 1;
W1 = constructW(R11',options);  
W2 = constructW(R22',options); 
r=10;
lamda1 = 0.04;
lamda2=lamda1;

%% NMFNA function
tic
[S11,S22,G1,G2,res,errorx] = NMFNA(X1, X2, X3,W1,W2, r, maxiter,lamda1,lamda2);
toc

%% identification of community
xita=2;
[ communityFinal1,mSim1] = communitySelection( G1, xita );
[ communityFinal2,mSim2 ] = communitySelection( G2, xita );

mSim1(isnan(mSim1)) = 0;
mSim2(isnan(mSim2)) = 0;
mSim=mSim1+mSim2;


    