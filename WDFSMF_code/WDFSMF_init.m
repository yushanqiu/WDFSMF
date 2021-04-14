clear all;
funspath=[pwd,filesep,'data',filesep];
addpath(funspath);
%% load data to construct the relation matrix R
load('RbpGene');
R13 = Re_13;
load('RbpAS');
R14 = Re_14;
load('MiDOs')
R25 = miDOs;
load('GeneAS')
R34 = Re_34;
load('miRNAGene')
R23 = MGasso;
load('GeneDisease')
R35 = GDasso;
[npro,ndi] = size(GDasso);
load('GeneDrug')
R63 = GDrgasso;
nDrug = size(R63,1);
fun_stat = sum(R14,1);
sel_do_idx = find(fun_stat>0);
R14 = R14(:,sel_do_idx);
R34 = R34(:,sel_do_idx);
ndi = length(sel_do_idx);
Rcell={R13,R14,R23,R25,R34,R35,R63};

k1=100;
k2=100;
k3=170;
k4=100;
k5=170;
k6=50;
k = {k1,k2,k3,k4,k5,k6};   
   %low-rank dim
%R1 = [R13,R14,R61'];
R1 = [R13];

G1=abs(rand(size(R1,1),k1));
G2=abs(rand(size(R23,1),k2));
G3=abs(rand(size(R13,2),k3));
G4=abs(rand(size(R34,2),k4));
G5=abs(rand(size(R35,2),k5));
G6=abs(rand(size(R63,1),k6));

Gcell={G1,G2,G3,G4,G5,G6};





load('HumanPPI')

load('DrugInter')
thetaCell=cell(6,1);
thetaCell{3}=PPI.*-1;
thetaCell{6}=DDI.*-1;

%instanseIdx = {3,4,5,9,11,16,17,31};  
instanseIdx = {3,4,9,11,16,17,33};% instanse index for the corresponding position in the whole relation matrix R
fprintf('begin,time:%s\n',datestr(now));
t0=clock;
nTypes=6;
newF_RBP_AS_auto_sparse = WDFSMF_demo(nTypes,instanseIdx,Gcell,Rcell,thetaCell);
save newF_RBP_AS_auto_sparse newF_RBP_AS_auto_sparse;


