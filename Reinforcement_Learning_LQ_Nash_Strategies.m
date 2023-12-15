% Reinforcement Learning Final Term Paper
% [Applications of Nash Differential Games to Aerospace]

%Clear All Variables in our Workspace
clear all

%Set value of A
A = [-0.4125 -0.0248 0.0741 0.0089 0 0;
     101.5873 -7.2651 2.7608 2.8068 0 0;
     0.0704 0.0085 -0.0741 -0.0089 0 0.02;
     0.0878 0.2672 0 -0.3674 0.0044 0.3962;
     -1.8414 0.099 0 0 -0.0343 -0.033;
     0 0 0 -359 187.5364 -87.0316];

%Set value of B1
B1 = [-0.0042 -1.036 0.0042 0.1261 0 0;
      0.0064 1.5849 0 0 -0.0168 0]';

%Set value of B2
B2 = [-0.0153 -1.526 0.0032 0.3261 0 0;
      0.0364 1.7 0 0 -0.2 0]';

eig(A)

%Get values of Q matrices
Q1=4*eye(6);
Q2=10*eye(6);

%Get values of R matrices
R11=2*eye(2);
R12=zeros(2,2);
R21=zeros(2,2);
R22=5*eye(2);


S1=B1*inv(R11)*B1';
S2=B2*inv(R22)*B2';


% initialization
K10=are(A,S1,Q1);
K20=are(A-S1*K10,S2,Q2+K10*S1*K10);
K100=K10;
K200=K20;


i=0;
J10=0.5*trace(K100);
J20=0.5*trace(K200);
% eig(A-S1*K10-S2*K20)

for i=1:10
    K1i=lyap2((A-S1*K10-S2*K20)',Q1+K10*S1*K10+K20*S2*K20);
    K2i=lyap2((A-S1*K10-S2*K20)',Q2+K20*S2*K20+K10*S1*K10);
    K10=K1i;
    K20=K2i;
    i
    J1i=0.5*trace(K1i);
    J2i=0.5*trace(K2i);
end


% Performance found under the assumption that the initial conditions are
% uniformly distributed on the unit sphere so that they do not affect the
% optimal performance
