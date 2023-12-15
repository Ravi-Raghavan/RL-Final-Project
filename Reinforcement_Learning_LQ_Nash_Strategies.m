% Reinforcement Learning Final Term Paper
% [Applications of Nash Differential Games to Aerospace]

%Clear All Variables in our Workspace
clear all

%Initialize the state value
x = ones(1, 6)';

%Note: Currently, values of A, B1, B2, B3, Q Matrices, R matrices were set
%to bullshit values. This needs to be verified b4 submission!

%Set value of A
A = [-0.4125 -0.0248 0.0741 0.0089 0 0;
     101.5873 -7.2651 2.7608 2.8068 0 0;
     0.0704 0.0085 -0.0741 -0.0089 0 0.02;
     0.0878 0.2672 0 -0.3674 0.0044 0.3962;
     -1.8414 0.099 0 0 -0.0343 -0.033;
     0 0 0 -359 187.5364 -87.0316];

%Set value of B1
B1 = [-0.0042 -1.036 0.0042 0.1261 0 0;
      0.0064 1.5849 0 0 -0.0168 0;
      0.0064 1.5849 0 0 -0.0168 0]';

%Set value of B2
B2 = [-0.0153 -1.526 0.0032 0.3261 0 0;
      0.0364 1.7 0 0 -0.2 0;
      0.0064 1.5849 0 0 -0.0168 0]';

%Set value of B3
B3 = [-0.0153 -1.526 0.0032 0.3261 0 0;
      0.0364 1.7 0 0 -0.2 0;
      0.0064 1.5849 0 0 -0.0168 0]';

%Get values of Q matrices
Q1=4*eye(6);
Q2=10*eye(6);
Q3=10*eye(6);

%Get values of R matrices
R11=2*eye(3);
R12=eye(3);
R13=eye(3);
R21=eye(3);
R22=5*eye(3);
R23=eye(3);
R31=5*eye(3);
R32=eye(3);
R33=eye(3);

%Set values of S based on B and R matrices
S1=B1*inv(R11)*B1';
S12=B1*inv(R11)*R21*inv(R11)*B1';
S13=B1*inv(R11)*R31*inv(R11)*B1';
S2=B2*inv(R22)*B2';
S21=B2*inv(R22)*R12*inv(R22)*B2';
S23=B2*inv(R22)*R32*inv(R22)*B2';
S3=B3*inv(R33)*B3';
S31=B3*inv(R33)*R13*inv(R33)*B3';
S32=B3*inv(R33)*R23*inv(R33)*B3';

%Solve for Initial Iterative Matrices P1, P2, and P3
P1 = are(A, S1, Q1);
P2 = are(A-S1*P1, S2, Q2 + P1*S12*P1);
P3 = are(A-S1*P1-S2*P2, S3, Q3 + P1*S13*P1 + P2*S23*P2);

iterations = 10;

%Store P1, P2, P3 Values to see if there is convergence
P1_values = cell(1, iterations + 1);
P2_values = cell(1, iterations + 1);
P3_values = cell(1, iterations + 1);

P1_values{1} = P1;
P2_values{1} = P2;
P3_values{1} = P3;

for i = 1:iterations
    P1_updated = lyap2((A-S1*P1-S2*P2-S3*P3)',Q1+P1*S1*P1+P2*S21*P2+P3*S31*P3);
    P2_updated = lyap2((A-S1*P1-S2*P2-S3*P3)',Q2+P1*S12*P1+P2*S2*P2+P3*S32*P3);
    P3_updated = lyap2((A-S1*P1-S2*P2-S3*P3)',Q2+P1*S13*P1+P2*S23*P2+P3*S3*P3);

    P1 = P1_updated;
    P2 = P2_updated;
    P3 = P3_updated;

    P1_values{i + 1} = P1_updated;
    P2_values{i + 1} = P2_updated;
    P3_values{i + 1} = P3_updated;
end


