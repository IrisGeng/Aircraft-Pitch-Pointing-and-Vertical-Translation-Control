clear all;
% Define state
A=[-1.3,0.98,0,-0.165,-0.248;42.81,-0.785,0,-17.3,-1.58;1.25,0.007,0,0.165,0.248;0,0,0,-18,0;0,0,0,0,-18];
B=[0,0;0,0;0,0;18,0;0,18];
C=[0,1,0,0,0;46.5,-0.256,0,-4.25,4.15;0,0,1,0,0];
D=[0,0;0,0;0,0];
[U,V]=eig(A)

% observable,controllable for MIMO
Mc=[B,A*B,A^2*B,A^3*B,A^4*B];
rank(Mc);
Mo=[C',A'*C',A'^2*C',A'^3*C',A'^4*C'];
rank(Mo);

%% 
% Mc for input 1
B1=B(:,1);
C1=C(1,:);
Mc1=[B1,A*B1,A^2*B1,A^3*B1,A^4*B1];
rank(Mc1); %=4

% Mc for input 2
B2=B(:,2);
Mc2=[B2,A*B2,A^2*B2,A^3*B2,A^4*B];
rank(Mc2); %=5

% Mo for output 1
Mo1=[C1',A'*C1',A'^2*C1',A'^3*C1',A'^4*C1'];
rank(Mo1); %=3

% Mo for output 2
C2=C(2,:);
Mo2=[C2',A'*C2',A'^2*C2',A'^3*C2',A'^4*C2'];
rank(Mo2); %=3

% Mo for output 3
C3=C(3,:);
Mo3=[C3',A'*C3',A'^2*C3',A'^3*C3',A'^4*C3'];
rank(Mo3); %=4



