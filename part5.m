% part 5

AN=[-1.3,0.98,0,-0.165,0;42.81,-0.785,0,-17.3,0;1.25,0.007,0,0.165,0;0,0,0,-18,0;0,0,0,0,0];
BN=[0,0;0,0;0,0;18,0;0,0];
CN=[0,1,0,0,0;46.5,-0.256,0,-4.25,0;0,0,1,0,0];
DN=[0,0;0,0;0,0];

syms f1 f2 f3 s f4 f5
FN=[f1,f2,f3,f4,f5;0,0,0,0,0];

%GN=C*inv(s*eye(5)-A+B*F)*B;
%GN=ss(AN-BN*F,BN,CN,DN);
%t=0:0.1:10;
%x=sin(transpose(t));
%[Y,T]=lsim(GN,[zeros(101,1) x],t');

PO=det(s*eye(5)-A+B*FN);
PPO=coeffs(PO,s);
TT=(s-a1)*(s-a2)*(s-a3)*(s-a4)*(s-a5);
TTO=coeffs(TT,s);
% set f5=1
[f1,f2,f3,f4]=solve(PPO(1,1)==TTO(1,1),PPO(1,2)==TTO(1,2),PPO(1,3)==TTO(1,3),PPO(1,4)==TTO(1,4));
f1=double(f1);f2=double(f2);f3=double(f3);f4=double(f4);f5=1
FN=[f1,f2,f3,f4,f5;0,0,0,0,0];
G=ss(AN-BN*FN,BN,CN,DN);