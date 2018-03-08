
clear all;
clear u X T t tt
A=[-1.3,0.98,0,-0.165,-0.248;42.81,-0.785,0,-17.3,-1.58;1.25,0.007,0,0.165,0.248;0,0,0,-18,0;0,0,0,0,-18];
B=[0,0;0,0;0,0;18,0;0,18];
C=[0,1,0,0,0;46.5,-0.256,0,-4.25,4.15;0,0,1,0,0];
D=[0,0;0,0;0,0];

% step response
figure(1)
sys=ss(A,B,C,D);
t=0:0.001:10;
tt=transpose(t);

u=[ones(size(tt)) ones(size(tt))];
[X,T]=lsim(sys,u,t);
XX=X(:,1);
plot(T,X)

% impulse disturbance
figure(2)
BB=[B,[1;0;0;0;0]];
D=[D,[0;0;0]]
sys2=ss(A,BB,C,D);
dis=[1;zeros(size(transpose(0:0.001:9.999)))]
uu=[zeros(size(tt)) zeros(size(tt)) dis];
[Y,TT]=lsim(sys2,uu,t)
plot(TT,Y)
