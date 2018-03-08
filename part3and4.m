clear all;
syms p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 q1 q2 q3 q4 q5 q6 q7 q8 q9 q10

% A,B,C,D
A=[-1.3,0.98,0,-0.165,-0.248;42.81,-0.785,0,-17.3,-1.58;1.25,0.007,0,0.165,0.248;0,0,0,-18,0;0,0,0,0,-18];
B=[0,0;0,0;0,0;18,0;0,18];
C=[0,1,0,0,0;46.5,-0.256,0,-4.25,4.15;0,0,1,0,0];
CC=[0,1,0,0,0;0,0,1,0,0];
D=[0,0;0,0;0,0];
DD=[0,0;0,0]

% assign eigenvector and Q
P1=[1;p1;0;p2;p3];
P2=[0;1;p4;p5;p6];
P3=[p7;0;1;p8;p9];
P4=[p10;p11;p12;1;0];
P5=[p13;p14;p15;0;1];
Q1=[q1;q2];
Q2=[q3;q4];
Q3=[q5;q6];
Q4=[q7;q8];
Q5=[q9;q10];

% assign eigenvalues
a1=-4;a2=-5;a3=-5;a4=-19;a5=-19.5;

% solve P AND Q
[p1,p2,p3,q1,q2]=solve(P1==inv(a1*eye(5)-A)*B*Q1);
[p4,p5,p6,q3,q4]=solve(P2==inv(a2*eye(5)-A)*B*Q2);
[p7,p8,p9,q5,q6]=solve(P3==inv(a3*eye(5)-A)*B*Q3);
[p10,p11,p12,q7,q8]=solve(P4==inv(a4*eye(5)-A)*B*Q4);
[p13,p14,p15,q9,q10]=solve(P5==inv(a5*eye(5)-A)*B*Q5);
p1=double(p1);p2=double(p2);p3=double(p3);p4=double(p4);p5=double(p5);p6=double(p6);p7=double(p7);
p8=double(p8);p9=double(p9);p10=double(p10);p11=double(p11);p12=double(p12);p13=double(p13);
p14=double(p14);p15=double(p15);q1=double(q1);q2=double(q2);q3=double(q3);q4=double(q4);q5=double(q5);q6=double(q6);q7=double(q7);
q8=double(q8);q9=double(q9);q10=double(q10);

P=[1,0,p7,p10,p13;p1,1,0,p11,p14;0,p4,1,p12,p15;p2,p5,p8,1,0;p3,p6,p9,0,1];
Q=[q1,q3,q5,q7,q9;q2,q4,q6,q8,q10]

F=-Q*inv(P);


%% part 4 SVD

AA=A-B*F;
G=ss(AA,B,C,D);
GG=tf(G);
y=-C*inv(AA-B*F)*B;
y2=-CC*inv(AA-B*F)*B;

% scaling matrix before r to make steady state error equals to zero
K=-inv(CC*inv(A-B*F)*B);
G=ss(AA,B*K,C,D);
GG=tf(G)


sigma(GG); % obtain maximum sigular value for G is orcuured when the frequency is around 12
w=12;
GW=evalfr(GG,i*w);
[U,S,V]=svd(GW);


sx=S(1,1);
sn=S(2,2);
ux=U(:,1);
un=U(:,2);
vx=V(:,1);
vn=V(:,2);

[th_vmx1, amp_vmx1] = cart2pol(real(vx(1)),imag(vx(1)));
[th_vmx2, amp_vmx2] = cart2pol(real(vx(2)),imag(vx(2)));
[th_vmn1, amp_vmn1] = cart2pol(real(vn(1)),imag(vn(1)));
[th_vmn2, amp_vmn2] = cart2pol(real(vn(2)),imag(vn(2)));

[th_umx1, amp_umx1] = cart2pol(real(ux(1)),imag(ux(1)));
[th_umx2, amp_umx2] = cart2pol(real(ux(2)),imag(ux(2)));
[th_umx3, amp_umx3] = cart2pol(real(ux(3)),imag(ux(3)));
[th_umn1, amp_umn1] = cart2pol(real(un(1)),imag(un(1)));
[th_umn2, amp_umn2] = cart2pol(real(un(2)),imag(un(2)));
[th_umn3, amp_umn3] = cart2pol(real(un(3)),imag(un(3)));

% MAX
t=0:0.01:10;

Umx1=amp_vmx1*sin(w*t+th_vmx1);
Umx2=amp_vmx2*sin(w*t+th_vmx2);

Ymx1=sx*abs(amp_umx1)*sin(w*t+th_umx1);
Ymx2=sx*abs(amp_umx2)*sin(w*t+th_umx2);
Ymx3=sx*abs(amp_umx3)*sin(w*t+th_umx3);
figure(1)
subplot(2,1,1)
plot(t, Umx1,'b'); hold on
plot(t, Umx2,'r'); hold off
xlabel('time (sec)'); ylabel('u_{max}')
subplot(2,1,2)
plot(t, Ymx1,'b'); hold on
plot(t, Ymx2,'r'); hold on
plot(t, Ymx3);hold off
xlabel('time (sec)'); ylabel('y_{max}')

%min
Umn1=abs(amp_vmn1)*sin(w*t+th_vmn1);
Umn2=abs(amp_vmn2)*sin(w*t+th_vmn2);

Ymn1=sn*abs(amp_umn1)*sin(w*t+th_umn1);
Ymn2=sn*abs(amp_umn2)*sin(w*t+th_umn2);
Ymn3=sn*abs(amp_umn3)*sin(w*t+th_umn3);
figure(2)
subplot(2,1,1)
plot(t, Umn1,'b'); hold on
plot(t, Umn2,'r'); hold off
xlabel('time (sec)'); ylabel('u_{min}')
subplot(2,1,2)
plot(t, Ymn1,'b'); hold on
plot(t, Ymn2,'r'); 
plot(t, Ymn3,'r');hold off
xlabel('time (sec)'); ylabel('y_{min}')
