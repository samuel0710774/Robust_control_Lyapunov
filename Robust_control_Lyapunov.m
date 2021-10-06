%%
clc; clear all;
%%
A = [0 0  1 0;
     0 0  0 1;
     0 -1 0 0;
     1 0  1 0];
b = [1;
     0;
     0;
     1];
E = [1 0;
     1 0;
     0 1;
     0 1];
   
x0 = [5;
      0;
     -6;
      0];
  
rank([A A*b A^2*b A^3*b]); %controllable
eigen_value = eig(A); %not stable
     
%pole = [-1+j -1-j -2 -5]; % need to explain
pole = [-1 -5 -10+5j -10-5j];
K = place(A,b,pole);

G = null(b','r');
t = [0:0.001:10];

%% Part(D)
d1 = 0.1*cos(2*t);
d2 = -0.2*sin(3*t);

d = [d1; d2];
 

dmm = inv(G'*G)*G'*E*d;
dm = inv(b'*b)*b'*E*d;

CH = inv([b G])*E;
c = CH(1,:);
D = CH(2:4,:);

alpha = max(abs(dm));

%% Part(E)

%(1)
Ac = A-b*K
Q = eye(4);
P = lyap(Ac',Q);

%(2)
H = P*G*inv(G'*G)*G'*E;

%(3)
H_norm = norm(H);

%(4)
R = 2*norm(H)*max(norm(d));






