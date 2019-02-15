clc
clear all
 
M = input('Enter Mass : ');
K = input('Enter spring stiffener : ');
e = input('Enter damping coefficient: ');
g = input('Enter ground acceleration : ');
Ts = input('Enter time step: ');
td = input('Enter duration of loading: ');
N = input('Enter number of points: ');
Code = input('Type of loading? (0 for base accel, 1 for applied force): ');
 
theta = 1.42;
 
T = zeros(1);
T(1) = 0;
for i = 1:N - 1
   T(i + 1) = i * Ts; 
end
 
x = zeros(1);
x(1) = 0;
for i = 1:N - 1
   x(i + 1) = i; 
end
 
F = xlsread('Duhamel','input','B2:B102');
j = (td / Ts) + 2;
F(j:N) = 0;
F = F(1:N);
 
if Code == 0
   F = -g * M * F; 
end
 
 
Wn = sqrt(K / M);
Period = 2 * pi / Wn;
C = 2 * e * sqrt(K * M);
Wd = Wn * (1 - e ^ 2);
 
V = zeros(0);
an = zeros(0);
U = zeros(0);
Ps = zeros(1);
Pd = zeros(1);
P = zeros(1);
U(1) = 0; V(1) = 0;
an(1) = (F(1) - C * V(1) - K * U(1)) / M;
kh = K + 3 * C / (theta * Ts) + 6 * M /(theta * Ts) ^ 2;
a = 6 * M / (theta * Ts) + 3 * C;
b = 3 * M + theta * Ts * C / 2;
 
for i = 2:N
   ww = (F(i) - F(i-1)) * theta + a * V(i-1) + b * an(i-1);
   xx = ww / kh;
   zz = (6 * xx / ((theta * Ts) ^ 2) - 6 * V(i-1) / (theta * Ts) - 3 * an(i-1)) / theta;
   yy = Ts * an(i-1) + Ts * zz / 2;
   V(i) = (V(i-1) + yy) / 1;
   an(i) = an(i-1) + zz;
   vv= Ts * V(i-1) + Ts *Ts * (3 * an(i-1) + zz) / 6;
   U(i) = U(i-1) + vv;
   Ps(i) = K * U(i);
   Pd(i) = C * V(i);
   P(i) = sqrt(Ps(i) ^ 2 + Pd(i)^2);   
end
 
   Y_max = max(abs(U));
   V_max = max(abs(V));
   a_max = max(abs(an));
   P_max = max(P);
   
format short g
disp('         Time        Force        Displ.    Velocity        Accel      Sup.Reac')
disp('-------------------------------------------------------------------------------')
out = [T' F U' V' an' P'];
disp(out)
 
disp(' ')
fprintf('Maximum Displacement          = %g\n', Y_max);
fprintf('Maximum Velocity              = %g\n', V_max);
fprintf('Maximum Acceleration          = %g\n', a_max);
fprintf('Maximum Support Reaction       = %g\n', P_max);
disp('------------------------------------------------------------------------')
 
figure(1);
    subplot(2, 2, 1);
    plot(T, U, 'b-', 'linewidth', 1.5);
    title('Response Displacement');
    xlabel('Time [sec]');
    ylabel('Displacement');
    grid on;
 
 
    subplot(2, 2, 2);
    plot(T, V, 'g-', 'linewidth', 1.5);
    title('Risponse Velocity');
    xlabel('Time [sec]');
    ylabel('Velocity');
    grid on;
    
    subplot(2, 2, 3);
    plot(T, an, 'r-', 'linewidth', 1.5);
    title('Response Acceleration');
    xlabel('Time [sec]');
    ylabel('Acceleration');
    grid on;
    
    subplot(2, 2, 4);
    plot(T, -Ps, 'k-', 'linewidth', 1.5);
    title('Response Reaction');
    xlabel('Time [sec]');
    ylabel('Support Reaction');
    grid on;
