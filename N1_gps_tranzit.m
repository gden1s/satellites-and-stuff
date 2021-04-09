%% ���������� 
w_sun = 72 / (73 * 1440);  % ����/���
w_earth = 0.25;  % ����/���
r_earth = 6300;  % �� 

%% �������� ���������� � �������. ������. � �������. �� 
r0 = r_earth;  % ��
psi0  = 80;  % �������� ������
lyam0 = 20;  % ��������� �������

% ������� � �������. �������.
x0 = r0 * cos(psi0 * pi / 180) * cos(lyam0 * pi / 180);
y0 = r0 * cos(psi0 * pi / 180) * sin(lyam0 * pi / 180);
z0 = r0 * sin(psi0 * pi / 180);
CO = [x0; y0; z0];

%% ����� 
month = 2;
date = 11;
hour = 11;
min = 30;
t0 = 3; % ���
% ������ ����� � �������
b = [0 31 28 31 30 31 30 31 31 30 31 30];
eye_ = [1 1 1 1 1 1 1 1 1 1 1 1]; 
bk = zeros(12, 1); bk(1 : month) = eye_(1 : month);
t = (24 * mod(b * bk + date - 81, 365) + hour - t0) * 60 + min; 

%% ��������� ��������� � �� ���������� ����������
% ������� �������� �� �������. � �������. �� (�.� A^(-1))
gam = (w_sun + w_earth) * pi / 180 * t; % ���
A_inv = [ cos(gam), sin(gam), 0;...
         -sin(gam), cos(gam), 0;...
          0,        0,        1];

% ������� ��������� TRANZIT � ������������ ��
i_tranzit = 90 * pi / 180; % ���
R_tranzit = 7500;
w_tranzit = 3 * pi / 180; % ���/���
tau_tranzit = [40 10 100 80 50 110];
OMEGA_tranzit = [30 60 120 210 270 330] * pi / 180; % ���;

% ������� ��������� TRANZIT � �������. ��
u_tranzit = w_tranzit * (t - tau_tranzit); % ���
CO_tranzit = A_inv * R_tranzit * ...
    [cos(OMEGA_tranzit) .* cos(u_tranzit) - sin(OMEGA_tranzit) .* sin(u_tranzit) * cos(i_tranzit);...
     sin(OMEGA_tranzit) .* cos(u_tranzit) + cos(OMEGA_tranzit) .* sin(u_tranzit) * cos(i_tranzit);...
     sin(u_tranzit) * sin(i_tranzit)];
% ����� ����������
x_tranzit = CO_tranzit(1, :);
y_tranzit = CO_tranzit(2, :);
z_tranzit = CO_tranzit(3, :);


% ������� ��������� GPS � ������������ ��
i_gps = 60 * pi / 180; % ���
R_gps = 15000;
w_gps = 2 * pi / 180; % ���/���
tau_gps = [30 100 70 0 150 120];
OMEGA_gps = [0 45 90 135 210 300] * pi / 180; % ���

% ������� ��������� GPS � �������. ��
u_gps = w_gps * (t - tau_gps); % ���
CO_gps = A_inv * R_gps * ...
    [cos(OMEGA_gps) .* cos(u_gps) - sin(OMEGA_gps) .* sin(u_gps) * cos(i_gps);...
     sin(OMEGA_gps) .* cos(u_gps) + cos(OMEGA_gps) .* sin(u_gps) * cos(i_gps);...
     sin(u_gps) * sin(i_gps)];
% ����� ����������
x_gps = CO_gps(1, :);
y_gps = CO_gps(2, :);
z_gps = CO_gps(3, :);

%% �������� ��������� ���������
% TRANZIT
for i = 1:length(tau_tranzit)
    R_ = CO_tranzit(:, i) - CO;
    if CO' * R_ > 0
        disp(['TRANZIT: ', num2str(i)])
    end
end
% GPS
for i = 1:length(tau_gps)
    R_ = CO_gps(:, i) - CO;
    if CO' * R_ > 0
        disp(['GPS: ', num2str(i)])
    end
end