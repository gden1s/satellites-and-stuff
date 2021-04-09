%% ������� ������
% ������� ������� GPS
num_gps = 2;
% ������� ����� (��� �������)
T = fix(clock); month = T(2); date = T(3); hour = T(4); min = T(5); t0 = 3;
%% ��������� ���������
% omega = 0;
% ������� ��������� TRANZIT � ������������ ��
i_tranzit = 90 * pi / 180;  % ���
R_tranzit = 7500;
w_tranzit = 3 * pi / 180;  % ���/���
tau_tranzit = [40 10 100 80 50 110];
OMEGA_tranzit = [30 60 120 210 270 330] * pi / 180;  % ���;
% ������� ��������� GPS � ������������ ��
i_gps = 60 * pi / 180;  % ���
R_gps = 15000;
w_gps = 2 * pi / 180;  % ���/���
tau_gps = [30 100 70 0 150 120];
OMEGA_gps = [0 45 90 135 210 300] * pi / 180;  % ���

%% ���������� 
w_sun = 72 / (73 * 1440);  % ����/���
w_earth = 0.25;  % ����/���
r_earth = 6300;  % �� 

%% ������ ����� � �������
b = [0 31 28 31 30 31 30 31 31 30 31 30];
eye_ = [1 1 1 1 1 1 1 1 1 1 1 1]; 
bk = zeros(12, 1); bk(1 : month) = eye_(1 : month);
t = (24 * mod(b * bk + date - 81, 365) + hour - t0) * 60 + min; 

%% ���������� ��������� ���������
% ������� �������� �� �������. � �������. �� (�.� A^(-1))
gam = (w_sun + w_earth) * pi / 180 * t;  % ���
A_inv = [ cos(gam), sin(gam), 0;...
         -sin(gam), cos(gam), 0;...
          0,        0,        1];

% ������� ��������� TRANZIT � �������. ��
u_tranzit = w_tranzit * (t - tau_tranzit);  % ���
CO_tranzit = A_inv * R_tranzit * ...
    [cos(OMEGA_tranzit) .* cos(u_tranzit) - sin(OMEGA_tranzit) .* sin(u_tranzit) * cos(i_tranzit);...
     sin(OMEGA_tranzit) .* cos(u_tranzit) + cos(OMEGA_tranzit) .* sin(u_tranzit) * cos(i_tranzit);...
     sin(u_tranzit) * sin(i_tranzit)];

% ������� ��������� GPS � �������. ��
u_gps = w_gps * (t - tau_gps);  % ���
CO_gps = A_inv * R_gps * ...
    [cos(OMEGA_gps) .* cos(u_gps) - sin(OMEGA_gps) .* sin(u_gps) * cos(i_gps);...
     sin(OMEGA_gps) .* cos(u_gps) + cos(OMEGA_gps) .* sin(u_gps) * cos(i_gps);...
     sin(u_gps) * sin(i_gps)];

%% ���������� ������� ���������
% ���� ����������� ���� ������������ (��� �����!)
obser_up = @(CO, CO_sp)((CO - CO_sp)' * CO_sp) / norm(CO - CO_sp) / norm(CO_sp) + ...
                         sqrt(1 - (r_earth / norm(CO_sp))^2);
% % ���� ����������� ���� ������������
% obser_down = @(CO, CO_sp)((CO_sp - CO)' * CO) / norm(CO_sp - CO) / norm(CO) + ...
%                          sqrt(1 - (r_earth / norm(CO))^2);

%% �������� ��������� ���������
if month < 10
    mon = ['0', num2str(month)];
else
    mon = num2str(month);
end
disp(['����:              ', num2str(date), '.', mon]);
disp(['�����:             ', num2str(hour), ':', num2str(min)]);
disp(['�����������:       GPS-', num2str(num_gps)]);  % ��������� � GPS
disp('����������� �������� ������� TRANZIT: ');
% TRANZIT
CO = CO_gps(:, num_gps);
for k = 1:length(tau_tranzit)
    if obser_up(CO, CO_tranzit(:, k)) > 0
        disp(['    TRANZIT-', num2str(k)]);
    end
end