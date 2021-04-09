%% ������� ������
% ���������� ���������
CO_sp1 = [6; 60] * pi / 180;
CO_sp2 = [-2; 50] * pi / 180;
% ������� �� ���������
az1 = -84.2467 * pi / 180;
az2 = 85.8152 * pi / 180;
% �������������� ���� ���������� 
psi = - 1.1;
lyam = 3;
CO = [psi; lyam] * pi / 180;

%% ��������������� �������
% ������������� ������� � �� ������� ����������� 
navig_fun = @(psi, lyam, psi_s, lyam_s)tan(psi_s) * cos(psi) / sin(lyam - lyam_s) - ...
            sin(psi) * cot(lyam - lyam_s);        
diff_psi  = @(psi, lyam, psi_s, lyam_s)-(sin(psi_s).*sin(psi) + cos(psi_s).*cos(psi).*cos(lyam - lyam_s))./(sin(lyam - lyam_s).*cos(psi_s));
diff_lyam = @(psi, lyam, psi_s, lyam_s)(cos(psi_s).*sin(psi) - cos(psi).*sin(psi_s).*cos(lyam - lyam_s))./(cos(psi_s).*(sin(lyam - lyam_s)).^2);
 
% ������ ������� �
B_ = @(psi, lyam, psi_s, lyam_s)[diff_psi(psi, lyam, psi_s, lyam_s), diff_lyam(psi, lyam, psi_s, lyam_s)];
% ������� B
B = @(CO, CO_sp1, CO_sp2)[B_(CO(1), CO(2), CO_sp1(1), CO_sp1(2));...
                          B_(CO(1), CO(2), CO_sp2(1), CO_sp2(2))];
% ������� ����������� ��������
delta_d = @(CO, CO_sp1, CO_sp2, az1, az2)[cot(az1) - navig_fun(CO(1), CO(2), CO_sp1(1), CO_sp1(2));...
                                          cot(az2) - navig_fun(CO(1), CO(2), CO_sp2(1), CO_sp2(2))];

% ������� ���������
delta_q = @(CO, CO_sp1, CO_sp2, az1, az2)B(CO, CO_sp1, CO_sp2) \ delta_d(CO, CO_sp1, CO_sp2, az1, az2);
 
% ����� ����������
new_q = @(CO, CO_sp1, CO_sp2, az1, az2)CO + delta_q(CO, CO_sp1, CO_sp2, az1, az2);

%% ��������� ���������
disp(['C�������� ���������� � ��������� ������:  [', num2str(psi), ', ', num2str(lyam), ']']);
for i = 1:3
    CO = new_q(CO, CO_sp1, CO_sp2, az1, az2);
    psi = CO(1) * 180 / pi;
    lyam = CO(2) * 180 / pi;
    disp(['    ', num2str(i), '-� �����������:                      [', num2str(psi, '%10.4f'), ', ', num2str(lyam, '%10.4f'), ']']);
end