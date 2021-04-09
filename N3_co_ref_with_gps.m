%% координаты и время
% истинные координаты в гринвич. сферич. СК
psiTrue  = 0.7;  % северной широты
lyamTrue = 0.55;  % восточной долготы
% предполагаемые координаты в гринвич. сферич. СК
psi = 0.5;  % с.ш
lyam = 1;   % в.д
% выбранное время 
month = 2;
date = 17;
hour = 13;
min = 50;
t0 = 3;  % мск
%  % текущее время (при желании)
% T = fix(clock); month = T(2); date = T(3); hour = T(4); min = T(5); t0 = 3;
%% параметры спутников
omega = 0;
% система спутников TRANZIT в оскулирующей СК
i_tranzit = 90 * pi / 180;  % рад
R_tranzit = 7500;
w_tranzit = 3 * pi / 180;  % рад/мин
tau_tranzit = [40 10 100 80 50 110];
OMEGA_tranzit = [30 60 120 210 270 330] * pi / 180;  % рад;
% система спутников GPS в оскулирующей СК
i_gps = 60 * pi / 180;  % рад
R_gps = 15000;
w_gps = 2 * pi / 180;  % рад/мин
tau_gps = [30 100 70 0 150 120];
OMEGA_gps = [0 45 90 135 210 300] * pi / 180;  % рад

%% постоянные 
w_sun = 72 / (73 * 1440);  % град/мин
w_earth = 0.25;  % град/мин
r_earth = 6300;  % км 

%% нахождение истинных гринвич. прямоуг. коорд. объекта
CO = zeros(3, 1);
CO(1) = r_earth * cos(psiTrue * pi / 180) * cos(lyamTrue * pi / 180);
CO(2) = r_earth * cos(psiTrue * pi / 180) * sin(lyamTrue * pi / 180);
CO(3) = r_earth * sin(psiTrue * pi / 180);

%% точное время в минутах
b = [0 31 28 31 30 31 30 31 31 30 31 30];
eye_ = [1 1 1 1 1 1 1 1 1 1 1 1]; 
bk = zeros(12, 1); bk(1 : month) = eye_(1 : month);
t = (24 * mod(b * bk + date - 81, 365) + hour - t0) * 60 + min; 

%% вычисление координат спутников
% матрица перехода от экватор. к гринвич. СК (т.е A^(-1))
gam = (w_sun + w_earth) * pi / 180 * t;  % рад
A_inv = [ cos(gam), sin(gam), 0;...
         -sin(gam), cos(gam), 0;...
          0,        0,        1];

% система спутников TRANZIT в гринвич. СК
u_tranzit = w_tranzit * (t - tau_tranzit);  % рад
CO_tranzit = A_inv * R_tranzit * ...
    [cos(OMEGA_tranzit) .* cos(u_tranzit) - sin(OMEGA_tranzit) .* sin(u_tranzit) * cos(i_tranzit);...
     sin(OMEGA_tranzit) .* cos(u_tranzit) + cos(OMEGA_tranzit) .* sin(u_tranzit) * cos(i_tranzit);...
     sin(u_tranzit) * sin(i_tranzit)];

% система спутников GPS в гринвич. СК
u_gps = w_gps * (t - tau_gps);  % рад
CO_gps = A_inv * R_gps * ...
    [cos(OMEGA_gps) .* cos(u_gps) - sin(OMEGA_gps) .* sin(u_gps) * cos(i_gps);...
     sin(OMEGA_gps) .* cos(u_gps) + cos(OMEGA_gps) .* sin(u_gps) * cos(i_gps);...
     sin(u_gps) * sin(i_gps)];

%% проверка видимости спутников и выбор двух наиболее наблюдаемых
if month < 10
    mon = ['0', num2str(month)];
else
    mon = num2str(month);
end
disp(['Дата:                     ', num2str(date), '.', mon]);
disp(['Время:                    ', num2str(hour), ':', num2str(min)]);
disp(['Истинные коорд. об.:      [', num2str(psiTrue), ', ', num2str(lyamTrue), ']']);
disp(['Счисленные коорд. об.:    [', num2str(psi), ', ', num2str(lyam), ']']);
disp('Наблюдаемые спутники, угол наблюдения: ');
min_angle1 = pi/2; min_angle2 = pi/2;
num1_str = ''; num1 = 0; 
% TRANZIT
for k = 1:length(tau_tranzit)
    R_ = CO_tranzit(:, k) - CO;
    if CO' * R_ > 0
        disp(['                          TRANZIT-', num2str(k), ': ', num2str(acos(CO' * R_ / norm(CO) / norm(R_)) * 180 / pi)])
        ang = acos(CO' * R_ / (6300 * norm(R_)));
        if ang < min_angle1
            min_angle2 = min_angle1;
            num2_str = num1_str;
            num2 = num1; 
            min_angle1 = ang;
            num1_str = ['TRANZIT-', num2str(k)];
            num1 = [1, k]; 
        elseif ang < min_angle2
            min_angle2 = ang;
            num2_str = ['TRANZIT-', num2str(k)];
            num2 = [1, k]; 
        end
    end
end
% GPS
for k = 1:length(tau_gps)
    R_ = CO_gps(:, k) - CO;
    if CO' * R_ > 0
        disp(['                          GPS-', num2str(k), ': ', num2str(acos(CO' * R_ / norm(CO) / norm(R_)) * 180 / pi)])
        ang = acos(CO' * R_ / (6300 * norm(R_)));
        if ang < min_angle1
            min_angle2 = min_angle1;
            num2_str = num1_str;
            num2 = num1; 
            min_angle1 = ang;
            num1_str = ['GPS-', num2str(k)];
            num1 = [1, k]; 
        elseif ang < min_angle2
            min_angle2 = ang;
            num2_str = ['GPS-', num2str(k)];
            num2 = [1, k]; 
        end
    end
end
if num2 ~= 0
    disp(['Используются ', num1_str, ' и ', num2_str, '.']);
else
    error('Наблюдается только один спутник! Вычисления не произведены.');
end

%% подготовка к уточнению счисленных координат
psi = psi * pi / 180; 
lyam = lyam * pi / 180;
% гринвич. прямоуг. координаты объекта
x = @(psi, lyam)r_earth .* cos(psi) .* cos(lyam);
y = @(psi, lyam)r_earth .* cos(psi) .* sin(lyam);
z = @(psi, lyam)r_earth .* sin(psi);

% счислимое расстояние до n-го спутника
dist_sp = @(psi, lyam, sp)sqrt((x(psi, lyam) - sp(1)).^2 + ...
                               (y(psi, lyam) - sp(2)).^2 + ...
                               (z(psi, lyam) - sp(3)).^2);

% част. производные
ddist_psi  = @(psi, lyam, sp)r_earth ./ dist_sp(psi, lyam, sp) .* (sp(1) * cos(lyam) * sin(psi) ...
                - sp(2) * sin(lyam) * sin(psi) - sp(3) * cos(psi));
ddist_lyam = @(psi, lyam, sp)r_earth ./ dist_sp(psi, lyam, sp) .* cos(psi) .* ...
                 (sp(1) * sin(lyam) - sp(2) * cos(lyam));
 
% строка матрицы В
B_ = @(psi, lyam, sp)[ddist_psi(psi, lyam, sp), ddist_lyam(psi, lyam, sp)];
% невязка расстояния до спутника
delta_d = @(psi, lyam, sp, dist)dist - dist_sp(psi, lyam, sp);

% пересчет координат
i = zeros(1,2); R = zeros(1,2); tau = zeros(1,2); OMEGA = zeros(1,2); w = zeros(1,2);
if num1(1) == 1
    i(1) = i_tranzit;
    R(1) = R_tranzit;
    tau(1) = tau_tranzit(num1(2));
    OMEGA(1) = OMEGA_tranzit(num1(2));
    w(1) = w_tranzit;
elseif num1(1) == 2
    i(1) = i_gps;
    R(1) = R_gps;
    tau(1) = tau_gps(num1(2));
    OMEGA(1) = OMEGA_gps(num1(2));
    w(1) = w_gps;
end
if num2(1) == 1
    i(2) = i_tranzit;
    R(2) = R_tranzit;
    tau(2) = tau_tranzit(num2(2));
    OMEGA(2) = OMEGA_tranzit(num2(2));
    w(2) = w_tranzit;
elseif num2(1) == 2
    i(2) = i_gps;
    R(2) = R_gps;
    tau(2) = tau_gps(num2(2));
    OMEGA(2) = OMEGA_gps(num2(2));
    w(2) = w_gps;
end
CO_sp1 = zeros(3, 10); CO_sp2 = zeros(3, 10);
dist1 = zeros(1, 10); dist2 = zeros(1, 10); 
B = zeros(20, 2);
dd = zeros(20, 1);

%% уточнение счисленных координат 
K_inv = [8, -4; -4, 8];  % K = [1/6 1/12; 1/12 1/6]
n_iter = 3;  % количество интераций
B_cell = cell(1, n_iter);
d_cell = cell(1, n_iter);
for k = 1 : n_iter
    B_cell{k} = zeros(2, 2);
    d_cell{k} = zeros(2, 1);
end
k = 1;
err = 1;
while err > 0.01
    disp(['Количество пар измерений: ', num2str(k)]);
    gam = (w_sun + w_earth) * pi / 180 * t;
    A_inv = [ cos(gam), sin(gam), 0;...
             -sin(gam), cos(gam), 0;...
              0,        0,        1];
    % общая ошибка измерения 
    err_0 = (rand(1) - 0.5);
    % ошибки измерения для спутников
    err_1 = (rand(1) - 0.5); err_2 = (rand(1) - 0.5);

    % истинные координаты спутников и расстояния (c ошибкой) до них
    u = w .* (t - tau);
    CO_sp1(:, k) = A_inv * R(1) * ...
        [cos(OMEGA(1)) * cos(u(1)) - sin(OMEGA(1)) * sin(u(1)) * cos(i(1));...
         sin(OMEGA(1)) * cos(u(1)) + cos(OMEGA(1)) * sin(u(1)) * cos(i(1));...
         sin(u(1)) * sin(i(1))];
    dist1(k) = norm(CO_sp1(:, k) - CO) + err_0 + err_1;
    CO_sp2(:, k) = A_inv * R(2) * ...
        [cos(OMEGA(2)) * cos(u(2)) - sin(OMEGA(2)) * sin(u(2)) * cos(i(2));...
         sin(OMEGA(2)) * cos(u(2)) + cos(OMEGA(2)) * sin(u(2)) * cos(i(2));...
         sin(u(2)) * sin(i(2))];
    dist2(k) = norm(CO_sp2(:, k) - CO) + err_0 + err_2;
     % начальная точка - наши первые счисленные координаты 
    psi_ = psi; lyam_ = lyam;
    for j = 1 : n_iter
        B_k = [B_(psi_, lyam_, CO_sp1(:, k));...
               B_(psi_, lyam_, CO_sp2(:, k))];
        B_cell{j} = B_cell{j} + B_k' * K_inv * B_k;
        d_k = [delta_d(psi_, lyam_, CO_sp1(:, k), dist1(k)); ...
               delta_d(psi_, lyam_, CO_sp2(:, k), dist2(k))];
        d_cell{j} = d_cell{j} + B_k' * K_inv * d_k;

        % невязка координат объекта
        delta_q = B_cell{j} \ d_cell{j}; 
        % новые координаты
        q = [psi_; lyam_] + delta_q;
        psi_ = q(1); lyam_ = q(2);
        disp(['    ', num2str(j), '-е приближение:      [', num2str(psi_ * 180 / pi), ', ', num2str(lyam_ * 180 / pi), ']']);
    end
    % новое время, новые вычисления координат спутников и расстояний до них
    t = t + 0.1;
    k = k + 1;
    err = norm([psiTrue; lyamTrue] - q * 180 / pi);
    disp(['    Текущая ошибка:       ', num2str(err)]);
end