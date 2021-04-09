%% входные данные
R_obj = 42000;
% предполагаемые координаты объекта
psi = 0; lyam = 71 * pi / 180;
% время
date = 14; month = 2; hour = 14; min = 11; t0 = 3;
% уточняем долготу
C = [0; 1];

%% постоянные 
w_sun = 72 / (73 * 1440);  % град/мин
w_earth = 0.25;  % град/мин
r_earth = 6300;  % км 

%% обобщенное условие видимости
% если объект выше спутников
obser_up = @(CO, CO_sp)((CO - CO_sp)' * CO_sp) / norm(CO - CO_sp) / norm(CO_sp) + ...
                         sqrt(1 - (r_earth / norm(CO_sp))^2);
% если объект ниже спутников
% obser_down = @(CO, CO_sp)((CO_sp - CO)' * CO) / norm(CO_sp - CO) / norm(CO) + ...
%                          sqrt(1 - (r_earth / norm(CO))^2);

%% параметры спутников
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
 
%% вычисление предполагаемых гринвич. прямоуг. коорд. объекта
x = @(psi, lyam)R_obj .* cos(psi) .* cos(lyam);
y = @(psi, lyam)R_obj .* cos(psi) .* sin(lyam);
z = @(psi, lyam)R_obj .* sin(psi);
CO = [x(psi, lyam); y(psi, lyam); z(psi, lyam)];

%% формирование строки матрицы B
% счислимое расстояние до n-го спутника
dist_sp = @(psi, lyam, sp)sqrt((x(psi, lyam) - sp(1)).^2 + ...
                               (y(psi, lyam) - sp(2)).^2 + ...
                               (z(psi, lyam) - sp(3)).^2);
% част. производные
ddist_psi  = @(psi, lyam, sp)R_obj ./ dist_sp(psi, lyam, sp) .* (sp(1) * cos(lyam) * sin(psi) ...
                - sp(2) * sin(lyam) * sin(psi) - sp(3) * cos(psi));
ddist_lyam = @(psi, lyam, sp)R_obj ./ dist_sp(psi, lyam, sp) .* cos(psi) .* ...
                 (sp(1) * sin(lyam) - sp(2) * cos(lyam));
% строка матрицы В
B_ = @(psi, lyam, sp)[ddist_psi(psi, lyam, sp), ddist_lyam(psi, lyam, sp)];
                     
%% проверка видимости спутников и подготовка к симплекс-методу
if month < 10
    mon = ['0', num2str(month)];
else
    mon = num2str(month);
end
disp(['Дата:              ', num2str(date), '.', mon]);
disp(['Время:             ', num2str(hour), ':', num2str(min)]);
disp('Наблюдаемые спутники: ');

H = zeros(2, 24);
num_V = zeros(1, 2);  % номера базисных спутников
% TRANZIT
for k = 1:length(tau_tranzit)
    if obser_up(CO, CO_tranzit(:, k)) > 0
        disp(['                   TRANZIT-', num2str(k)]);
        H(:, k) = B_(psi, lyam, CO_tranzit(:, k))';
        if num_V(1) == 0
            num_V(1) = k;
        elseif num_V(2) == 0
            num_V(2) = k;
        end
    end
end
% GPS
for k = 1:length(tau_gps)
    if obser_up(CO, CO_gps(:, k)) > 0
        disp(['                   GPS-', num2str(k)]);
        H(:, k + 6) = B_(psi, lyam, CO_gps(:, k))';
        if num_V(1) == 0
            num_V(1) = k + 6;
        elseif num_V(2) == 0
            num_V(2) = k + 6;
        end
    end
end
H(:, 13 : 24) = -H(:, 1 : 12);

%% симплекс-метод
while true
    % находим G
    G = [H(:, num_V(1)), H(:, num_V(2))];
    % находим вектор V
    V = G \ C;
    % проверяем на положительность
    if V(1) < 0
        num_V(1) = mod(num_V(1) + 11, 24) + 1;  % 24-ый становится 12-ым, а не нулевым
    end
    if V(2) < 0
        num_V(2) = mod(num_V(2) + 11, 24) + 1;
    end
    % пересчитываем G и V
    G = [H(:, num_V(1)), H(:, num_V(2))];
    V = G \ C;
    % вычисляем матрицу коэффициентов gamma
    gamma = G \ H;
    % находим сумму коэфф по столбцам матрицы gamma
    delta =  gamma(1, :) + gamma(2, :);
    % условие остановки
    cond = true;
    for k = 1 : 24
        if delta(k) > 1
            cond = false;
        end
    end
    if cond
        break
    end
    % пересчитываем номера строк, если условие остановки не выполнено
    opt_num_V = zeros(1, 2); min_frac = Inf;
    % находим дроби и ищем минимальную
    fracs = [V(1) ./ gamma(1, :); V(2) ./ gamma(2, :)];
    for i = 1 : 2
        for k = 1 : 24
            if (delta(k) > 1) && (gamma(i, k) > 0)
                if fracs(i,k) < min_frac
                    min_frac = fracs(i,k);
                    opt_num_V = [i, k];
                end
            end
        end
    end
    num_V(opt_num_V(1)) = opt_num_V(2);
end

%% итого
disp('Оптимальная пара для уточнения долготы:')
opt = mod(num_V - [1, 1], 12) + [1, 1];
if opt(1) <= 6
    disp(['                   TRANZIT-', num2str(opt(1))]);
else
    disp(['                   GPS-', num2str(opt(1) - 6)]);
end
if opt(2) <= 6
    disp(['                   TRANZIT-', num2str(opt(2))]);
else
    disp(['                   GPS-', num2str(opt(2) - 6)]);
end