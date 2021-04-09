%% входные данные
R = 6300;
% координаты спутников
sp1 = [8757.75, -2183.55,  4305.11];
sp2 = [9061.25,  2598.27, -3338.07];
% расстояние до спутников
d1 = 5381.18;
d2 = 5072.63;
% предполагаемые наши координаты 
psi = 0.5 * pi / 180;
lyam = 1 * pi/ 180;

%% вспомогательные функции
% гринвич. прямоуг. координаты объекта
x = @(psi, lyam)R .* cos(psi) .* cos(lyam);
y = @(psi, lyam)R .* cos(psi) .* sin(lyam);
z = @(psi, lyam)R .* sin(psi);
 
% счислимое расстояние до n-го спутника
dist_sp = @(psi, lyam, sp)sqrt((x(psi, lyam) - sp(1)).^2 + ...
                               (y(psi, lyam) - sp(2)).^2 + ...
                               (z(psi, lyam) - sp(3)).^2);

% част. производные
ddist_psi  = @(psi, lyam, sp)R ./ dist_sp(psi, lyam, sp) .* (sp(1) * cos(lyam) * sin(psi) ...
                - sp(2) * sin(lyam) * sin(psi) - sp(3) * cos(psi));
ddist_lyam = @(psi, lyam, sp)R ./ dist_sp(psi, lyam, sp) .* cos(psi) .* ...
                 (sp(1) * sin(lyam) - sp(2) * cos(lyam));
 
% матрица В
B = @(psi, lyam, sp1, sp2)[ddist_psi(psi, lyam, sp1) ddist_lyam(psi, lyam, sp1); ...
                           ddist_psi(psi, lyam, sp2) ddist_lyam(psi, lyam, sp2)]; 

% уточнение расстояния до спутников
delta_d = @(psi, lyam, sp1, sp2, d1, d2)[d1 - dist_sp(psi, lyam, sp1); ...
                                         d2 - dist_sp(psi, lyam, sp2)];
% уточнение координат
delta_q = @(psi, lyam, sp1, sp2, d1, d2)B(psi, lyam, sp1, sp2) \ delta_d(psi, lyam, sp1, sp2, d1, d2);
 
% новые координаты
new_q = @(psi, lyam, sp1, sp2, d1, d2)[psi; lyam] + delta_q(psi, lyam, sp1, sp2, d1, d2);

%% уточнение координат
for i = 1:3
    q = new_q(psi, lyam, sp1, sp2, d1, d2);
    psi = q(1);
    lyam = q(2);
    disp([num2str(i), '-е приближение:']);
    disp(q * 180 / pi);
end