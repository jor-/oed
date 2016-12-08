% Kostanten des Modells
R = 8.3144;
b = 2.857*10^(-10);
alpha =  1.0;
C_T = 38.45;
T1 = 2975;
Gref = 25500.;
c = [alpha b R C_T T1 Gref];

% Parameter des Modells
sigma0 = 30;
rho0 = (sigma0/alpha/Gref/b)^2;
m0 = 5.762100e-001;
U_0 = 1.728264e+008;
Omega0 = 1.183352e+001;
C = 3.374056e+005;
Qv = 1.1917*10^5;
p0 = [rho0; m0; U_0; Omega0; C; Qv];

% Erzeugung des Modell-Objektes
model = model_bergstroem(c);

% Mögliche Messbedingungen
T = [25 50 100 200 300 400];        % Temperatur
e = 0:0.01:0.8;                     % Dehnung
e_dot = [0.01 0.1 1 10];            % Dehnrate

% Exemplarischer Plot (T = 50,  e_dot = 1, e = 0:0.01:0.8)
n = length(e);
plot_res = zeros(n, 1);
for i = 1:n
    plot_res(i) = model.get_M(p0, [50, 1, e(i)]);
end
plot(e, plot_res);

% Solver Objekt erstellen
x = make_x(T, e, e_dot);            % alle Kombinationen von Temperatur, Dehnung und Dehnrate
m = length(x);
v = ones(m,1);                      % Varianz der einzelnen Messungen
sol = solver(model, p0, x, v);

% Erzielte Qualität beim Durchführen aller Messungen
w = ones(m,1);
best_quality = sol.get_quality(w);

% Erzielte Qualität beim Durchführen aller Messungen mit T >= 200
w(1 : length(w)/2) = 0;
quality_T_greater_equal_200 = sol.get_quality(w);

% Erzielte Qualität beim Durchführen aller Messungen mit T < 200
w(length(w)/2+1 : length(w)) = 0;
quality_T_less_200 = sol.get_quality(w);


% Funktion um x zu konstruieren aus allen möglichen Werten von T, e und e_dot
function x = make_x(T, e, e_dot)
    T_len = length(T);
    e_len = length(e);
    e_dot_len = length(e_dot);

    total_length = T_len * e_len * e_dot_len;

    x = ones(total_length, 3);

    for i = 1:T_len
        for j = 1:e_dot_len
            k = (i-1) * e_dot_len + j;
            indices = (k-1) * e_len + 1 : k * e_len;
            x(indices, 1) = T(i);
            x(indices, 2) = e_dot(j);
            x(indices, 3) = e;
        end
    end
end