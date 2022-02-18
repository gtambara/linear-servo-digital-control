%%% questao 1

% declarações de constantes e vetores

T = 0.5;
k = 60;
xi = 0.5;
wn = 1;

y = zeros(1, k);
u = ones(1, k);
t = 1:1:k;

% iterador para equações de diferenças
for i=3:1:k
    y(i) = ((wn*T)^2)*u(i-2) - (2*xi*wn*T-2)*y(i-1) - ((wn*T)^2 - 2*xi*wn*T + 1)*y(i-2);
    %y(i) = ((wn*T)^2/((wn*T)^2 + 2*xi*wn*T + 1))*u(i)  + ((2*xi*wn*T+2)/((wn*T)^2 + 2*xi*wn*T + 1))*y(i-1) - (((wn*T)^2 + 2*xi*wn*T + 1)^(-1))*y(i-2);
end

stairs(t, y)

%%% questao 2

% declarações de constantes e vetores

k2 = [27.825; 59.93];
num = [0.003679 0.002642];
den = [1 -1.3679 0.3679];

% o -1 significa dominio Z(discretizado)
G = tf(num, den, -1);

%rlocus(H, [27.825, 59.93])

H = feedback(G, 1);

rlocusplot(H, k2);
G
H%%% questao 1

% declarações de constantes e vetores

T = 0.5;
k = 60;
xi = 0.5;
wn = 1;

y = zeros(1, k);
u = ones(1, k);
t = 1:1:k;

% iterador para equações de diferenças
for i=3:1:k
    y(i) = ((wn*T)^2)*u(i-2) - (2*xi*wn*T-2)*y(i-1) - ((wn*T)^2 - 2*xi*wn*T + 1)*y(i-2);
    %y(i) = ((wn*T)^2/((wn*T)^2 + 2*xi*wn*T + 1))*u(i)  + ((2*xi*wn*T+2)/((wn*T)^2 + 2*xi*wn*T + 1))*y(i-1) - (((wn*T)^2 + 2*xi*wn*T + 1)^(-1))*y(i-2);
end

stairs(t, y)

%%% questao 2

% declarações de constantes e vetores

k2 = [27.825; 59.93];
num = [0.003679 0.002642];
den = [1 -1.3679 0.3679];

% o -1 significa dominio Z(discretizado)
G = tf(num, den, -1);

%rlocus(H, [27.825, 59.93])

H = feedback(G, 1);

rlocusplot(H, k2);
G
H