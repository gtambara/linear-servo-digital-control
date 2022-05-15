%  Universidade de Brasilia - Faculdade de Tecnologia
%  Laboratorio Controle Digital - semestre 2021/2
%  Alunos:
%       Gabriel Tambara Rabelo - 18/0017021
%       Pedro Pereira Nunes - 17/00

%  tendo obtido a matriz de dados x da simulação:

%  separando vetor temporal das medições
time = x(:,1);

%  separando vetor sinal de saída das medições
output = x(:,2);

% separando vetor de input das medições
input = x(:,3);
%% 

%  para corrigir os vetores tempo e saida para iniciarem com t=0:
%
aux = [0:0.001:0.032];
time = [aux time.'];
time = time.';
output = [zeros([1 33]) output.'];
output = output.';
aux = 15+zeros(1, 33);
input = [aux input.'];
input = input.';

%  buscando o valor máximo da saída
saida_max = max(output);

%  calculado o maximo sobressinal dada uma entrada de 15
maxss = (max(output)-15)/15;

%  buscando tempo de pico
tp = time(find(output==saida_max, 1));

%  calculando o fator de amortecimento
xi = -log(maxss)/(((pi^2)+((log(maxss))^2))^(0.5));

%  calculando a frequencia amortecida
wd = pi/tp;

%  calculando frequencia natural
wn = wd/((1-xi)^(1/2));

%  definindo a variavel s de frequencia complexa
s = tf('s');

%  modelando o formato geral para um sistema de segunda ordem se o erro
%sistema = (wn^2)/(s^2 + 2*s*xi*wn + wn^2);

%  como o vetor de dados inicia em 0.033 segundos, percebe-se que o modelo
%  tem um atraso, o seu sistema representando esse atraso seria na forma:
%
%  alternativa para sistema com vetores de entrada corrigidos (iniciados a
%  partir de t=0):
% 
sistema = ((wn^2)*exp(-s*0.033))/(s^2 + 2*s*xi*wn + wn^2);

%  definindo a entrada degrau para o vetor temporal
%t = 0.033:0.001:1
u = zeros(length(time),1);
u(time>=0) = 15;

%  declarando saida do sistema
output_teorico = lsim(sistema,u,time);

%  plotando o resultado
figure
plot(time, output, time, output_teorico)
hold on
title('Sistema modelado x real');
legend('medido','teorico');
xlabel('tempo');
ylabel('posi��o');

%  fazendo adaptações heurísticas para melhorar a modelagem, usa-se:
wn_adpt = 24.5;
xi_adpt = 0.253;
%sistema_adpt = (wn_adpt^2)/(s^2 + 2*s*xi_adpt*wn_adpt + wn_adpt^2);
%  alternativa para sistema com vetores de entrada corrigidos (iniciados a
%  partir de t=0)

sistema_adpt = ((wn_adpt^2)*exp(-s*0.033))/(s^2 + 2*s*xi_adpt*wn_adpt + wn_adpt^2);
output_teorico_adpt = lsim(sistema_adpt,u,time);

%  plotando o resultado
figure
plot(time, output, time, output_teorico_adpt)
hold on
title('Sistema modelado adaptado x real');
legend('medido','teorico adaptado');
xlabel('tempo');
ylabel('posi��o');

%  declarando sistema com modelagem obtida por fitting via metodo dos minimos
%  quadrados
T = 0.001;
np = 2;
nz = 0;
data = iddata(output, input, T);
iodelay = 0.033; %atraso conhecido do sistema
tfmmq = tfest(data,np,nz,iodelay);
output_mmq = lsim(tfmmq,u,time);

%  plotando para comparacao das 3 curvas
figure
plot(time, output, time, output_teorico_adpt, time, output_mmq)
hold on
title('Sistema modelado adaptado x por minimos quadrados x real');
legend('medido','teorico adaptado','minimos quadrados');
xlabel('tempo');
ylabel('posi��o');

%  avaliacao de qual metodo possui menor erro em relacao a curva medida
fitmmq = goodnessOfFit(output_mmq,output,'NRMSE');
fit = goodnessOfFit(output_teorico_adpt,output,'NRMSE');

%  realização canonica controlável de sistema_adpt
T = T*11;
Ganho = 420;
b = (wn_adpt^2)/Ganho;

%  definindo realização com estados x e x':
A = [0 1; 0 -2*xi_adpt*wn_adpt];
B = [0 ;b];
C = [1 0];
D = 0;

% sistema_adpt_base = ss(A, B, C, D);
sistema_adpt_base = (b*exp(-s*0.033))/(s^2 + 2*s*xi_adpt*wn_adpt);
sistema_adpt_d1 = c2d(sistema_adpt_base, T,'zoh');
sistema_adpt_d = canon(sistema_adpt_d1, 'companion');
%sistema_adpt_d = canon(sistema_adpt_d, 'companion'); %retorna na forma canonica observável

G = sistema_adpt_d.A;
H = sistema_adpt_d.B;
C = sistema_adpt_d.C;
D = sistema_adpt_d.D;

%  checar controlabilidade (verdadeira se diferente de 0)
dc = det(ctrb(sistema_adpt_d));

%  checar observabilidade (verdadeira se diferente de 0)
do = det(obsv(sistema_adpt_d));

%DETERMINAR OS P�LOS DESEJADOS COM BASE NAS CARACTER�STICAS DO SISTEMA)
Ts = 1; %Settling Time;
Tp = 0.2; %Peak Time;
Mp = 0.25; %Sobressinal;
i = 3; %Seleciona qual dois par�metros ser�o utilizados
if i==1 %Ts e Tp
    sigma = 4/Ts;
    omega_d = pi/Tp;
elseif i==2 %Mp e Tp
    a = -(pi^2)-(log(Mp))^2;
    c = (log(Mp))^2;
    zeta = abs(roots([a 0 c]));
    omega_d = pi/Tp;
    omega_n = omega_d/(sqrt(1-zeta(1)^2));
    sigma = zeta(1)*omega_n;
else %Ts e Mp
    a = -(pi^2)-(log(Mp))^2;
    c = (log(Mp))^2;
    zeta = roots([a 0 c]);
    sigma = 4/Ts;
    omega_n = sigma/zeta(1);
    omega_d = omega_n*sqrt(1-zeta(1)^2);   
end

r = exp(-sigma*T);
teta = omega_d*T;
polo1 = r*cos(teta)+r*sin(teta)*j;
polo2 = r*cos(teta)-r*sin(teta)*j;

%P�LOS E MATRIZ DO OBSERVADOR
V = 0.1; %Qu�o mais r�pido ser� o p�lo do observador
r_obs = V*r;
polo1_obs = r_obs*cos(teta)+r_obs*sin(teta)*j;
polo2_obs = r_obs*cos(teta)-r_obs*sin(teta)*j;

polos_obs = [polo1_obs polo2_obs];
L = place(G', C', polos_obs);
L = L';

%DADOS OS P�LOS, CALCULAR MATRIZES K1 e K2
%PARA CONTROLADOR PI

G_chapeu = [G H; 0 0 0];
H_chapeu = [0; 0; 1];
polos = [polo1 polo2 0];
K_chapeu = place(G_chapeu, H_chapeu, polos);
I = eye(2);
M_aux = [G-I H;C*G C*H];
M_I = [0 0 1];
K = (K_chapeu+M_I)*inv(M_aux);
K1 = K(3);
K2 = [K(1) K(2)];

GCon = G;
VMAX_AMP = 8;

sim('Teorico_plus_estimado');

%  plotando o resultado
figure
plot(sistema_teorico);
hold on
plot(sistema_modelado);
hold on
plot(acao_controle);
title('Controlador PI');
legend('Planta Te�rica com Observador','Planta identificada', 'A��o de Controle');
xlabel('tempo');
ylabel('posi��o');