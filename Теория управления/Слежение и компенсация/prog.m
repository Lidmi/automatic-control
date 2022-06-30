% Исходные данные
A=[0 5/4; 1 6/4]
B=[-1; 3/4]
C=[0 1]
Bf=[-2; 5/4];

% Расширенная модель ошибок
G=[1;0];
Ge=[1;0];
Gamma_f=[0 1;-81 0];
H_f=[1 0];
Gamma_g=[0 1;-4 0];
H_g=[1 0];
A_bar=[Gamma_f zeros(2,2) G*C; zeros(2,2) Gamma_g -Ge*C; zeros(2,4) A];
B_bar=[zeros(4,1); B];

% Вычисление матрицы линейных стационарных обратных связей
Gamma_bar=[-5.2 1 0 0 0 0;
           -1 -5.2 1 0 0 0;
           0 0 -5.1 1 0 0;
           0 0 0 -5.9 1 0;
           0 0 0 0 -5.8 1;
           0 0 0 0 0 -5.5];
H_bar=[1 0 1 1 1 1];
M_bar=sylv(-A_bar, Gamma_bar, B_bar*H_bar);
K_bar=-H_bar*inv(M_bar)

%Наблюдатель
Gamma=[-5.25 0; 0 -5.95];
H=[1 1];
M=sylv(-A', Gamma, C'*H);
L=(-H*inv(M))';
F_n=A-L*C;
eig(F_n)

% Проверка корней
F_bar=A_bar-B_bar*K_bar;
eig(F_bar)

A_bar_2=[Gamma_f zeros(2,2) G*C zeros(2,2); zeros(2,2) Gamma_g -Ge*C zeros(2,2); zeros(2,4) A zeros(2,2); zeros(2,4) L*C A-L*C];
B_bar_2=[zeros(4,1); B; B];
% K_bar_2=[K_bar zeros(1,2)];
K_bar_2=[K_bar(1) K_bar(2) K_bar(3) K_bar(4) 0 0 K_bar(5) K_bar(6)]

F_bar_2=A_bar_2-B_bar_2*K_bar_2;
eig(F_bar_2)