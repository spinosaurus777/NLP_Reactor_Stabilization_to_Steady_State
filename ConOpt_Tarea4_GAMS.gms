*
*CONTROL ÓPTIMO DE PROCESOS
*12/09/2023
*TAREA #4

*Conjuntos
Set
    i Discretización del Tiempo /0*20/;
    
*Parametros
Parameters
    N Discretización /20/
    C_init Concetración Inicial /0.1367/
    T_init Temperatura Incial /0.7293/
    alpha Constante - Balance de Energía /1.95e-4/
    theta Constante - Balances /20/
    Tf Temperatura de Entrada /0.3947/
    Tc Temperatura del Intercambiador /0.3916/
    k10 Constante de Velocidad de Reacción /10/
    eta Constante /1/
    thao Tiempo Adimensional /10/
    C_des Concetración en EE /0.0194/
    T_des Temperatura en EE /0.8/
    u_des Flujo del Intercambaidor en EE /340/
    alpha_1 Constante para la concetracion FO /1e4/
    alpha_2 Constante para la temperatura FO /2e4/
    alpha_3 Constante para el flujo de refrigerante FO /1e-3/;

*Variables
Variables
    C(i) Concetración en el tiempo i
    T(i) Temperatrua en el tiempo i
    w(i) Coso en el tiempo i
    C_1(i) Diferencias en el tiempo i
    T_1(i) Diferencias en el tiempo i
    u(i) Flujo del Intercambiador en el tiempo i
    fobj(i) Función Objetivo a cada paso
    sum_fobj FO;
    
Positive Variables u(i);

*Restricciones
Equations
    Bal_Masa(i) Balance de Masa
    Bal_Energia(i) Balance de Energia
    diff_C(i) Diferencias de C
    diff_T(i) Diferencias de T
    eq1(i) Ecuacion de w
    cal_fobj(i) Funcion objetivo a cada paso
    cal_sum_fobj Sumatoria funcion objetivo;
 
u.lo(i)=0;
u.up(i)=500;
C_1.l(i)=1;
T_1.l(i)=1;
u.l(i)=250;
C.l(i)=C_init+(((C_des-C_init)*(ord(i)-1))/N);
T.l(i)=T_init+(((T_des-T_init)*(ord(i)-1))/N);
T.fx("0")=T_init;
C.fx("0")=C_init;

eq1(i).. w(i)*T(i)+eta=e=0;
diff_C(i)$(ord(i)>1).. (C(i)-(C(i-1)))/(thao/N)=e=C_1(i);
diff_T(i)$(ord(i)>1).. (T(i)-(T(i-1)))/(thao/N)=e=T_1(i);
Bal_Masa(i).. C_1(i)=e=((1/theta)*(1-C(i)))-(k10*exp(w(i))*C(i));
Bal_Energia(i).. T_1(i)=e=((1/theta)*(Tf-T(i)))+(k10*exp(w(i))*C(i))-(alpha*u(i)*(T(i)-Tc));
cal_fobj(i).. fobj(i)=e=(alpha_1*power((C_des-C(i)),2))+(alpha_2*power((T_des-T(i)),2))+(alpha_3*power((u_des-u(i)),2));
cal_sum_fobj.. sum_fobj=e=sum(i,fobj(i));

*Solucion
model CSTR /all/;
options nlp=ipopt;
solve CSTR using nlp minimizing sum_fobj

