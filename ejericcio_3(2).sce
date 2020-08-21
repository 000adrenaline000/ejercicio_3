// CONTROL AEREOGENERADOR
// Alumno: MAmani Villanueva Jhojan Felipe
// Author: juan C. Cutipa-Luque
// Docente: juan C. Cutipa-Luque

//Codigo Reutilizado de un ejemplo Mostrado en Clase
clf();         // close current figure
clear          // clear all pasta variables
xdel(winsid()) // close all windows
Ap=[-50.424 -7.4588    416.4347   0        -926.3487   -450.407;
    -0.1845 -2.794     -18.8738   0        0           -29.2377;
    0       1          -0.2266   -0.0067   0           -0.219;
    1       0          0         0         0           0;
    0       0          0         0         0           0;
    0       0          0         0         0           0];
Bp=[0 0;
    0 0;
    0 0;
    0 0;
    1 0;
    0 1];
Cp=[0         0          0            1     0            0;
    0         0          -26.8638     0     0            -2.5965;
    1         0          0            0     0            0;
    0         1          0            0     0            0;
   -50.424    -7.4588    416.4347     0     -926.3487   -450.407;
   -0.1845    -2.794     -18.8738     0     0           -29.2377];
Dp=[0 0;
    0 0;
    0 0;
    0 0;
    0 0;
    0 0];
// 1. Represente la planta GN en matriz de transferencia indicando las magnitudes físicas de sus variables y el valor de sus parámetros

s=poly(0,"s")
S=s*eye(6,6);
SIAinv= inv(S-Ap)
// Calculate the transfer-function
H=Cp*SIAinv*Bp +Dp

ap_s=Ap;  
bp_s=Bp;
cp_s=Cp;
dp_s=Dp;

// Controllability and Observability
// Cc=[B, AB, A^2 B,..., A^(n-1) B]
Cc = cont_mat(ap_s,bp_s)
rankCc=rank(Cc)
//
// O=[C; CA; CA^2;...; CA^(n-1) ]
O = obsv_mat(ap_s, cp_s)
rankO=rank(O)
// verify if the rank of Cc is n, dimension of a
// verify if the rank of O is n, dimension of a

/*             Plot singular values of LTI the model                      */
G = syslin('c', ap_s, bp_s, cp_s, dp_s);
poles=spec(ap_s)
tr  = trzeros(G)

w = logspace(-3,3);
sv = svplot(G,w);

scf(1);
plot2d("ln", w, 20*log(sv')/log(10))
xgrid(12)
xtitle("Singular values plot","Frequency (rad/s)", "Amplitude (dB)");
ms=1.4;// 0.3;%1.5;    % guarantee overshot Mp < 6dB = 20*log10(2) 
wbs=0.23;//0.05;%0.23;
ee=1e-3;//1e-4
ki=1; // used to give more accurate adjustment to the cut-off frequency wbs
      // by default set it to 1
//           --------     WT Data    ------------
mt=1.3;//1.00;    % guarantee overshot Mp < 2dB = 20*log10(1.26)
wbt=4.1;//9.1;%4.1;
ee=1e-3;//1e-4
//           --------     WS     ------------
s=poly(0,'s');
ws1=(s/ms+wbs)/(s+wbs*ee),
ws2=ws1;
ws=[ws1,0;0,ws2]
Ws=syslin('c',ws)
//Ws=blockdiag(ws1,ws2)
//           --------     WT     ------------
s=poly(0,'s');
wt1=(s+wbt/mt)/(ee*s+wbt),
wt2=wt1;
wt=[wt1,0;0,wt2]
Wt=syslin('c',wt)
//Wt=blockdiag(wt1,wt2)
//           --------     WR     ------------
s=poly(0,'s');
wr1=s/s,
wr2=wr1;
wr=[wr1,0;0,wr2]
// ------------------ Plot weighting functions
svs = svplot(Ws,w);
svt = svplot(Wt,w);

//A partir de las barreras de estabilidad y desempeño, defina las funciones de ponderación WS y WT y plotee en relación a la frecuencia.

scf(2);
plot2d("ln", w, [-20*log(svs')/log(10) -20*log(svt')/log(10)])
xgrid(12)
xtitle("Singular values plot inv(Ws) and inv(Wt)","Frequency (rad/s)", "Amplitude (dB)");
