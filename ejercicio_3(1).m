clear

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

sys=ss(Ap,Bp,Cp,Dp);

[ns,nc]=size(Bp)  %ns= number of inputs; nc=number of controls
%checking controllability and observability
[i,j] = size(Ap);
%e=[B, AB, A^2 B,..., A^(n-1) B]
e = ctrb(Ap,Bp);
rankC=rank(e);

if i == rankC 
    disp('Continuous System is Controllable')
end
%o=[C; CA; CA^2;...; CA^(n-1) ]
o = obsv(Ap, Cp);
rankO=rank(o);
if j == rankO 
    disp('Continuous System is Observable')
end

% 2. Calculo de los polos y ceros de tranmision.
    [p,z]=pzmap(sys)
    % plot the poles and zeros
    figure; pzmap(sys)


% 3. Calculo de las barreras de estabilidad y desempeño
    w=logspace(-1,3,400);
    a=50;b=200;c=150;
    x1=20*ones(1,a); x2=60*zeros(1,b); x3=-20*ones(1,c);
    xt=[x1 x2 x3];
    x4=[10*ones(1,80) 0*zeros(1,320)];
    % barreras de estabilidad

    semilogx(w, xt,'r')
    ylim([-60,60])
    grid;
    hold on;
    semilogx(w,x4,'b')
    title('Barreras de Estabilidad')
    xlabel('w(rad/s)')
    ylabel('dB')

    hold off;

% 4. Calculo de los valores singulares de la planta
    w=logspace(-2,2,100);
    sv=sigma(sys,w);
    vsmax(1)=max(sv(1,:));    vsmin(1)=min(sv(1,:));
    vsmax(2)=max(sv(2,:));    vsmin(2)=min(sv(2,:));
    vsmax     %vector de valores singulares maximos
    vsmin     %vector de valores singulares min

% 5. Plotee los valores singulares máximos y mínimos de GN en relación a la frecuencia. 
    % compute singular values and plot
    w=logspace(-2,2,100);
    svalues=sigma(sys,w);
    V=svalues(:,51); %VALORES SINGULARES EN 1 RAD/S
    svalues=20*log10(sv);
    figure;
    semilogx(w,svalues)

% 6. Diseñe un controlador óptimo LQR (lineal cuadrático), asumiend que conoce todos los estados
    %QL=[5 0;0 2];
    QL=Cp'*Cp;
    %Q=eye(4,4);
    rho=1e-1; %Cheap control recovery parameter
    %RL=[4 0;0 1]
    RL = rho*eye(nc); %Control Weigthing Matrix
    [G, ~, ~] = lqr(sys,QL,RL)
%G=[-1.085 0.1660 -0.3516 0.5581 -0.007 0.1214;
%-0.5397 -1.3119 -0.1294 -7.1731 -0.0001 -0.7172];
% 7. Diseñe un observador de estados determinístico considerando que solamente es posible medir algunos estados de la planta

    QO=eye(6);
    RO=0.01*eye(6);
    [H, ~, ~] = lqe(sys,QO,RO)
    
% 8. Con el controlador y el observador proyectados anteriormente, forme la estructura del compensador K(s)

    K11=[Ap-Bp*G];  
    K12=[-Bp*G];    
    K21=[H*Cp]; 
    K22=[Ap-H*Cp-Bp*G];   

    AK=[K11 K12;K21 K22];  
    BK=[zeros(6,6);-H];    
    CK=[Cp zeros(6,6)];     
    DK=0;
 
    sysK=ss(AK,BK,CK,DK)
    
    
% 9. Plotee las curvas de sensibilidad S y T y de la función de malla GNK. Evalúe si estas respetan las barreras de desempeño (realice gráficos superpuestos)

    al=[zeros(2,4) zeros(2,4); Bp Ap];
    bl=[eye(2,2) zeros(2,4) ;zeros(6,6)];
    cl=[zeros(6,2) Cp];
    dl=0;

    w=logspace(-2,3,100);
    sv = sigma(ss(al-bl*cl, bl, -cl,dl),w);
    sv = 20*log10(sv);
    figure
    semilogx(w, sv)
    title('Sensitivity Singular Values')
    grid
    xlabel('Frequency (rad/sec)')
    ylabel('Singular Values (dB)')


    svc = sigma(ss(al-bl*cl, bl, -cl, 0*dl),w);
    svc = 20*log10(sv);
    figure
    semilogx(w,svc)
    title('Complementary Sensitivity Singular Values')
    grid
    xlabel('Frequency (rad/sec)')
    ylabel('Singular Values (dB)')

% 10. Plotee la respuesta del sistema controlado en función del tiempo frente a una entrada escalón unitario, quiere decir el compensador K y la planta GN en realimentación
    OL=sysK*sys;
    temp=eye(6,6) +[ OL zeros(6,4)];
    CL=[zeros(6,4) OL]*inv(temp);
    t=(0:0.1:10);
    [NC DC]=ss2tf(CL);
    TF1=tf(NC(1,1),DC(1,1));
    TF2=tf(NC(1,2),DC(1,2));
    TF3=tf(NC(2,1),DC(2,1));
    TF4=tf(NC(2,2),DC(2,2));
    t=[0:0.1:5];

    figure
    subplot(2,2,1)
    step(TF1,t)
    title('Funcion Transferencia INPUT 1 OUTPUT 1')
    subplot(2,2,2)
    step(TF2,t)
    title('Funcion Transferencia INPUT 1 OUTPUT 2')
    subplot(2,2,3)
    step(TF3,t)
    title('Funcion Transferencia INPUT 2 OUTPUT 1')
    subplot(2,2,4)
    step(TF4,t)
    title('Funcion Transferencia INPUT 2 OUTPUT 2')