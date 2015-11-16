%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                     TER: Commande avancée                       %%%% 
%%%%             Groupe 1 : Maingot Bouissiere Bedu                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Paramètre simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_echant=10^-2;
start=0;
stop=200;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Configuration simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SIMULATION
BF=1; %%%% 1 =BO
      %%%% 0 =BF
PI=0; %%%% 1:PI
      %%%% 0:P

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Définition des paramètres du système
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ka=0.1;                  % Gain CVS [A/V]
Kc=5.7e-2;               % Constante couple courant du moteur [N.m/I]
pas=5*10^-3;             % pas de vis [m]
Umax=15;%15       % Saturation commande [V]
Umin=-Umax;              % 
Bmmax=11111;%0.35        % Butée mécanique du système [m]
Bmmin=-Bmmax;            %
f=1.52*10^-4;            % Frotemment []
m=1.3;                   % Masse plateau + vis [kg]
Jr=4.5*10^-5;            % Inertie du moteur et de la vis raménée [] 
M=30;                    % Masse posée sur le plateau         
% Calcul inertie totale du système ramenée au moteur
J=Jr+(m+M)*(pas/(2*pi))^2;
% Schema bloc equivalent
Tau=J/f;
Kred=pas/(2*pi);         % Coefficient de réduction
K=1/f;
Kt=Kred*Ka*Kc*K;        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Consigne
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ech=0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Correcteur proportionnel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kprop=100;
Ki=1;
Ti=1;
%% Algorithme de détermination du gain proportionnel idéal 
%%%%% Nombre d'itération fixé. On trouve donc le meilleur gain pour ce nombre d'itération
Kprop_max=100;
Kprop_min=0;
DE=[];
n=0;
while (n<10 )
K_current=(Kprop_max+Kprop_min)/2;
Kprop=K_current;
sim 'modele'

%calcul depassement
D=(max(dep_plat)-dep_plat(length(dep_plat)))/dep_plat(length(dep_plat));
    if D > 0
        Kprop_max=K_current;
    else
        Kprop_min=K_current;
        n=n+1;
    end
 DE=[DE D]
 

end
%% Autre methode
%Kprop=1/(4*Kt*Tau);
%% Correcteur PI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Marge de phase °
Mphi=89;
phim=pi/180*Mphi;
a=(1+sin(phim))/(1-sin(phim));
Ti=a*Tau;
Ki=1/(sqrt(a)*Tau*Kt);
%% Algorithme %%%%%%%%%%%%%%%%% MARCHE PAS ENCORE %%%%%%%%%%%%%%%%%%%%%%%
%% Probleme dépacement jms nul.
%%%% Nombre d'itération fixé. On trouve donc le meilleur gain pour ce nombre d'itération
% phim_max=pi/2
% phim_min=0;
% DE=[];
% n=0;
% while (n<10 )
% phim_current=(phim_max+phim_min)/2;
% a=(1+sin(phim_current))/(1-sin(phim_current));
% Ti=a*Tau;
% Ki=1/(sqrt(a)*Tau*Kt);
% sim 'modele'
% 
% %calcul depassement
% D=(max(dep_plat)-dep_plat(length(dep_plat)))/dep_plat(length(dep_plat));
%     if D > 0
%         phim_max=phim_current;
%     else
%         phim_min=phim_current;
%         n=n+1;
%     end
%  DE=[DE D]
%  
% 
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Essais
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Robustesse Correcteur proportionnel
% for n=1:1:5
% M=[0 10 15 20 30];
% J=Jr+(m+M(n))*(pas/(2*pi))^2;
% Tau=J/f;
% Kprop=1/(4*Kt*Tau);  
% figure(n)
% coul=['r','g','y','b','m'];
% d_tab=[];
% for i=1:1:5
% J=Jr+(m+M(i))*(pas/(2*pi))^2;
% Tau=J/f;
% sim 'modele'
% 
% hold on
%   plot(t,dep_plat,coul(i)); 
%   grid on;
% end
% legend('M=0','M=10','M=15','M=20','M=30'); 
% xlabel('Temps (sec)')
% ylabel('Déplacement (m)')
% 
% str=sprintf('Robustesse du correcteur vis à vis d''variation de masse (Correcteur synthétisé pour M=%d kg)',M(n));
% title(str)
% end
%% Robustesse Correcteur PI
% for n=1:1:5
% M=[0 10 15 20 30];
% J=Jr+(m+M(n))*(pas/(2*pi))^2;
% Tau=J/f;
% Ti=a*Tau;
% Ki=1/(sqrt(a)*Tau*Kt); 
% figure(n)
% coul=['r','g','y','b','m'];
% d_tab=[];
% for i=1:1:5
% J=Jr+(m+M(i))*(pas/(2*pi))^2;
% Tau=J/f;
% sim 'modele'
% 
% hold on
%   plot(t,dep_plat,coul(i)); 
%   grid on;
% end
% legend('M=0','M=10','M=15','M=20','M=30'); 
% xlabel('Temps (sec)')
% ylabel('Déplacement (m)')
% 
% str=sprintf('Robustesse du correcteur vis à vis d''variation de masse (Correcteur synthétisé pour M=%d kg)',M(n));
% title(str)
% end


%% Essai Correcteur PI
%  sim 'modele'
% figure(1)
% plot(t,U,'m',t,dep_plat,'b')
% grid on
 
%figure(1)
% subplot (2,2,1);
% 
% plot (t,U);
% grid on;
% legend ('Commande');
% 
% subplot(2,2,3);
% 
% % plot(t,dep_pot);
% % grid on;
% % legend ('Deplacement pot');
% tm5=[0 stop]
% vm5=0.95*[Omega_mot1(stop/t_echant) Omega_mot1(stop/t_echant)]
% 
% plot (t,Omega_mot1,'r')
% grid on;
% legend ('rotation_moteur');
% 
% 
% subplot (2,2,2);
% 
% plot (t,dep_plat,'r')
% grid on;
% legend ('Deplacement chariot');
% 
% subplot (2,2,4);
% 
% plot (t,I_mot);
% grid on;
% legend ('Courant');
% 
% figure(2)
% plot (t,U,'--b');
% 
% hold on
% plot (t,dep_plat,'r')
% grid on;
% legend ('Commande','Deplacement pot');
% 
% % figure(3)
% % plot (t,Omega_mot1,'r')
% % hold on
% % plot(tm5,vm5,'r')
% % hold off
% % grid on;
% % legend ('ref 95%','rotation_{moteur}');
% % title('recherche cte temps systeme pur' )
% 
% % %Bode
% % Ktot=Kprop*Kred*Ka*(1/f)*Kc;
% % F=tf([Ktot],[J/f 1 0]);
% % bode(F);
% % grid on;