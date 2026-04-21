% RHS of the bacteria-phage-immune ODE system.
% Implements Eqs. (1)-(3) of Leung & Weitz, J. Theor. Biol. 429, 241 (2017).
%   dB/dt : logistic bacterial growth, phage lysis (mass action),
%           saturating immune killing (half-max at B = KD)
%   dP/dt : phage production (burst size beta) minus first-order decay
%   dI/dt : logistic immune activation, gated by bacterial antigen (KN)
% State vector y = [B; P; I] in units of ml^-1.
% Parameters passed via struct `para`; see fig3_synergy.m or fig5_heat.m
% for values and units.

function dy=infection(t,y,para)
    dy=zeros(3,1);
    bpop=y(1);
    bpop(bpop<para.thres)=0;
    vpop=y(2);
    vpop(vpop<para.thres)=0;
    imm=y(3);
    
    % Bacteria population change
    dy(1)=para.r*bpop*(1-bpop/para.KC)- ...
    bpop*para.phi*vpop- ...
    para.eps*imm*bpop/(1+bpop/para.KD); 
    
    % Viruses population change
    dy(2)=para.beta*para.phi*vpop*bpop ...
    -para.omg*vpop;
    
    % Immune response change
    dy(3)=para.alpha*imm*(1-imm/para.KI)* ...
    bpop/(bpop+para.KN);
