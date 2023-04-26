function [sys, x0, str, ts] = sliding_function(t,x,u,flag)
    switch flag
    case 0
        [sys, x0, str,ts] = mdlInitializeSizes;
    case 1
        sys = mdlDerivatives(t,x,u);
    case 3
        sys = mdlOutputs(t,x,u);
    case {2,4,9}
        sys = [];
    otherwise
        error(['Unhandled flag = ', num2str(flag)]);
    end
end
function [sys,x0,str,ts]=mdlInitializeSizes
    sizes = simsizes;
    sizes.NumContStates = 0;
    sizes.NumDiscStates = 0;
    sizes.NumOutputs = 4;
    sizes.NumInputs =18 ;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 0;
    sys = simsizes(sizes);
    x0 = [];
    str = [];
    ts = [];
end


function sys = mdlOutputs(t,x,u)
 
% ------------------------- Parameters Setup ------------------------------
  alphaa = diag([1,1]); alphau = diag([0.1,0.1]) ; lamdaa = diag([0.5,0.5]) ; lamdau = diag([-2.5,-2.5]);
% ------------------------- Reference ------------------------------
  xr = u(13);
  dxr = u(14);
  ddxr = u(15);
  gammar = u(16); dgammar = u(17) ; ddgammar = u(18);
  a = u(1) ; da = u(2) ; dda = u(3);
  gamma = u(4) ; dgamma = u(5) ; ddgamma = u(6);
  theta = u(7) ; dtheta = u(8) ; ddtheta = u(9);
  phi = u(10) ; dphi = u(11) ; ddphi = u(12);
  qa = [a; gamma] ; qu = [theta; phi]; dqa = [da;dgamma]; dqu = [dtheta;dphi]; ddqa = [dda;ddgamma]; ddqu = [ddtheta;ddphi];
  qar = [xr;gammar]; qur = [0;0] ; dqar = [dxr;dgammar]; dqur = [0;0]; ddqar = [ddxr;ddgammar]; ddqur = [0;0];
  
  dqr = alphaa*dqar - lamdaa*(qa - qar) - lamdau*(qu - qur);
  ddqr = alphaa*ddqar - lamdaa*(dqa - dqar) - lamdau*(dqu - dqur);
  
%------------------------Sliding function and its derivative -------------
  s = alphaa*dqa + alphau*dqu - dqr;
  ds = alphaa*ddqa + alphau*ddqu - ddqr;

%-----------------------------------------------------------------------  
  sys(1) = s(1);
  sys(2) = s(2);
  sys(3) = ds(1);
  sys(4) = ds(2);


end
