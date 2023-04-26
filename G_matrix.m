function [sys, x0, str, ts] = G_matrix(t,x,u,flag)
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
    sizes.NumOutputs = 6;
    sizes.NumInputs =18 ;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 0;
    sys = simsizes(sizes);
    x0 = [];
    str = [];
    ts = [];
end

function sys = mdlOutputs(t,x,u)
 %------------------------- Parameters Setup ------------------------------
  alphaa = diag([1,1]); alphau = diag([0.1,0.1]) ; lamdaa = diag([0.5,0.5]) ; lamdau = diag([-2.5,-2.5]);
  xr = u(13);dxr = u(14);ddxr = u(15);
  gammar = u(16); dgammar = u(17) ; ddgammar = u(18);
  a = u(1) ; da = u(2) ; dda = u(3);gamma = u(4) ; dgamma = u(5) ; ddgamma = u(6);
  theta = u(7) ; dtheta = u(8) ; ddtheta = u(9);phi = u(10) ; dphi = u(11) ; ddphi = u(12);
  qa = [a; gamma] ; qu = [theta; phi]; dqa = [da;dgamma]; dqu = [dtheta;dphi]; ddqa = [dda;ddgamma]; ddqu = [ddtheta;ddphi];
  qar = [xr;gammar]; qur = [0;0] ; dqar = [dxr;dgammar]; dqur = [0;0]; ddqar = [ddxr;ddgammar]; ddqur = [0;0];
  
  dqr = alphaa*dqar - lamdaa*(qa - qar) - lamdau*(qu - qur);
  ddqr = alphaa*ddqar - lamdaa*(dqa - dqar) - lamdau*(dqu - dqur);
  
  
 %-------------------------------G Matrix Calculation----------------------------------------------------------
sth = sin(theta) ; cth = cos(theta);
sph = sin(phi) ; cph = cos(phi);
%--------------------------------------------------------------------------------------
mc = 2.5; mt =10; J =6; l=3;
m11 = mt+ mc;
m22 = J +(mt+mc)*a^2 + mc*l^2*sth^2 + mc*l^2*sph^2*cth^2 + 2*mc*l*sph*cth*a;
m33 = mc*l^2*cth^2; m44 = mc*l^2;m21 = -mc*l*sth; m12 = -mc*l*sth;
m23 = -mc*l^2*cph*cth*sth; m32 = -mc*l^2*cph*cth*sth;
m24 = mc*l*(a*cth+l*sph); m42 =  mc*l*(a*cth+l*sph); 
m13 = mc*l*cph*cth; m31 = mc*l*cph*cth;m14 = -mc*l*sph*sth; m41 = -mc*l*sph*sth;
M11 = [m11 m12;m21 m22];M12 = [m13 m14 ;m23 m24];
M21 = [m13 m23 ;m14 m24];M22 = [m33 0;0 m44];
M1222 = M12/(M22);M2111 = M21/(M11);M1_ = M11 - M1222*M21;M2_ = M22 - M2111*M12;
Da = inv(M1_);Du = -inv(M2_)*M21*inv(M11);
G = Da+ alphau*Du;
%-------------------------------------------------------------------------------------


%-----------------------------------------------------------------------  
  sys(1) = G(1,1);
  sys(2) = G(1,2);
  sys(3) = G(2,1);
  sys(4) = G(2,2);
  sys(5) = ddqr(1);
  sys(6) = ddqr(2);

end

