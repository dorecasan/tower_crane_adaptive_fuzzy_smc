function [sys, x0, str, ts] = plant_tc_F2(t,x,u,flag)
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
sizes.NumContStates = 8 ;
sizes.NumDiscStates = 0;
sizes.NumOutputs = 10;
sizes.NumInputs = 2;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 0;
sys=simsizes(sizes); 
x0 = [zeros(8,1)];
str = [];
ts = [];
end

function sys = mdlDerivatives(t,x,u)
dq1 = [x(2);x(4)] ; dq2 = [x(6); x(8)];

%--------------------------------------------------------------------------------------
mc = 2.5; mt =10; J =6; l=3;
bx =0; bgamma = 0; btheta = 0; bphi = 0; g =9.8;

a = x(1) ; gamma = x(3) ; phi = x(5) ; theta = x(7);
da = x(2) ; dgamma = x(4) ; dphi = x(6) ; dtheta = x(8);
sth = sin(theta) ; cth = cos(theta);
sph = sin(phi) ; cph = cos(phi);
%--------------------------------------------------------------------------------------
m11 = mt+ mc;
m22 = J +(mt+mc)*a^2 + mc*l^2*sth^2 + mc*l^2*sph^2*cth^2 + 2*mc*l*sph*cth*a;
m33 = mc*l^2*cth^2; m44 = mc*l^2;
m21 = -mc*l*sth; m12 = -mc*l*sth;
m23 = -mc*l^2*cph*cth*sth; m32 = -mc*l^2*cph*cth*sth;
m24 = mc*l*(a*cth+l*sph); m42 =  mc*l*(a*cth+l*sph); 
m13 = mc*l*cph*cth; m31 = mc*l*cph*cth;
m14 = -mc*l*sph*sth; m41 = -mc*l*sph*sth;

M11 = [m11 m12;m21 m22];
M12 = [m13 m14 ;m23 m24];
M21 = [m13 m23 ;m14 m24];
M22 = [m33 0;0 m44];
M = [M11 M12; M21 M22];
%--------------------------------------------------------------------------------------
c22 = bgamma+ mc*l*(sph*sth*da + cph*cth*a*dphi - sph*sth*a*dtheta) + mc*l^2*(sph*cph*cth^2*dphi + cph^2*sth*cth*dtheta)+(mc+mt)*a*da;
c21 = dgamma*((mt+mc)*a + mc*l*sph*cth);
c23 = mc*l*a*cph*cth*dgamma+mc*l^2*sph*sth*cth*dphi+mc*l^2*sph*cph*cth^2*dgamma+mc*l^2*cph*sth^2*dtheta;
c24 = -mc*l*a*sph*sth*dgamma - mc*l*a*sth*dtheta + mc*l^2*cph^2*sth*cth*dgamma + mc*l^2*cph*sth^2*dphi;
c12 = -(mc*l*cth*(dtheta + sph*dgamma) + (mt+mc)*a*dgamma);
c11 = bx ;c13 = -mc*l*(sph*cth*dphi + cph*sth*dtheta);
c14 = -mc*l*(cph*sth*dphi + sph*cth*dtheta + cth*dgamma);
c32 = -mc*l*cph*cth*(a+l*sph*cth)*dgamma - mc*l^2*cph*cth^2*dtheta;
c34 = -mc*l^2*(cph*cth^2*dgamma + sth*cth*dphi);
c42 = mc*l*(a*sph*sth - l*cph^2*sth*cth)*dgamma + mc*l^2*cph*cth^2*dphi + mc*l*cth*da;
c43 = mc*l^2*(cph*cth^2*dgamma + sth*cth*dphi);
c41 = mc*l*cth*dgamma;c33 = bphi - mc*l^2*sth*cth*dtheta;
c44 = btheta;
C11 = [c11 c12; c21 c22]; C12 = [c13 c14; c23 c24]; C21 = [0 c32; c41 c42];
C22 = [c33 c34; c43 c44]; 
C = [C11 C12;C21 C22];

%--------------------------------------------------------------------------------------
g1 = mc*g*l*cth*sph;
g2 = mc*g*l*sth*cph;
G2 = [g1; g2];
G = [0;0;g1;g2];
%--------------------------------------------------------------------------------------
M1222 = M12/(M22);
M2111 = M21/(M11);
%--------------------------------------------------------------------------------------
M1_ = M11 - M1222*M21;
C11_ = C11 - M1222*C21;
C12_ = C12 - M1222*C22;
G1 = -M1222*G2;
M2_ = M22 - M2111*M12;
C21_ = C21 - M2111*C11;
C22_ = C22 - M2111*C12;
%--------------------------------------------------------------------------------------
invM1_ = inv(M1_);
invM2_ = inv(M2_);
%--------------------------------------------------------------------------------------

U = [u(1);u(2);0;0];

nx = 0.1*sin(t); ngamma = 0.2*sin(t) ; nphi = 0.1*sin(t) ; ntheta = 0.1*sin(t);
N = [nx;ngamma;nphi;ntheta];


%ddq1 = invM1_*(F1 - C11_*dq1 - C12_*dq2 - G1);
%ddq2 = invM2_*(F2 - C21_*dq1 - C22_*dq2 - G2);
ddq = inv(M)*( U - C*[dq1;dq2] -G - N);


sys(1) = da;
sys(2) = ddq(1);
sys(3) = dgamma;
sys(4) = ddq(2);
sys(5) = dphi;
sys(6) = ddq(3);
sys(7) = dtheta;
sys(8) = ddq(4);
end

function sys = mdlOutputs(t,x,u)

alphaa = diag([1,1]); alphau = diag([0.1,0.1]) ; lamdaa = diag([0.5,0.5]) ; lamdau = diag([-2.5,-2.5]);
dq1 = [x(2);x(4)] ; dq2 = [x(6); x(8)];
%--------------------------------------------------------------------------------------
mc = 2.5; mt =10; J =6; l=3;
bx =0; bgamma = 0; btheta = 0; bphi = 0; g =9.8;
a = x(1) ; gamma = x(3) ; phi = x(5) ; theta = x(7);
da = x(2) ; dgamma = x(4) ; dphi = x(6) ; dtheta = x(8);
sth = sin(theta) ; cth = cos(theta);
sph = sin(phi) ; cph = cos(phi);
%--------------------------------------------------------------------------------------
m11 = mt+ mc;
m22 = J +(mt+mc)*a^2 + mc*l^2*sth^2 + mc*l^2*sph^2*cth^2 + 2*mc*l*sph*cth*a;
m33 = mc*l^2*cth^2; m44 = mc*l^2;m21 = -mc*l*sth; m12 = -mc*l*sth;
m23 = -mc*l^2*cph*cth*sth; m32 = -mc*l^2*cph*cth*sth;
m24 = mc*l*(a*cth+l*sph); m42 =  mc*l*(a*cth+l*sph); 
m13 = mc*l*cph*cth; m31 = mc*l*cph*cth;m14 = -mc*l*sph*sth; m41 = -mc*l*sph*sth;

M11 = [m11 m12;m21 m22];M12 = [m13 m14 ;m23 m24];
M21 = [m13 m23 ;m14 m24];M22 = [m33 0;0 m44];
M = [M11 M12; M21 M22];
%--------------------------------------------------------------------------------------
c22 = bgamma+ mc*l*(sph*sth*da + cph*cth*a*dphi - sph*sth*a*dtheta) + mc*l^2*(sph*cph*cth^2*dphi + cph^2*sth*cth*dtheta)+(mc+mt)*a*da;
c21 = dgamma*((mt+mc)*a + mc*l*sph*cth);
c23 = mc*l*a*cph*cth*dgamma+mc*l^2*sph*sth*cth*dphi+mc*l^2*sph*cph*cth^2*dgamma+mc*l^2*cph*sth^2*dtheta;
c24 = -mc*l*a*sph*sth*dgamma - mc*l*a*sth*dtheta + mc*l^2*cph^2*sth*cth*dgamma + mc*l^2*cph*sth^2*dphi;
c12 = -(mc*l*cth*(dtheta + sph*dgamma) + (mt+mc)*a*dgamma);
c11 = bx ;
c13 = -mc*l*(sph*cth*dphi + cph*sth*dtheta);
c14 = -mc*l*(cph*sth*dphi + sph*cth*dtheta + cth*dgamma);
c32 = -mc*l*cph*cth*(a+l*sph*cth)*dgamma - mc*l^2*cph*cth^2*dtheta;
c34 = -mc*l^2*(cph*cth^2*dgamma + sth*cth*dphi);
c42 = mc*l*(a*sph*sth - l*cph^2*sth*cth)*dgamma + mc*l^2*cph*cth^2*dphi + mc*l*cth*da;
c43 = mc*l^2*(cph*cth^2*dgamma + sth*cth*dphi);
c41 = mc*l*cth*dgamma;
c33 = bphi - mc*l^2*sth*cth*dtheta;
c44 = btheta;
C11 = [c11 c12; c21 c22]; C12 = [c13 c14; c23 c24]; C21 = [0 c32; c41 c42];
C22 = [c33 c34; c43 c44]; 
C = [C11 C12;C21 C22];

%--------------------------------------------------------------------------------------
g1 = mc*g*l*cth*sph;g2 = mc*g*l*sth*cph;
G2 = [g1; g2];G = [0;0;g1;g2];
%--------------------------------------------------------------------------------------
M1222 = M12/(M22);M2111 = M21/(M11);
%--------------------------------------------------------------------------------------
M1_ = M11 - M1222*M21;C11_ = C11 - M1222*C21;
C12_ = C12 - M1222*C22;G1 = -M1222*G2;
M2_ = M22 - M2111*M12;C21_ = C21 - M2111*C11;
C22_ = C22 - M2111*C12;
%--------------------------------------------------------------------------------------
invM1_ = inv(M1_);invM2_ = inv(M2_);
%--------------------------------------------------------------------------------------

%------------------------------------- Exact value of Fxp--------------------------------------------------
Fxp1 = invM1_*(-C11_*dq1 - C12_*dq2 - G1);
Fxp2 = invM2_*(-C21_*dq1 - C22_*dq2 - G2);
Fxp = Fxp1 + alphau*Fxp2;


sys(1)=x(1);
sys(2)=x(2);
sys(3)=x(3);
sys(4)=x(4);
sys(5)=x(5);
sys(6)=x(6);
sys(7)=x(7);
sys(8)=x(8);
sys(9) = Fxp(1);
sys(10) = Fxp(2);
end









