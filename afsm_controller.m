function [sys, x0, str, ts] = afsm_controller(t,x,u,flag)
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
    sizes.NumOutputs = 2;
    sizes.NumInputs =12;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 0;
    sys = simsizes(sizes);
    x0 = [];
    str = [];
    ts = [];
end

function sys = mdlOutputs(t,x,u)

 %------------------------------ Silding function and its derivative ----------------------
  s = [u(1); u(2)];
  ds = [u(3); u(4)];
  
 %-------------------------------F , G Matrix----------------------------------------------------------
Fxp = [u(5) ; u(6)];
G = [u(7) u(8);u(9) u(10)];
ddqr = [u(11) ; u(12)];
%--------------------------------- COntroller----------------------------------------------------
  
  Ke = diag([2,2]);
  ep  = 0.1;
if abs(s(1)/ep) >1
    s(1) = sign(s(1));
else 
    s(1) = s(1)/ep;
end
if abs(s(2)/ep) >1
 s(2) = sign(s(2));
else 
 s(2) = s(2)/ep;
end

 us  = -Ke*s;
 u = inv(G)*(-Fxp+ddqr+us);

%-----------------------------------------------------------------------  
  sys(1) = u(1);
  sys(2) = u(2);


end
