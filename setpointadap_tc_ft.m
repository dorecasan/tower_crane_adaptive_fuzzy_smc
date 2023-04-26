function [sys,x0,str,ts] = setpointadap_tc_ft(t,x,u,flag)
switch flag
case 0
[sys,x0,str,ts]=mdlInitializeSizes;
case 1
sys=mdlDerivatives(t,x,u);
case 3
sys=mdlOutputs(t,x,u);
case {2,4,9}
sys=[];
otherwise
error(['Unhandled flag = ',num2str(flag)]);
end
end
function [sys,x0,str,ts]=mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates = 0;
sizes.NumDiscStates = 0;
sizes.NumOutputs = 6;
sizes.NumInputs = 0;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;
sys = simsizes(sizes);
x0 = [];
str = [];
ts = [0 0];
end
function sys = mdlOutputs(t,x,u)
 T = 5;sr0 = 0.001;ar0 = 0.001;st = 2;at = pi/4;
 if t < T
     xr = sr0+(35*(t/T)^4-84*(t/T)^5+70*(t/T)^6-20*(t/T)^7)*(st-sr0);
     phir = ar0+(35*(t/T)^4-84*(t/T)^5+70*(t/T)^6-20*(t/T)^7)*(at-ar0);
     dxr = (140*(t/T)^3-420*(t/T)^4+420*(t/T)^5-140*(t/T)^6)*(st-sr0)/T;
     dphir = (140*(t/T)^3-420*(t/T)^4+420*(t/T)^5-140*(t/T)^6)*(at-ar0)/T;
     ddxr = (420*(t/T)^2-1680*(t/T)^3+2100*(t/T)^4-840*(t/T)^5)*(st-sr0)/(T^2);
     ddphir = (420*(t/T)^2-1680*(t/T)^3+2100*(t/T)^4-840*(t/T)^5)*(at-ar0)/(T^2);
 else
     xr =st;
     dxr =0; ddxr=0;
     phir =at;dphir=0;ddphir=0;
 end
  sys(1) = xr; sys(2) = dxr; sys(3) = ddxr;
  sys(4) = phir; sys(5) = dphir; sys(6) = ddphir;
end









