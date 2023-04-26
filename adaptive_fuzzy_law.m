function [sys, x0, str, ts] = adaptive_fuzzy_law(t,x,u,flag)
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
    sizes.NumContStates = 2592 ;
    sizes.NumDiscStates = 0;
    sizes.NumOutputs = 2;
    sizes.NumInputs =4 ;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 0;
    sys = simsizes(sizes);
    x0 = [1*ones(1296,1); 1*ones(1296,1)];
    str = [];
    ts = [];
end

function sys = mdlDerivatives(t,x,u)
% ------------------------- Parameters Setup ------------------------------
  gammaf = 1000; 
  thetaf_up = [200;200] ; thetaf_down = [-200;-200];
  
%------------------------Sliding function and its derivative -------------
  s = [u(1); u(2)];
  ds = [u(3); u(4)];
  
%--------------------------------------Fuzzy System-----------------------------------------------
     u11 = exp(-1/2*((s+2.5)/0.5).^2); 
     u12 = exp(-1/2*((s+1.5)/0.5).^2);
     u13 = exp(-1/2*((s+0.5)/0.5).^2);
     u14 = exp(-1/2*((s-0.5)/0.5).^2);
     u15 = exp(-1/2*((s-1.5)/0.5).^2);
      u16 = exp(-1/2*((s-2.5)/0.5).^2); 

     u21 = exp(-1/2*((ds+2.5)/0.5).^2); 
     u22 = exp(-1/2*((ds+1.5)/0.5).^2);
     u23 = exp(-1/2*((ds+0.5)/0.5).^2);
     u24 = exp(-1/2*((ds-0.5)/0.5).^2);
     u25 = exp(-1/2*((ds-1.5)/0.5).^2);
     u26 = exp(-1/2*((ds-2.5)/0.5).^2); 
   u1 = [u11,u12,u13,u14,u15,u16]; u2 = [u21,u22,u23,u24,u25,u26];
  
  for i = 1:1:1296
      thetaf(i,1) = x(i); 
      thetaf(i+1296,1) = x(i+1296);
  end

  FS1 = 0;
  for l11 = 1:1:6
      for l12 = 1:1:6
          for l21 = 1:1:6
              for l22 = 1:1:6
                  idx = 216*(l11-1)+36*(l12-1)+6*(l21-1)+l22;
                  FS2(idx) = u1(1,l11)*u1(2,l12)*u2(1,l21)*u2(2,l22);
                  FS1 = FS1 + u1(1,l11)*u1(2,l12)*u2(1,l21)*u2(2,l22);
              end
          end
      end
  end
  FS = FS2/FS1;
  W = [FS , zeros(1,1296);zeros(1,1296), FS];

%--------------------------------- Adaptive law --------------------------  
  
  tempf = gammaf*W'*s;

  for i = 1:1:1296
      sys(i) = project(tempf(i), thetaf(i,1),thetaf_down(1), thetaf_up(1));
      sys(i+1296) = project(tempf(i+1296),thetaf(i+1296),thetaf_down(2), thetaf_up(2));
  end
  

end

function sys = mdlOutputs(t,x,u)
 %------------------------------ Silding function and its derivative ----------------------
  s = [u(1); u(2)];
  ds = [u(3); u(4)];
  
  
%--------------------------------------Fuzzy System----------------------
     u11 = exp(-1/2*((s+2.5)/0.5).^2); 
     u12 = exp(-1/2*((s+1.5)/0.5).^2);
     u13 = exp(-1/2*((s+0.5)/0.5).^2);
     u14 = exp(-1/2*((s-0.5)/0.5).^2);
     u15 = exp(-1/2*((s-1.5)/0.5).^2);
     u16 = exp(-1/2*((s-2.5)/0.5).^2); 

     u21 = exp(-1/2*((ds+2.5)/0.5).^2); 
     u22 = exp(-1/2*((ds+1.5)/0.5).^2);
     u23 = exp(-1/2*((ds+0.5)/0.5).^2);
     u24 = exp(-1/2*((ds-0.5)/0.5).^2);
     u25 = exp(-1/2*((ds-1.5)/0.5).^2);
     u26 = exp(-1/2*((ds-2.5)/0.5).^2); 
   u1 = [u11,u12,u13,u14,u15,u16]; u2 = [u21,u22,u23,u24,u25,u26];
  
  for i = 1:1:1296
      thetaf(i,1) = x(i); 
      thetaf(i+1296,1) = x(i+1296);
  end

  FS1 = 0;
  for l11 = 1:1:6
      for l12 = 1:1:6
          for l21 = 1:1:6
              for l22 = 1:1:6
                  idx = 216*(l11-1)+36*(l12-1)+6*(l21-1)+l22;
                  FS2(idx) = u1(1,l11)*u1(2,l12)*u2(1,l21)*u2(2,l22);
                  FS1 = FS1 + u1(1,l11)*u1(2,l12)*u2(1,l21)*u2(2,l22);
              end
          end
      end
  end
  FS = FS2/FS1;
  W = [FS , zeros(1,1296);zeros(1,1296), FS];
  
  Fxp = W*thetaf;

%-----------------------------------------------------------------------  
  sys(1) = Fxp(1);
  sys(2) = Fxp(2);

end

function p = project(x,y,a,b)
 if (y <= a) && ( x<0)
     p = 0;
 elseif (y>=b) && (x>0)
     p = 0;
 else
     p = x; 
 end
end