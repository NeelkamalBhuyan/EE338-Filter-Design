% Parameters for Chebyshev Analog Low Pass Filter

% order of filter
Nmin = 4;
D1 = 1/(0.85*0.85)-1;      
epsilon = sqrt(D1);

% Poles of Chebyshev analog LPF (4 roots with Re(s)<0 and 4 roots with
% Re(s)>0
root1 = -sin(pi/(2*Nmin))*sinh(asinh(1/epsilon)/Nmin)+1i*cos(pi/(2*Nmin))*cosh(asinh(1/epsilon)/Nmin);
root2 = -sin(pi/(2*Nmin))*sinh(asinh(1/epsilon)/Nmin)-1i*cos(pi/(2*Nmin))*cosh(asinh(1/epsilon)/Nmin);
root3 = -sin(3*pi/(2*Nmin))*sinh(asinh(1/epsilon)/Nmin)+1i*cos(3*pi/(2*Nmin))*cosh(asinh(1/epsilon)/Nmin);
root4 = -sin(3*pi/(2*Nmin))*sinh(asinh(1/epsilon)/Nmin)-1i*cos(3*pi/(2*Nmin))*cosh(asinh(1/epsilon)/Nmin);        

% Transfer function with  poles as root1 to root4
polynomial1 = [1 -root1-root2 root1*root2];
polynomial2 = [1 -root3-root4 root3*root4];
denominator = conv(polynomial1,polynomial2);          
numerator = [denominator(5)*sqrt(1/(1+epsilon*epsilon))];  

% Un-normalised speifications of Bandstop filter
SB1 = 31;
PB1 = 27;
PB2 = 55;
SB2 = 51;

% Normalised specifications of Bandstop filter
F_sampling = 260;
omega_s1 = tan(SB1/F_sampling*pi);          
omega_p1 = tan(PB1/F_sampling*pi);
omega_p2 = tan(PB2/F_sampling*pi);
omega_s2 = tan(SB2/F_sampling*pi);

% Omega_0 and B for frequency Transformation
omega_0 = sqrt(omega_p1*omega_p2);
B = omega_p2-omega_p1;

% Frequency response of analog lpf, then analog bandstop filter, then
% discrete time bandstop filter
syms s z;

% analog lowpass filter Transfer Function
analog_lpf(s) = poly2sym(numerator,s)/poly2sym(denominator,s);

% analog bandstop filter Transfer Function
analog_bsf(s) = analog_lpf((B*s)/(s*s + omega_0*omega_0));

% discrete bandstop filter Transfer Function
discrete_bsf(z) = analog_bsf((z-1)/(z+1));          

% coefficients of analog bandstop filter
[num_s, den_s] = numden(analog_bsf(s));                   
num_s = sym2poly(expand(num_s));                          
den_s = sym2poly(expand(den_s));                          
DC = den_s(1);
% normalising the coefficients with highest power of denominator = 1
den_s = den_s/DC;
num_s = num_s/DC;

%coefficients of discrete bandstop filter
[num_z, den_z] = numden(discrete_bsf(z));                 
num_z = sym2poly(expand(num_z));                          
den_z = sym2poly(expand(den_z));                          
DC = den_z(1);     
% normalising the coefficients with highest power of denominator = 1
den_z = den_z/DC;
num_z = num_z/DC;

% visualising frequency response
fvtool(num_z,den_z)                                       

% magnitude plot
[H,f] = freqz(num_z,den_z,1024*1024, 260e3);
plot(f,abs(H))
grid