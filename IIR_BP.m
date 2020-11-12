% Parameters for Butterworth Analog Low Pass Filter

%cut-off frequency
omega_cutoff = sqrt(1.054599*1.089411);    
% order of filter
N = 9;                   

% Poles of Butterworth analog LPF (9 roots with Re(s)<0 and 9 roots with
% Re(s)>0
root1 = omega_cutoff*cos(pi/2 + pi/18) + 1i*omega_cutoff*sin(pi/2 + pi/18);
root2 = omega_cutoff*cos(pi/2 + pi/18) - 1i*omega_cutoff*sin(pi/2 + pi/18);
root3 = omega_cutoff*cos(pi/2 + pi/18+pi/9) + 1i*omega_cutoff*sin(pi/2 + pi/18+pi/9);
root4 = omega_cutoff*cos(pi/2 + pi/18+pi/9) - 1i*omega_cutoff*sin(pi/2 + pi/18+pi/9);
root5 = omega_cutoff*cos(pi/2 + pi/18+2*pi/9) + 1i*omega_cutoff*sin(pi/2 + pi/18+2*pi/9);
root6 = omega_cutoff*cos(pi/2 + pi/18+2*pi/9) - 1i*omega_cutoff*sin(pi/2 + pi/18+2*pi/9);
root7 = omega_cutoff*cos(pi/2 + pi/18+3*pi/9) + 1i*omega_cutoff*sin(pi/2 + pi/18+3*pi/9);
root8 = omega_cutoff*cos(pi/2 + pi/18+3*pi/9) - 1i*omega_cutoff*sin(pi/2 + pi/18+3*pi/9);
root9 = omega_cutoff*cos(pi);


% Un-normalised speifications of Bandpass filter
PB1 = 32.8;
SB1 = 28.8;
SB2 = 56.8;
PB2 = 52.8;

% Normalised specifications of Bandpass filter
F_sampling = 330;         
omega_p1 = tan(PB1/F_sampling*pi);
omega_s1 = tan(SB1/F_sampling*pi); 
omega_s2 = tan(SB2/F_sampling*pi);
omega_p2 = tan(PB2/F_sampling*pi);

% Omega_0 and B for frequency Transformation
omega_0 = sqrt(omega_p1*omega_p2);
B = omega_p2-omega_p1;

% Transfer function with no zeroes, numerator as omega_cutoff^N (so that DC value is 1)and poles
% as root1 to root9
[numerator,denenominator] = zp2tf([],[root1 root2 root3 root4 root5 root6 root7 root8 root9],omega_cutoff^N);   

% Frequency response of analog lpf, then analog bandpass filter, then
% discrte time bandpass filter
syms s z;

% analog lowpass filter Transfer Function
analog_lpf(s) = poly2sym(numerator,s)/poly2sym(denenominator,s);

% analog bandpass filter Transfer Function
analog_bpf(s) = analog_lpf((s*s+omega_0*omega_0)/(B*s));

% discrete bandpass filter Transfer Function
discrete_bpf(z) = analog_bpf((z-1)/(z+1));

% coefficients of analog bandpass filter
[num_s, den_s] = numden(analog_bpf(s));                   
num_s = sym2poly(expand(num_s));                          
den_s = sym2poly(expand(den_s));                          
DC = den_s(1);
% normalising the coefficients with highest power of denominator = 1
den_s = den_s/DC;
num_s = num_s/DC;

%coefficients of discrete bandpass filter
[num_z, den_z] = numden(discrete_bpf(z));                                         
num_z = sym2poly(expand(num_z));
den_z = sym2poly(expand(den_z));                              
DC = den_z(1);
% normalising the coefficients with highest power of denominator = 1
den_z = den_z/DC;
num_z = num_z/DC;

% visualising frequency response
fvtool(num_z,den_z)                                           

% magnitude plot
[H,f] = freqz(num_z,den_z,1024*1024, 330e3);
plot(f,abs(H))
grid