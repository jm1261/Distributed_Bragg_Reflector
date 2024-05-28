% Multilayer Bragg Reflector Program %
% Version 2, created 07/02/2020 %

% Clear all previous data before setup %
clear all;

%mkdir ("Results_"+date)
stamp="Remake_TFK_Results_200207";
mkdir ("Results_"+stamp);
results_path = pwd+"/Results_"+stamp;
set(0,'DefaultFigureVisible','off')

% Set Parameters %
n_cladding = 1.0;
n_substrate = 3.6;
n_effective = 0.0;

% Now we set up the mirrors, these are layers (x) of thickness (nm) %
for x=1:2:60
    n(x)=3.6;
    t(x)=54;
    n(x+1)=3.34;
    t(x+1)=58;
end

% Set up cavity thickness and index %
n_cavity=3.65;
target_lambda=780;
lambda_n = target_lambda/n_cavity;
cavity_thickness=round(6*(lambda_n/2));

% Set extinction coefficient for cavity layer, 5% AlGaAs %
pumped=0.0;
unpumped=0.05;

% Cavity thickness scalar %
d=0;
%% Pumped State %%
% Set cavity thickness d (nm) at some cavity x=cavity %
cavity=7;
t(cavity)=cavity_thickness+(d*cavity_thickness);
k=pumped;
n(cavity)=n_cavity-k*i; % Positive imaginary part is gain; k=0.05?1e4/cm %
n(x+2)=sqrt(n_substrate); % Lambda/4 layer %
t(x+2)=780/4/n(x+2);
max=length(n);

% Loop through the wavelength and calculate the matrix at each interface %
for m=1:1500
    lambda=700+0.1*m;
    wavelength=lambda;
    lam(m)=lambda;
    k0=2*pi/lambda;
    M0=[1 0;0 1];
   
    % Calculate the matrix at each layer %
    for x=1:max
        test=n(x)*n(x)-n_effective*n_effective;
        kappa=k0*sqrt(test);
        arg=kappa*t(x);
        M=[cos(arg) sin(arg)/kappa*i;sin(arg)*kappa*i cos(arg)];
        M1=M0*M;
        M0=M1;
    end
    
    % Using the matrix, calculate the transmission and reflection %
    kc=k0*sqrt(n_cladding*n_cladding-n_effective*n_effective);
    ks=k0*sqrt(n_substrate*n_substrate-n_effective*n_effective);
    transmission=2*ks/(ks*M0(1,1)+kc*M0(2,2)+ks*kc*M0(1,2)+M0(2,1));
    reflection=(ks*M0(1,1)-kc*M0(2,2)+ks*kc*M0(1,2)-M0(2,1))*transmission/(2*ks);
    phase=atan(real(transmission)/imag(transmission));
    %T1=transmission*conj(transmission);
    %R1=reflection*conj(reflection);
    T(m)=transmission*conj(transmission);
    R(m)=reflection*conj(reflection);
    P(m)=phase;
end

fig = plot(lam,R);
title("Reflection Spectrum - pumped");
%xlim([700 850]);
ylim([0 1.5]);
legend_string="thickness = "+t(cavity)+"nm";
legend(legend_string);
xlabel("Wavelength [nm]");
ylabel("Reflection [au]");
set(gca,'FontSize',16);
% save out
file_name="n"+n_cavity+"_x"+cavity+"_t"+t(cavity)+"_k"+k+".png";
file_path=results_path+"/"+file_name;
saveas(fig, file_path);

%% Unpumped State %%
% Set cavity thickness d (nm) at some cavity x=cavity %
cavity=7;
t(cavity)=cavity_thickness+(d*cavity_thickness);
k=unpumped;
n(cavity)=n_cavity-k*i; % Positive imaginary part is gain; k=0.05?1e4/cm %
n(x+2)=sqrt(n_substrate); % Lambda/4 layer %
t(x+2)=780/4/n(x+2);
max=length(n);

% Loop through the wavelength and calculate the matrix at each interface %
for m=1:1500
    lambda=700+0.1*m;
    wavelength=lambda;
    lam(m)=lambda;
    k0=2*pi/lambda;
    M0=[1 0;0 1];
   
    % Calculate the matrix at each layer %
    for x=1:max
        test=n(x)*n(x)-n_effective*n_effective;
        kappa=k0*sqrt(test);
        arg=kappa*t(x);
        M=[cos(arg) sin(arg)/kappa*i;sin(arg)*kappa*i cos(arg)];
        M1=M0*M;
        M0=M1;
    end
    
    % Using the matrix, calculate the transmission and reflection %
    kc=k0*sqrt(n_cladding*n_cladding-n_effective*n_effective);
    ks=k0*sqrt(n_substrate*n_substrate-n_effective*n_effective);
    transmission=2*ks/(ks*M0(1,1)+kc*M0(2,2)+ks*kc*M0(1,2)+M0(2,1));
    reflection=(ks*M0(1,1)-kc*M0(2,2)+ks*kc*M0(1,2)-M0(2,1))*transmission/(2*ks);
    phase=atan(real(transmission)/imag(transmission));
    %T1=transmission*conj(transmission);
    %R1=reflection*conj(reflection);
    T(m)=transmission*conj(transmission);
    R(m)=reflection*conj(reflection);
    P(m)=phase;
end

fig = plot(lam,R);
title("Reflection Spectrum - unpumped");
%xlim([700 850]);
ylim([0 1.5]);
legend_string="thickness = "+t(cavity)+"nm";
legend(legend_string);
xlabel("Wavelength [nm]");
ylabel("Reflection [au]");
set(gca,'FontSize',16);
% save out
file_name="n"+n_cavity+"_x"+cavity+"_t"+t(cavity)+"_k"+k+".png";
file_path=results_path+"/"+file_name;
saveas(fig, file_path);