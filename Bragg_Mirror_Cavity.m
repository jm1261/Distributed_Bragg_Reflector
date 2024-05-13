%% Multilayer Bragg Reflector Program %%

clear all; % Clear all previous data %

if 0==exist ("Results_"+date)
    mkdir ("Results_"+date)
end
results_path = "/media/mass_storage/josh/Documents/PhD_Project/Optalysys/Bragg_Mirror/Results_"+date;

% Set the refractive indices of surround layers: cladding, substrate and %
% effective index %
n_cladding=1.0;
n_substrate=1.0;
n_effective=0;

% Set extinction coefficient for cavity layer, 5% AlGaAs %
pumped=0.0;
unpumped=0.05;

% Set up the mirrors, these are layers (x) of thickness (nm) %
for x=1:2:28 % x=start:step:finish %
    n(x)=3.6; % Layer 1 is  15% AlGaAs %
    t(x)=54;
    n(x+1)=3.34; % Layer 2 is 50% AlGaAs %
    t(x+1)=58;
end
%% Pumped State %%

% Set cavity thickness d (nm) at some cavity x=cavity %
cavity=7; 
for d=100:10:1000 % d=start:step:end %
    t(cavity)=d;
    a=3.65;
    k=pumped;
    n(cavity)=a-(k*i); % Positive imaginary part is gain %
    max=length(n);
    
    lambda_i=600; % Set initial wavelength%
    for m=1:4000 % Loop wavelength %
        lambda=lambda_i+0.1*m;
        %wavelength=lambda;
        lam(m)=lambda;
        k0=2*pi/lambda; % Set k0 vector %
        M0=[1 0;0 1]; % Set initial matrix %
        for x=1:max
            % Do matrix calculation for each interface %
            test=n(x)*n(x)-n_effective*n_effective;
            kappa=k0*sqrt(test);
            arg=kappa*t(x);
            M=[cos(arg) sin(arg)/kappa*i;sin(arg)*kappa*i cos(arg)];
            M1=M0*M;
            % Update initial matrix %
            M0=M1;
        end
        
        % Using calculate matrix, determine transmission and reflection %
        kc=k0*sqrt(n_cladding*n_cladding-n_effective*n_effective); 
        ks=k0*sqrt(n_substrate*n_substrate-n_effective*n_effective);
        transmission=2*ks/(ks*M0(1,1)+kc*M0(2,2)+ks*kc*M0(1,2)+M0(2,1));
        reflection=(ks*M0(1,1)-kc*M0(2,2)+ks*kc*M0(1,2)-M0(2,1))*transmission/(2*ks);
        phase=atan(real(transmission)/imag(transmission));
        T(m)=transmission*conj(transmission);
        R(m)=reflection*conj(reflection);
        P(m)=phase;
    end
    
    fig = plot(lam,R);
    title("Reflection Spectrum - pumped");
    xlim([700 850]);
    ylim([0 1.5]);
    legend_string="thickness = "+t(cavity)+"nm";
    legend(legend_string);
    xlabel("Wavelength [nm]");
    ylabel("Reflection [au]");
    set(gca,'FontSize',16);
    % save out
    file_name =  "n"+a+"_x"+cavity+"_t"+t(cavity)+"_k"+k+".png";
    saveas(fig, fullfile(results_path, file_name));
end
%% Unpumped State %%

% Set cavity thickness d (nm) at some cavity x=cavity %
cavity=7; 
for d=100:10:1000 % d=start:step:end %
    t(cavity)=d;
    a=3.65;
    k=unpumped;
    n(cavity)=a-(k*i); % Positive imaginary part is gain %
    max=length(n);
    
    lambda_i=600; % Set initial wavelength%
    for m=1:4000 % Loop wavelength %
        lambda=lambda_i+0.1*m;
        %wavelength=lambda;
        lam(m)=lambda;
        k0=2*pi/lambda; % Set k0 vector %
        M0=[1 0;0 1]; % Set initial matrix %
        for x=1:max
            % Do matrix calculation for each interface %
            test=n(x)*n(x)-n_effective*n_effective;
            kappa=k0*sqrt(test);
            arg=kappa*t(x);
            M=[cos(arg) sin(arg)/kappa*i;sin(arg)*kappa*i cos(arg)];
            M1=M0*M;
            % Update initial matrix %
            M0=M1;
        end
        
        % Using calculate matrix, determine transmission and reflection %
        kc=k0*sqrt(n_cladding*n_cladding-n_effective*n_effective); 
        ks=k0*sqrt(n_substrate*n_substrate-n_effective*n_effective);
        transmission=2*ks/(ks*M0(1,1)+kc*M0(2,2)+ks*kc*M0(1,2)+M0(2,1));
        reflection=(ks*M0(1,1)-kc*M0(2,2)+ks*kc*M0(1,2)-M0(2,1))*transmission/(2*ks);
        phase=atan(real(transmission)/imag(transmission));
        T(m)=transmission*conj(transmission);
        R(m)=reflection*conj(reflection);
        P(m)=phase;
    end
    
    fig = plot(lam,R);
    title("Reflection Spectrum - pumped");
    xlim([700 850]);
    ylim([0 1.5]);
    legend_string="thickness = "+t(cavity)+"nm";
    legend(legend_string);
    xlabel("Wavelength [nm]");
    ylabel("Reflection [au]");
    set(gca,'FontSize',16);
    % save out
    file_name =  "n"+a+"_x"+cavity+"_t"+t(cavity)+"_k"+k+".png";
    saveas(fig, fullfile(results_path, file_name));
end