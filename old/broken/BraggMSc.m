% Multilayer Bragg Reflector Program %% Clear all previous data before setup %clear all;% Set the refractive indices of the surrounding layers %% These layers are the cladding, substrate and the effective index %n_clad = 1.0;n_sub = 1.0;n_eff = 0.0;% Set the number of mirror pairsnumber_pairs = 28;% Set the high and low refractive indicesn_high = 3.6;n_low = 3.34;t_high = 54;t_low = 58;% Cavity wavelengthcav_lam = 750;% Now we set up the mirrors, these are layers (x) of thickness (nm)for x=1:2:number_pairs    n(x)=n_high;	t(x)=t_high;	n(x+1)=n_low;	t(x+1)=t_low;end%for p=-100:1:100    % Set the k of the complex refractive index%    k = 0.01*p    for d=1000    % Set the cavity to be at x=c, where i is some integer.    t(7)=d;    a=3.65;    k = 0.00;    n(7)=a-k*i; % Positive imaginary part is gain    %n(k+2)=sqrt(nsub);% Lambda/4 layer    %t(k+2)=650/4/n(k+2);    max=length(n);    % Scroll through the wavelength and calculate the matrix at each interface    for m=1:4000        lambda=600+0.1*m;        wavel=lambda;        lam(m)=lambda;        k0=2*pi/lambda;        M0=[1 0;0 1];        for x=1:max            test=n(x)*n(x)-n_eff*n_eff;            kappa=k0*sqrt(test);            arg=kappa*t(x);            M=[cos(arg) sin(arg)/kappa*i;sin(arg)*kappa*i cos(arg)];             M1=M0*M;            M0=M1;        end        % Using the matrix, calculate the transmission and reflection        kc=k0*sqrt(n_clad*n_clad-n_eff*n_eff);         ks=k0*sqrt(n_sub*n_sub-n_eff*n_eff);	        trans=2*ks/(ks*M0(1,1)+kc*M0(2,2)+ks*kc*M0(1,2)+M0(2,1));        ref=(ks*M0(1,1)-kc*M0(2,2)+ks*kc*M0(1,2)-M0(2,1))*trans/(2*ks);        phase=atan(real(trans)/imag(trans));        %T1=trans*conj(trans);        %R1=ref*conj(ref);        T(m)=trans*conj(trans);        R(m)=ref*conj(ref);        P(m)=phase;    end    %semilogy (lam,T);     %plot (lam, R, lam, P);    %ylim([-0.1 1.1]);    %plot (lam, P);    %ylim([-2 2]);    fig = plot (lam,R);    title('Reflection Spectrum - unpumped')    % legend("thickness = 641nm")    grid on;    ylim([0 1]);    xlim([600 1000]);    set(gca,'FontSize',20)    % save out    %fpath = '/media/mass_storage/josh/Documents/PhD_Project/Optalysys/Bragg_Mirror/Results_200204';    %fname =  "n"+a+"_x16"+"_t"+t(16)+"_k"+k+".png"    %saveas(fig, fullfile(fpath, fname));end