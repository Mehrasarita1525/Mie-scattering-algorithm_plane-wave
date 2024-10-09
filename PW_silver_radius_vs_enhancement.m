clc
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plane Wave %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_m=1.52;                   %optical index of surrounding medium(glass/air)
 
N=15;                     %order of multipole


C1=0;
for a=1:1:100                %Radius of nanoparticle
    C1=C1+1;
    C2=0;
for Inter_lambda=300:1:900
    C2=C2+1;

    lambda=[301 311 320 331 342 354 368 381 397 413 430 451 471 496 521 548 582 616 659 704 755 821 891];     %wavelength of incident light(nm)  
    R=[1.34 1.13 0.81 0.17 0.14 0.10 0.07 0.05 0.05 0.05 0.04 0.04 0.05 0.05 0.05 0.06 0.05 0.06 0.05 0.04 0.03 0.04 0.04]; 
    I=[0.964 0.616 0.392 0.829 1.142 1.419 1.657 1.864 2.070 2.275 2.462 2.657 2.869 3.093 3.324 3.586 3.858 4.152 4.483 4.838 5.242 5.727 6.312];
    
    Agr=interp1(lambda,R,Inter_lambda);
    Agi=interp1(lambda,I,Inter_lambda);
    
    KrP=0;
   
% Input Parameters  
   n_Ag=Agr+(1i*Agi);          %optical index of particle(gold/silver)
   k=2*pi*n_m/Inter_lambda;    %wavenumber
   m=n_Ag/n_m;                 %m=relative refractive index
   x=k*a;                      %x=size parameter
   z=m*x;
  
      
for n=1:N
   sx=sqrt(pi*x/2);
   sz=sqrt(pi*z/2);

   Phinx=sx.*besselj(n+0.5,x);                 
   Phin1x=sx.*besselj(n-0.5,x);   
   Phinz=sz.*besselj(n+0.5,z);                       
   Phin1z=sz.*besselj(n-0.5,z);
     
   Xinx=sx.*((besselj(n+0.5,x))+1i*(bessely(n+0.5,x)));       
   Xin1x=sx.*((besselj(n-0.5,x))+1i*(bessely(n-0.5,x))); 
      
   Phindx=Phin1x-(n.*Phinx./x);
   Phindz=Phin1z-(n.*Phinz./z);
   Xindx=Xin1x-(n.*(Xinx)./x);
   
an=((m*Phinz.*Phindx)-(Phinx.*Phindz))./((m*Phinz.*Xindx)-(Xinx.*Phindz));
bn=((Phinz.*Phindx)-(m*Phinx.*Phindz))./((Phinz.*Xindx)-(m*Xinx.*Phindz));

W=((abs(an))^2.*(abs(Xinx))^2)+(abs(Phinx)^2)-(2*real(an.*Xinx).*Phinx);
X=((abs(an)^2).*(abs(Xindx)^2))+((abs(bn)^2).*(abs(Xinx)^2))+(abs(Phinx)^2)+(abs(Phindx)^2)-(2*real(an.*Xindx).*Phindx)-(2*real(bn.*Xinx).*Phinx);
Y=(2*(n^3)*(n+1)^3)./((2*n)+1);
Z=(2*(n^2)*(n+1)^2)./((2*n)+1);

    E1=3/2*1i;
    En=1i^n.*(((2*n)+1)./(n.*(n+1)));
    EnPW=(abs(En)^2)./(abs(E1)^2);
    
    Kr=EnPW.*Y.*W;
    KrP=Kr+KrP;
    
end
    KrPW=9.*(KrP)./(16.*x^4);
      
d(C1,C2)=KrPW;

end
end

a=1:1:100;
Inter_lambda=300:1:900;
figure(1);imagesc(Inter_lambda,a,d);
colormap(hot)
xlabel('Wavelength(nm)')
ylabel('Radius of nanoparticle(nm)')


