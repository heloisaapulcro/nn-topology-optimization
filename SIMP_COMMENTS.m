%%%% A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLE SIGMUND, JANUARY 2000 %%%
%%%% CODE MODIFIED FOR INCREASED SPEED, September 2002, BY OLE SIGMUND %%%
function SIMP(nelx,nely,volfrac,penal,rmin,newF,problem) 

%%%%%%%%%% CANTILEVER: SIMP(64,40,0.40,3,1.5,-1,1)  %%%%%% FORÇA NO PONTO MÉDIO 
%%%%%%%%%% CANTILEVER: SIMP(64,40,0.40,3,1.5,-1,2)  %%%%%% FORÇA NA EXTREMIDADE SUPERIOR
%%%%%%%%%% CANTILEVER: SIMP(64,40,0.40,3,1.5,-1,3)  %%%%%% FORÇA NA EXTREMIDADE INFERIOR
%%%%%%%%%% VIGA - MBB: SIMP(120,20,0.40,3,1.5,-1,4) %%%%%% VIGA BIAPOIADA COM FORÇA NA EXTREMIDADE SUPERIOR
%%%%%%%%%% MICHELL: SIMP(100,50,0.40,3,1.5,-1,5)    %%%%%% VIGA BIAPOIADA COM FORÇA NA EXTREMIDADE INFERIOR
tic
% INITIALIZE
x(1:nely,1:nelx) = volfrac; 
%%% definição de uma malha homogenea do elemento
loop = 0; 
change = 1.;
i=1;
% START ITERATION
while change > 0.01
  %% mede em relação ao volfrac as alterações maior que 1%
  i=i+1;
  loop = loop + 1;
  xold = x;

% FE-ANALYSIS
  [U]=FE(nelx,nely,x,penal,newF,problem);
  
% FUNCAO OBJETIVO E ANALISE DE SENSIBILIDADE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [KE] = lk;
  c(i) = 0.;
  vol(i)=0.;
  for ely = 1:nely
    for elx = 1:nelx
      n1 = (nely+1)*(elx-1)+ely; 
      n2 = (nely+1)* elx   +ely;
      Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
      c(i) = c(i) + x(ely,elx)^penal*Ue'*KE*Ue;
      dc(ely,elx) = -penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue;
    end
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % CHAMADA PARA FILTRO DE SENSIBILIDADE
 [dc]   = check(nelx,nely,rmin,x,dc);  
 vol(i)=sum(sum(x))/(nelx*nely);
%%%%%%%% ATUALIZA O PROJETO PELO CRITERIO DE OTIMALIDADE %%%%%%%%%%%%%%%%%%%%
 [x]    = OC(nelx,nely,x,volfrac,dc);
%%%%%%%%%%%%%%%%%%%%IMPRIME OS RESULTADOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  change = max(max(abs(x-xold)));
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c(i)) ...
       ' Vol.: ' sprintf('%6.3f',vol(i)) ...
        ' ch.: ' sprintf('%6.3f',change )])
%%%% BESO A ESTRUTURA e PLOTA GRAFICO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(1); contourf(x,[0,0]);
        colormap(gray);imagesc(-x); axis equal;  axis tight; axis off; pause(1e-6);

%%%%%%%%%%%% GERA O GIF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256); % Converte para formato de GIF

        gif_filename = 'cantilever3.gif'; % Mudar nome do arquivo para não sobrepor os outros

        if loop == 1
            imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.2); % Primeiro frame
        else
            imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.2);
        end
        figure(2); surf(x); caxis([-12,12]); 
        axis equal; axis([0,nelx,0,nely,-12,12]); view(3);
        figure(3); subplot(2,1,1); plot(c(1:i),'-'); title('Compliance');
                   subplot(2,1,2); plot(vol(1:i),'-'); title('Volume fraction');
end 
toc
%%%%%%%%%% FUNCAO (CRITERIO DE OTIMALIDADE)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnew]=OC(nelx,nely,x,volfrac,dc)  
l1 = 0; l2 = 100000; move = 0.2;
while (l2-l1 > 1e-4)
  lmid = 0.5*(l2+l1);
  xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(-dc./lmid)))));
  if sum(sum(xnew)) - volfrac*nelx*nely > 0;
    l1 = lmid;
  else
    l2 = lmid;
  end
end
%%%%%%%%%% FUNCAO (FILTRO DE INDEPENDENCIA DE MALHA) %%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn]=check(nelx,nely,rmin,x,dc)
dcn=zeros(nely,nelx);
for i = 1:nelx
  for j = 1:nely
    sum=0.0; 
    for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
      for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
        fac = rmin-sqrt((i-k)^2+(j-l)^2);
        sum = sum+max(0,fac);
        dcn(j,i) = dcn(j,i) + max(0,fac)*x(l,k)*dc(l,k);
      end
    end
    dcn(j,i) = dcn(j,i)/(x(j,i)*sum);
  end
end
%%%%%%%%%% FUNCAO (ANALISE POR ELEMENTOS FINITOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U]=FE(nelx,nely,x,penal,newF,problem)
%% define a malha de elementos finitos e as condicoes de contorno
[KE] = lk; 
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
F = sparse(2*(nely+1)*(nelx+1),1); U = zeros(2*(nely+1)*(nelx+1),1);
for elx = 1:nelx
  for ely = 1:nely
    n1 = (nely+1)*(elx-1)+ely; 
    n2 = (nely+1)* elx   +ely;
    edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;
  end
end
%%%%%%%%%%%%%%%%%%%% CANTILEVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if problem == 1
        Name='CANTILEVER';
        F(2*(nelx+1)*(nely+1)-nely,1) = newF; %FORÇA NO PONTO MÉDIO
        fixeddofs   = (1:2*(nely+1)); 
        alldofs     = (1:2*(nely+1)*(nelx+1));
        freedofs    = setdiff(alldofs,fixeddofs);
        % SOLVING
        U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
        U(fixeddofs,:)= 0;
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%FIM DE CANTILEVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% MÃO FRANCESA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if problem == 2
        Name='MÃO FRANCESA INVERTIDA';
        F(2*(nelx)*(nely+1)+2,1) = newF;  %FORÇA NA EXTREMIDADE SUPERIOR
        fixeddofs   = (1:2*(nely+1)); 
        alldofs     = (1:2*(nely+1)*(nelx+1));
        freedofs    = setdiff(alldofs,fixeddofs);
        % SOLVING
        U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
        U(fixeddofs,:)= 0;
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%FIM DE MAO FRANCESA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% MÃO FRANCESA INVERTIDA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if problem == 3
        Name='MÃO FRANCESA';
        F(2*(nelx+1)*(nely+1),1) = newF;  %FORÇA NA EXTREMIDADE INFERIOR
        fixeddofs   = (1:2*(nely+1)); 
        alldofs     = (1:2*(nely+1)*(nelx+1));
        freedofs    = setdiff(alldofs,fixeddofs);
        % SOLVING
        U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
        U(fixeddofs,:)= 0;
   end
   %%%%%%%%%%%%%%%%% VIGA MBB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if problem == 4
        Name='VIGA MBB';
        F(2*(nelx+1)*(nely+1)-nely,1) = newF; %FORÇA NO PONTO MÉDIO
        fixeddofs   = (1:2*(nely+1)); 
        alldofs     = (1:2*(nely+1)*(nelx+1));
        freedofs    = setdiff(alldofs,fixeddofs);
        % SOLVING
        U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
        U(fixeddofs,:)= 0;
   end
   %%%%%%%%%%%%%%%%%%%%CANTILEVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if problem == 5
        Name='MICHELL';
        F(2*(nely+1)*(nelx/2+1),1) = -1.;
        fixeddofs = union(2*nely+1:2*(nely+1), 2*(nely+1)*(nelx+1));
        alldofs = (1:2*(nely+1)*(nelx+1));
        freedofs = setdiff(alldofs,fixeddofs);
        % SOLVING
        U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
        U(fixeddofs,:)= 0;
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%FIM DE MAO FRANCESA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% MATRIZ DE RIGIDEZ DO ELEMENTO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE]=lk
E = 1.; 
nu = 0.3;
%% mudar o coeficiente de poisson?
k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ... 
   -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];
KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                  k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                  k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                  k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                  k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                  k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                  k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                  k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
%% rigidez de um elemento é k=E/1-nu² ajusta os coeficientes com as propriedades do material
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by Ole Sigmund, Department of Solid         %
% Mechanics, Technical University of Denmark, DK-2800 Lyngby, Denmark.     %
% Please sent your comments to the author: sigmund@fam.dtu.dk              %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "A 99 line topology optimization code written in Matlab"                 %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The author reserves all rights but does not guaranty that the code is    %
% free from errors. Furthermore, he shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
