%%%% A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLE SIGMUND, JANUARY 2000 %%%
%%%% CODE MODIFIED FOR INCREASED SPEED, September 2002, BY OLE SIGMUND %%%

function SIMP(nelx,nely,volfrac,penal,rmin,newF,problem,geometry) 
%{
    nelx ⇾ number of elements in direction X
    nely ⇾ number of elements in direction Y
    volfrac ⇾ volume fraction (minimum volume in % to be achieved)
    penal ⇾ SIMP penalty factor (usually between 3 and 4) x^penal
    rmin ⇾ filtering radius to smooth sensitivity (avoids irregular patterns)
    newF ⇾ value and direction of the applied force
    problem ⇾ type of problem that define boundary conditions and loading
    geometry ⇾ beam geometry type

%   ======================== PROBLEMS DEFINITION ==========================
1. CANTILEVER:  FORCE AT MIDPOINT 
2. CANTILEVER:  FORCE AT THE UPPER END
3. CANTILEVER:  FORCE AT THE LOWER END
4. VIGA-MBB:    DOUBLY SUPPORTED BEAM WITH FORCE AT THE UPPER END
5. MICHELL:     DOUBLY SUPPORTED BEAM WITH FORCE AT THE LOWER END
%}

%   =========================== INITIALIZE ================================
tic % Timer
x(1:nely,1:nelx) = volfrac;
% Density matrix (all elements start with the value of volfrac)

loop = 0;
change = 1.; % Initial maximum change
i=1;

% Geometry definition
% R
if geometry == 1
    x(1:nely,1:nelx) = volfrac;
end

% L
if geometry == 2
    % Lower right rectangular void
    x(ceil(nely/2):nely, ceil(nelx/2):nelx) = 0.001;
end

% T
if geometry == 3
    % Lower left rectangular void
    x(ceil(nely/2):nely, 1:ceil(nelx/4)) = 0.001; 
    % Lower right rectangular void
    x(ceil(nely/2):nely, nelx-ceil(nelx/4):nelx) = 0.001; 
end

% U
if geometry == 4
    % Lower rectangular void
    x(ceil(nely/2):nely, ceil(nelx/4):nelx-ceil(nelx/4)) = 0.001;
end
mask=x;

% Save the results with parameters in name
outputName = sprintf('r_nx%d_ny%d_V%.2f_p%d_rm%.2f_F%.2f_P%d_G%d', ...
    nelx, nely, volfrac, penal, rmin, newF, problem, geometry);

% Save initial matrix
txtNameInitial = sprintf('%s_inicial.txt', outputName);
writematrix(x, txtNameInitial, 'Delimiter', ',');

gifName = sprintf('%s.gif', outputName);

%   =================== START ITERATION ⇾ OPTIMIZATION ===================
while change > 0.01 
  i=i+1;
  loop = loop + 1;
  xold = x; % Save previous configuration

%   =========================== FE-ANALYSIS ===============================
  [U]=FE(nelx,nely,x,penal,newF,problem); % Calculates U displacements
  

%   ============= OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS =============
  [KE] = lk; % Element stiffness matrix
  c(i) = 0.; % Objective function (compliance) < compliance = > stiffness
  vol(i)=0.;
  for ely = 1:nely
    for elx = 1:nelx
      % Element nodes in the mesh
      n1 = (nely+1)*(elx-1)+ely; 
      n2 = (nely+1)* elx   +ely;
      
      % Degress of freedom of the elements (2* = 2 degressof freedom)
      Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
      
      % Compliance sum
      c(i) = c(i) + x(ely,elx)^penal*Ue'*KE*Ue;
      
      % Sensitivity of the objective function in relation to density
      dc(ely,elx) = -penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue;
    end
  end


%   ========================= SENSITIVITY FILTER ==========================
 [dc]   = check(nelx,nely,rmin,x,dc);  
 
% Updates density by optimality criterion
 vol(i) = sum(sum(x))/(nelx*nely);
 [x] = OC(nelx,nely,x,volfrac,dc,mask);
 
 
%   ====================== PLOTS AND RESULTS - BESO ======================= 
% Convergence check
  change = max(max(abs(x-xold)));
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c(i)) ...
       ' Vol.: ' sprintf('%6.3f',vol(i)) ...
        ' ch.: ' sprintf('%6.3f',change )])

        % Binary contour
        figure(1); contourf(x,[0,0]);
        colormap(gray);imagesc(-x); axis equal;  axis tight; axis off; pause(1e-6);
        exportgraphics(gcf, gifName, "Append", true);
          
        % 3D density plot
        figure(2); surf(x); caxis([-12,12]); 
        axis equal; axis([0,nelx,0,nely,-12,12]); view(3);
        
        % Objective function and volume history
        figure(3); subplot(2,1,1); plot(c(1:i),'-'); title('Compliance');
                   subplot(2,1,2); plot(vol(1:i),'-'); title('Volume fraction');
end 
toc
% ========================== EXPORT FINAL RESULTS =========================
finalMatrix = x;

txtName = sprintf('%s_final.txt', outputName);
gifName = sprintf('%s.gif', outputName);

% Salva matriz final
writematrix(x, txtName, 'Delimiter', ',');


disp('---------------------------------------------');
disp(['Arquivos salvos com sucesso:']);
disp([' - Matriz inicial: ' txtNameInitial]);
disp([' - Matriz final:   ' txtName]);
disp([' - Figura final:   ' gifName]);
disp('---------------------------------------------');


%   ======================== AUXILIARY FUNCTIONS ==========================
%%% Optimality Criteria (OC)
function [xnew]=OC(nelx,nely,x,volfrac,dc,mask)
    
% Updatss the density of each element to improve structural performance
l1 = 0; l2 = 100000; move = 0.2; % Bisection method
while (l2-l1 > 1e-4)
  lmid = 0.5*(l2+l1);
  
  % Updates densities respecting movement limits
  xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(-dc./lmid)))));
  xnew(mask < 0.1) = 0.001
  
  % Adjust Lagrange multipliers to meet volume constraints
  if sum(sum(xnew)) - volfrac*nelx*nely > 0;
    l1 = lmid;
  else
    l2 = lmid;
  end
end

%%% Mesh independence filter (check)
function [dcn]=check(nelx,nely,rmin,x,dc)

% Smooths sensitivity (avoids irregular patterns)
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

%%% Finite elements analysis
function [U]=FE(nelx,nely,x,penal,newF,problem)
[KE] = lk;

% Initialize global stiffness matrix and vectors
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
F = sparse(2*(nely+1)*(nelx+1),1); 
U = zeros(2*(nely+1)*(nelx+1),1);

% Set up the global stiffness matrix
for elx = 1:nelx
  for ely = 1:nely
    n1 = (nely+1)*(elx-1)+ely; 
    n2 = (nely+1)* elx   +ely;
    edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;
  end
end


%   ============================= PROBLEMS ================================ 
% Define boundary conditions and loading for each case

%%% CANTILEVER
   if problem == 1
        Name='CANTILEVER';
        F(2*(nelx+1)*(nely+1)-nely,1) = newF; % Force at midpoint
        fixeddofs = (1:2*(nely+1)); 
        alldofs = (1:2*(nely+1)*(nelx+1));
        freedofs = setdiff(alldofs,fixeddofs);
        % SOLVING
        U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
        U(fixeddofs,:)= 0;
   end
%   =======================================================================

%%% REVERSED FRENCH HAND
   if problem == 2
        Name='MÃO FRANCESA INVERTIDA';
        F(2*(nelx)*(nely+1)+2,1) = newF; % Force at the upper end
        fixeddofs = (1:2*(nely+1)); 
        alldofs = (1:2*(nely+1)*(nelx+1));
        freedofs = setdiff(alldofs,fixeddofs);
        % SOLVING
        U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
        U(fixeddofs,:)= 0;
   end
%   =======================================================================

%%% FRENCH HAND
   if problem == 3
        Name='MÃO FRANCESA';
        F(2*(nelx+1)*(nely+1),1) = newF; % Force at thelower end
        fixeddofs = (1:2*(nely+1)); 
        alldofs = (1:2*(nely+1)*(nelx+1));
        freedofs = setdiff(alldofs,fixeddofs);
        % SOLVING
        U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
        U(fixeddofs,:)= 0;
   end
%   =======================================================================
   
%%% BEAM MBB
   if problem == 4
        Name='VIGA MBB';
        F(2*(nelx+1)*(nely+1)-nely,1) = newF; % Force at midpoint
        fixeddofs = (1:2*(nely+1)); 
        alldofs = (1:2*(nely+1)*(nelx+1));
        freedofs = setdiff(alldofs,fixeddofs);
        % SOLVING
        U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
        U(fixeddofs,:)= 0;
   end
%   =======================================================================

%%% MICHELL
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
%   =======================================================================

%   ====================== ELEMENT RIGIDITY MATRIX ========================
function [KE]=lk
E = 1.; % Modulus of elasticity
nu = 0.3; % Poisson ratio

% 8-node quadrilateral element matrix coefficients
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
