function [flux,classes,DSE,tailles,varargout]=nb_flux(X,mini,maxi,pas,A,B,varargin)
% Conversion des matrices de nombre de particules par litre et par classe 
% de taille en flux
%
%           [flux,classes,milieu_classes,tailles]=nb_flux(X,mini,maxi,pas,A,B)
% Input:
%   X : Matrice des concentrations particulaires en Nb/l
%   mini : Borne inferieure de la premiere classe (biovolume)
%   maxi : Borne superieure de la derniere classe (biovolume)
%   pas : Pas pour passer d'une classe a l'autre (2^1/3 pour ECOpart)
%   A : Coefficient de la relation exponnentielle (defaut 109.5)
%   B : Exposent de la relation exponnentielle    (defaut 3.52)
% 
% Output:
%   flux : Matrice correspondante des flux en mg.m^(-2).j^(-1)
%   classes : Bornes des classes de taille utilisees en mm
%   milieu_classes : Valeur du milieu des classes en mm
%   tailles : Tailles de la classe en mm


[a,b]=size(X);
flux=zeros(a,b);

[nl,nc]=size(X);
if nargin==1            % Pour les 22 classes de taille
    mini=7.5000e-05;
    maxi=7000;
    pas=2;
end

if nargin<5
    A=109.5;
    B=3.52;
end

mini=esdEdges(1); %1*10^-3;       % in mm

maxi=esdEdges(end); %30;            % in mm

pas=2^(1/3);

% Calcul des classes de taille (classes de volume !!!!!)
[classe,taille,medi]=pasvar(mini,maxi,pas);
% classe = esdEdges;
% taille = esdBinWidth;
% medi = esdMiddle;

%%%%%%%%%%%%%%% Troisieme partie (calcul du poids sec et de la sedimentation) %%%%%%%%%%%%
% Transformation des medi0 en volume, en DSE (necessaire pour le poids sec et les flux)
ra=medi/((4/3)*pi);
ra=ra.^(1/3);
DSE=ra*2;

% Calcul du flux
flu=A*DSE.^B;        % Formule pour le calcul du flu basee sur les minimisation sur les pieges a sediment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fin de la Troisieme partie%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:nl
    flux(i,:)=X(i,:).*flu;                     % Matrice des flux (en mg.m^(-2).j^(-1))
end

if nargout>=4
    ra=classe/((4/3)*pi);
    ra=ra.^(1/3);
    classes=ra*2;
    
    tailles=classes(2:end)-classes(1:end-1);
end
