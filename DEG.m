function [Out_DispEn, npdf]=DEG(x,Adj,m,L,nc,MA)
% This function calculates dispersion entropy for graphs (DEG), using
% different mapping approaches (MA).
% Works for Adjacency matrix (directed and undirected graphs) and also
% weighted graphs
%
% Inputs:
%
% x: graph signal - a vector of size 1 x N (the number of vertices)
% Adj: square matrix N x N  The elements of the matrix indicate whether pairs of vertices are adjacent or not in the graph.
% m: embedding dimension
% nc: number of classes (it is usually equal to a number between 2 and 9. Worth noting that the maximum number of classes for the code is 9.
% MA: mapping approach, chosen from the following options:
% 'LM' (linear mapping),
% 'NCDF' (normal cumulative distribution function),
% 'TANSIG' (tangent sigmoid),


% Output:
% Out_mvMDE: a scalar value between 0 and 1 giving the Dispersion Entropy value of the graph signal.
% npdf: Normalised probability density function of graph dispersion
% patterns.
%
% Ref:
%
% Emails: john.fabila@ed.ac.uk / javier.escudero@ed.ac.uk
%  21-April-2021
%

%%Default to some input parameters for clarity and to control errors.
if nargin > 6
    error('too many inputs');
end
if nargin < 6
    MA = 'NCDF';
end
if nargin < 5
    nc = 5;
end
if nargin < 4
    L = 1;
end
if nargin < 3
    m = 2;
end

% Initial checks
N = length(x); % N=number of vertices

if or(size(Adj,1)~=N ,size(Adj,1)~=N )
        error('not square matrix');
else

% Step 1: Matrix Lap where first column is the signal, and the columns the
% j-average
    Lap = zeros(N,m);%Set_Matrix_store
    Lap(:,1)=x';%First Column is the signal
    Aux=eye(N);
    t=2;
    %Each column of Lap is an average of the j-neighbourhoods
    for j=L:L:(m-1)*L
        Aux=Aux*(Adj^(L));
        Lap(:,t)=sparse(diag(1./sum(Aux,2)))*(Aux*x');
        t=t+1;
    end

    %%Delete row, only it is used for directed graphs, otherwise, this step
    %%does not do nothing
    Lap(any(isnan(Lap),2),:) = [];
    N=size(Lap,1);

  % Step 2:  Normalization of the matrix X=Lap according several mappings
  % Mapping approaches
    X=Lap;
            sigma=std(X(:,1));
            mu=mean(X(:,1));
%X
%Lap
        for ii=1:m  
            h=X(:,ii);
            switch MA
            case   'LM'
                y=mapminmax(h',0,1);
                y(y==1)=1-1e-10;
                y(y==0)=1e-10;
                X(:,ii)=round(y*nc+0.5)';
                
            case 'NCDF'
                y=normcdf(h',mu,sigma);
                y=mapminmax(y,0,1);
                y(y==1)=1-1e-10;
                y(y==0)=1e-10;
                X(:,ii)=round(y*nc+0.5)';

            case 'LOGSIG'
                y=logsig((h'-mu)/sigma);
                y=mapminmax(y,0,1);
                y(y==1)=1-1e-10;
                y(y==0)=1e-10;
                X(:,ii)=round(y*nc+0.5)';

            end     
        end
%X
%N
% Step 3: Generation of all possible DisEn patterns (using a recursive process:)

all_patterns=(1:nc)';

for f=2:m
    temp=all_patterns;
    all_patterns=[];
    j=1;
    for w=1:nc
        [a,~]=size(temp);
        all_patterns(j:j+a-1,:)=[temp,w*ones(a,1)];
        j=j+a;
    end
end

pattern_num=nc^m;% obtain the number of all possible dispersion pattern
key = sum(all_patterns.*(ones(pattern_num,1)*(10.^flip(0:(m-1)))),2);

% Step 4: Array with the dispersion pattern of the graph
rescom = sum(X.*(ones(N,1)*(10.^flip(0:(m-1)))),2);


% Step 5: Calculation of dispersion pattern frequency in the graph
pdf=zeros(1,pattern_num);
for id=1:pattern_num
    pdf(id)=sum(rescom==key(id));
end

npdf=pdf/N;
p=npdf(npdf~=0);
%Out_DispEn = -sum(p .* log(p))
Out_DispEn = -sum(p .* log(p))/log(pattern_num);
end


