function Bcon=conload(X,model,mode,options);
%CONLOAD Congruence loadings for PARAFAC, TUCKER and NPLS
%  Determines congruence (earlier known as correlation) loadings 
%  for a specific mode of a model. Congruence loadings look at 
%  "non-average correlations", hence take differences in offset into 
%  account. 
%
%  Note that due to non-orthogonal loadings in PARAFAC, 
%  individual correlations can add to more than 1. Therefore, 
%  such loadings are not drawn with ellipses but squares 
%  added. Use options.force = 'ellipse' or 'square' to force 
%  one or the other on the plot. 
%
%  INPUTS:
%        X = nway data
%    model = parafac, tucker or npls model (plstoolbox structure)
%     mode = loading mode to investigate (i.e. mode = 1 for 
%                                             samples if they are in the 
%                                             first mode)
%
%  OPTIONAL INPUTS
%  options = modify options
%
%  OUTPUT:
%     Bcon = Congruence loadings
%
% I/O: Bcon = conload(X,model,mode); 

%Copyright Rasmus Bro, 2004
%Licensee shall not re-compile, translate or convert "M-files"
% for use with any software other than MATLAB®, without
% written permission

% Apr, 2006, RB, Ver 1.02, included PARAFAC2
% Feb, 2005, RB, Ver 1.01, fixed bug in normalization


if nargin == 0; X = 'io'; end
varargin{1} = X;
if ischar(varargin{1});
  options = [];

  options.plots = 'on';
  options.force = 'off';
  if nargout==0; 
    evriio(mfilename,varargin{1},options); 
  else; 
    Bcon = evriio(mfilename,varargin{1},options); 
  end
  return;
end

%if nargin < 4 | isempty(options);
%  options = conload('options');
%else
%  options = reconopts(options,'corrload');
%end


if isa(X,'dataset')
  inc = X.includ;
  X = X.data(inc{:});
end

if isa(X,'dataset')% Then it's a SDO
  if iscell(X.data)
    X = cell2array(X);
  end
else
  if iscell(X)
    X = cell2array(X);
  end
end

%if any(isnan(X(:)))
%    xhat = datahat(model);
%    X.data(isnan(X(:)))=xhat(isnan(X(:)));
%end


if isa(model,'cell') % if using Factors (cell format) from the nway toolbox
  model_loads = model;
end

xsize = size(X);
Bcon=repmat(NaN,size(model_loads{mode}));
Bcorr=repmat(NaN,size(model_loads{mode}));

% MAKE A CORE ARRAY
dotuck = 0;
if length(model_loads)==length(xsize)+1 % Then it's a tucker model with a core in the last cell element
  core = model_loads{end};
  dotuck = 1;
elseif isfield(model,'core') % Then NPLS
  core = model.core{end};
  dotuck = 1;
else % Then PARAFAC
  % Ideal superdiagonal array of ones
  Fac = size(model_loads{1},2);
  core = zeros(Fac*ones(1,length(xsize)));
  j=length(core(:));
  core(linspace(1,j,Fac))=1;
end

% rearrange so that the mode of interest is row mode
X = permute(X,[mode 1:mode-1 mode+1:length(xsize)]);
loads = model_loads([mode 1:mode-1 mode+1:length(xsize)]);
core = permute(core,[mode 1:mode-1 mode+1:length(xsize)]);


if isstruct(model_loads{1}) % Then it's PARAFAC2 - treat it separately
  if mode==1 % First mode is special
    Bcon = cell(1,size(model_loads{end},1));
    for k=1:size(model_loads{end},1) % Number of samples (last mode)
      for f = 1:size(model_loads{2},2) % For every factor
        m = kron(loads{end}(k,f),loads{end-1}(:,f));
        for o = length(loads)-2:-1:2
          m = kron(m,loads{o}(:,f));
        end
        M(:,f)=m;
        for i=1:size(X,1);
          xik = reshape(X(i,:),xsize(2:end));
          xik = permute(xik,[ndims(xik) 1:ndims(xik)-1]);
          xik = xik(k,:);
          xik = xik(:);
          b(i,f)= xik'*M(:,f)/sqrt( (xik'*xik)*(M(:,f)'*M(:,f)) );
        end
      end
      Bcon{k}=b;
    end
  else % All other modes than mode 1
    if mode ==length(xsize) % Last mode is also special!
      for k=1:size(model_loads{end},1) % Number of samples (last mode)
        xk = X(k,:)';
        PH = model_loads{1}.P{k}*model_loads{1}.H;
        for f = 1:size(model_loads{2},2) % For every factor
          m = model_loads{end-1}(:,f);
          for o = length(model_loads)-2:-1:2
            m = kron(m,model_loads{o}(:,f));
          end
          m = kron(m,PH(:,f));
          M(:,f)=m;
          Bcon(k,f)= xk'*M(:,f)/sqrt( (xk'*xk)*(M(:,f)'*M(:,f)) );
        end
      end
      Overall = nm(M);
    else % All 'middle' modes
      for f = 1:size(model_loads{2},2) % For every factor
        m1 = 1;
        for o = length(loads)-1:-1:3
          m1 = kron(m1,model_loads{o}(:,f));
        end
        m=[];
        for k=1:size(model_loads{end},1) % Number of samples (last mode)
          PH = model_loads{1}.P{k}*model_loads{1}.H;
          m = [m;kron(m1,PH(:,f)*loads{end}(k,f))];
        end
        M(:,f)=m;
        for i = 1:size(X,1)
          Bcon(i,f)= X(i,:)*M(:,f)/sqrt( (X(i,:)*X(i,:)')*(M(:,f)'*M(:,f)) );
        end
      end
      Overall = nm(M);
    end
  end

else

    % ALL OTHER MODELS THAN PARAFAC2

    Overall = [];
    for fac = 1:size(Bcon,2)
      % Make loads for component fac
    % M = outerm(L);
    if ndims(X) == 2,
      M = loads{2}(:,fac); % outerm(L) does not work for 2-way data (ndims(L) = 1;)
    else
      if dotuck
        z = kron(loads{end},loads{end-1});
        for j = length(xsize)-2:-1:2
          z = kron(z,loads{j});
        end
        M = core(fac,:)*z';
      else
        %%%% the mode L{1} is not to be included in the submodel (following the 2-mode correlation formula...)
        for m = 2:length(loads)
          L{m-1} = loads{m}(:,fac);
        end
        M = L{1} * L{2}'; %outerm(L);
      end
    end
    for i=1:xsize(mode)
      temp = (X(i,:)); Xi = temp(:);
      r = misscorrmatrix([M(:) Xi],0);
      % r = misscorrmatrix([vec(M(i,:)) vec(X(i,:))],0);
      Bcorr(i,fac) = r(1,2);
      r = misscorrmatrix([M(:) Xi],2);
      %     a = vec(M(i,:))'*vec(M(i,:));
      %     b = vec(X(i,:))'*vec(X(i,:));
      %     r = vec(M(i,:))'*vec(X(i,:))/(sqrt(a)*sqrt(b))
      Bcon(i,fac) = r(1,2);
    end
    Overall = [Overall M(:)];
  end
  Overall = nm(Overall);
end

% if strcmp(lower(options.plots),'on')
%   if iscell(Bcon) % PARAFAC2 mode 1
%     disp([' Plotting of PARAFAC2 mode 1 loadings not enabled. The set of ',num2str(length(Bcon)),' loading matrices are held in a cell array that you may plot yourself'])
%   else
%   if size(Bcon,2)>1
%     plot(Bcon(:,1),Bcon(:,2),'o')
%     text(Bcon(:,1),Bcon(:,2),num2str([1:size(Bcon,1)]'))
%     xlabel('Component 1')
%     ylabel('Component 2')
%     
%     
%     % Chk if loadings orthogonal, if so plot ellipse, else squares
%     if strcmp(lower(options.force),'ellipse')|strcmp(lower(options.force),'square')
%       % If known and given in options, save all the computations 
%       % sto figure out if loadings are orthogonal
%     else
%       LtL = Overall'*Overall;
%       I = size(LtL,1);
%       pI = ones(I,I)-eye(I);
%       LtLI = LtL.*pI;
%       m = max(abs(LtLI(:))); % Max of off-diagonal elements
%       if m<.001 % Then zero => ellipses
%         options.force = 'ellipse';
%       else
%         options.force = 'square';
%       end
%       
% %       % OLD version incorrectly checked the mode under consideration instead of the other modes      
% %       L = loads{1}*diag(sum(loads{1}.^2).^(-.5)); % Normalize
% %       LtL = L'*L;
% %       I = size(LtL,1);
% %       pI = ones(I,I)-eye(I);
% %       LtLI = LtL.*pI;
% %       m = max(abs(LtLI(:))); % Max of off-diagonal elements
% %       if m<eps*100 % Then zero => ellipses
% %         options.force = 'ellipse';
% %       else
% %         options.force = 'square';
% %       end
%     end
% 
%     if strcmp(lower(options.force),'ellipse')
%       h=ellps([0 0],[1 1],'-k');
%       set(h,'linewidth',1.3);
%       grid on
%     else
%       hold on
%       plot([-1 1],[1 1],'k')
%       plot([-1 1],[-1 -1],'k')
%       plot([1 1],[-1 1],'k')
%       plot([-1 -1],[-1 1],'k')
%       hold off
%       axis([-1.1 1.1 -1.1 1.1])
%       grid on
%     end
%   else
%     bar(Bcon(:,1))
%     xlabel('Variable')
%     ylabel('Correlation loading')
%   end
%   end
% end



function [B]=misscorrmatrix(A,type);
%B=misscorrmatrix(A,type);
%
% Claus A. Andersson, claus@andersson.dk, Copyrighted (C) 98
% Rasmus Bro, rb@kvl.dk, modified to enable congruence, 2004
%
% This file derives a correlation matrix of a matrix that has
% missing values. Note that all columns must have at least two
% non-missing observations in common with all other columns.
% Otherwise the returned correlation mapping has nans
% on the corresponding positions.
%
% A      :  Matrix (Samples x Variables) with missing values.
% type   :  Defines the type of operation
%           0: Correlation (std&mean), [-1:1] (default)
%           1: Covariances (mean), [-inf:inf]
%           2: Congruence (correlation but without centering)
%
% Note: If any variables are constant over the samples, the results
%       will have nans in them.
%
% Two steps to make a mapping:
%   B=misscorrmatrix(A); %Step 1
%   imagesc(B);colorbar;grid on; %Step 2


if ~exist('type'),
  type=0;
end;

i_A=size(A,1);
j_A=size(A,2);
isnanA=isnan(A);
sumisnanA=sum(isnanA);
B=nan*zeros(j_A);

nan_list=find(sumisnanA);
true_list=find(~sumisnanA);

for i=1:j_A,
  true_ivec=~isnanA(:,i);
  for j=i+1:j_A,
    true_jvec=~isnanA(:,j);
    list=true_ivec.*true_jvec;
    divis=1/sum(list);
    idx=find(list);

    ivec=A(idx,i);
    jvec=A(idx,j);

    if type~=2 % Congruence, hence don't subtract averages
      ivec=ivec-sum(ivec)*divis;
      jvec=jvec-sum(jvec)*divis;
    end

    if type==0|type==2,
      ivecsum=sum(ivec.^2);
      if ivecsum>eps,
        ivec=ivec/sqrt(ivecsum);
      else
        ivec=ivec*0;
      end;

      jvecsum=sum(jvec.^2);
      if jvecsum>eps,
        jvec=jvec/sqrt(jvecsum);
      else
        jvec=jvec*0;
      end;
    end;

    B(i,j)=ivec'*jvec;
  end;
end;


function Xn = nm(X,dontsignswitch);

if nargin<2
  dontsignswitch=0;
end
  
Xn=X;
for i=1:size(X,2)
   Xn(:,i)=Xn(:,i)/norm(Xn(:,i));
end

if ~dontsignswitch
Xn = Xn * diag(sign(sum(Xn.^3)));
end