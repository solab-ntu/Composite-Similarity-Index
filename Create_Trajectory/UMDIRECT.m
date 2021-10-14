function results = UMDIRECT(fileInfo,lb,ub,options,restartFile)
% results=UMDIRECT(fileInfo,lb,ub,options,restartFile)
%
% This is the University of Michigan's implementation of the DIRECT
% algorithm.
% It solves problems of the form:
%
%    min   f(x)
%   s.to:  g(x) <= 0
%          lb <= x <= ub
%
% It is based upon work from the following papers:
%
% 1) D.R. Jones, C.D. Perttunen and B.E. Stuckman (1993).
% "Lipschitzian Optimization Without the Lipschitz Constant",
% Journal of Optimization Theory and Application, 79(1): 157-181.
%
% 2) D.R. Jones (2001). "DIRECT", Entry in "Encyclopedia of Optimization",
% Kluwer Academic Publishers, 1:431-440.
% 
%
%INPUTS:
%-------
%
% fileInfo - structure array containing following file names:
%   fileInfo.fName           filename for objective function
%   fileInfo.gName           filename for constraint functions
%   fileInfo.fParams         row cell array of additional inputs for fName
%   fileInfo.gParams         row cell array of additional inputs for gName
%
% lb - [1 x nVar] vector of lower bounds on the design variables
%
% ub - [1 x nVar] vector of upper bounds on the design variables
%
% options - structure array of DIRECT specific options:
%   options.maxfCount    maximum number of function calls allowed (default = 200*nVar)
%   options.maxIter      maximum number of iterations allowed (default = 50)
%   options.conTol       [1 x nCon] vector of maximum allowable constraint violation (default = 1e-3)
%                        (if scalar, use same tolerance on all constraints)
%   options.display      scalar value for amount of output sent to screen (default = 0)
%   options.saveFile     save the results structure at each iteration to 'results_fName' (default = 0)
%   options.lgBalance    local / global balance parameter (default = 1e-4)
%   options.localSearch  allow for local optimization (0=no, 1=use SQP, 2=use local-DIRECT)
%   options.termType     choose the termination criterion:
%   						1 = fCount limit exceeded
%							2 = iter limit exceeded
%							3 = either fCount OR iter limit exceeded
%							4 = both fCount AND iter limit exceeded
%							5 = no improvement in last 100 evaluations
%							6 = no improvement in last 10 iterations
%
% restartFile - name of file to load in when restarting DIRECT
%              	(either the name of file or the structure array itself
%						are acceptable)
%
%
%OUTPUTS:
%--------
%
% results - structure array of results containing the following fields:
%   results.xBest        [1 x nVar] vector for best design point
%   results.fBest        scalar best feasible function value at xBest
%   results.gBest        [1 x nCon] vector of constraint values at xBest
%   results.fCount       total number of function evaluations
%   results.iter         number of iterations performed
%   results.lb				 [1 x nVar] vector of lower bounds on design variables
%   results.ub				 [1 x nVar] vector of upper bounds on design variables
%   results.fileInfo     structure array containing the fileInfo pased to UMDIRECT
%   results.options      structure array containing the options pased to UMDIRECT
%   results.fChangeRate  scalar rate of change in f (used for calculating constraint weights)
%   results.gChangeRate  [1 x nCon] vector of rate changes for each constraint
%   results.conWeight    [1 x nCon] vector of constraint weights (for auxiliary function)
%
%   results.hist - structure array containing history of search:
%           hist.fCount  column vector w/ function call number when beat previous fBest
%           hist.iter    column vector w/ iteration number when beat previous fBest
%
%   results.rect - structure array of rectangle properties:
%           rect.dist   		[fCount x 1] vector of vertex distances
%           rect.nDivAll    	[1 x nVar] vector w/ number of divisions along each dimensions for all rects.
%           rect.nDiv        	[fCount x nVar] matrix w/ num. of divisions along each dimension for each rect.
%           rect.x            [fCount x nVar] matrix of design points
%           rect.f            [fCount x 1] vector of objective function values
%           rect.g            [fCount x nCon] matrix of constraint values
%           rect.aux          column vector of Auxiliary Function values (constrained version only)
%           rect.auxCon       column vector of constraint violation portion of auxiliary function
%
%
%EXAMPLE:
%--------
%Here's an example of how one might start a run of DIRECT from scratch.
%
% fileInfo.fName='obj';         %name of file to run objective function
% fileInfo.gName='cons';        %name of file to run constraints
% fileInfo.fParams=[1e4 10];    %extra parameters that must be sent to 'obj.m'
%
% lb=[0 5 0];           		  %lower bounds on x
% ub=[10 25 50];        		  %upper bounds on x
%
% options.maxfCount=500;        %limit on # of function calls
% options.conTol=[1e-6 1e-2];   %constraint tolerance on the 2 constraints
% options.display=1;            %show progress every iteration
% options.saveFile='myfile';    %name of file to save to every iteration
%
% myresults=DIRECT(fileInfo,lb,ub,options);     %start DIRECT
%
%
%
%If the user then wanted to continue with the same problem,
%       they could restart DIRECT with the following.
%
% newoptions.maxfCount=1000;  %set a higher limit on the fCount
%
% mynewresults=DIRECT([],[],[],newoptions,'myfile');
%               OR
% mynewresults=DIRECT([],[],[],newoptions,myresults);
%
%
%
% copyright 2002
% Ryan A. Fellini, Michael J. Sasena, John W. Whitehead
% University of Michigan
% last modified 7/12/02
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Step 0 - Initialize DIRECT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<5,
   restartFile=[];
   if nargin<4
      options=[];
      if nargin<3
         error('Not enough inputs for DIRECT!')
      end
   end
end

%fill in any missing fileInfo arguments
if ~isempty(fileInfo)
   %allow for fileInfo to be empty in case of restart
   if ischar(fileInfo), fileInfo.fName=fileInfo; end
   if ~isfield(fileInfo,'gName'), fileInfo.gName=[]; end
   if ~isfield(fileInfo,'fParams'), fileInfo.fParams=[]; end
   if ~isfield(fileInfo,'gParams'), fileInfo.gParams=[]; end
   if ~isfield(fileInfo,'nArgOut'), fileInfo.nArgOut=[]; end
end

%check if lb and ub are properly assigned
if ~(isempty(lb) | isempty(ub))
   %allow for them to be empty in case of restart
   lb=lb(:)'; ub=ub(:)'; %force row vector
   %check for consistency
   nVar=length(lb);
   if length(ub)~=nVar
      %chech vector sizing
      error('lb and ub must be the same size!')
   end
end

%track which options settings to override if restart is used.
%Also tracks which options were user-specified.
newOptions=options;
if isempty(options), newOptionsNames={};
else, newOptionsNames=fieldnames(options);
end

%fill in any missing options arguments
if isempty(options), options.maxfCount=200*nVar; end %need to initialize it
if ~isstruct(options), error('Options must be a structure array'), end
if ~isfield(options,'maxfCount'), options.maxfCount=200*nVar; end
if ~isfield(options,'maxIter'), options.maxIter=50; end
if ~isfield(options,'conTol'), options.conTol=1e-10; end
if ~isfield(options,'fTol'), options.fTol=1e-10; end
if ~isfield(options,'xTol'), options.xTol=0.001; end
if ~isfield(options,'display'), options.display=1; end
if ~isfield(options,'saveFile'), options.saveFile='noSave'; end
if ~isfield(options,'lgBalance'), options.lgBalance=1e-4; end
if ~isfield(options,'localSearch'), options.localSearch=0; end
if ~isfield(options,'termType'), options.termType=3; end
if ~isfield(options,'termParams'), options.termParams=[]; end

%adjust termination type according to what user specified
if any(strcmp(newOptionsNames,'maxfCount')) & ~any(strcmp(newOptionsNames,'maxIter'))
   options.termType=1;
elseif ~any(strcmp(newOptionsNames,'maxfCount')) & any(strcmp(newOptionsNames,'maxfCount'))
   options.termType=2;
end

%fix any problems with file names in fileInfo or saveFile
if ~isempty(fileInfo)
   if strcmp(fileInfo.fName(end-1:end),'.m')
      fileInfo.fName=fileInfo.fName(1:end-2);
   end
   if ~isempty(fileInfo.gName)
      if strcmp(fileInfo.gName(end-1:end),'.m')
         fileInfo.gName=fileInfo.gName(1:end-2);
      end
   end
end
if ~isempty(options.saveFile)
   if strcmp(options.saveFile(end-3:end),'.mat')
      options.saveFile=options.saveFile(1:end-4);
   end
end

%initialize some variables
stopFlag=0;		%stop when stopFlag = 1
feasFlag=0;		%have we found a feasible point yet?
firstLocal=1;	%allow for an initial local search after 50 fn calls

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Internal Nomenclature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nVar                  number of design variables
% nCon                  number of constraint functions
% iter                  number of iterations
% fCount                number of function evaluations
% stopFlag              stop DIRECT and exit main "do while" loop
% feasFlag              has a feasible point been found
% nRect                 number of rectangles
% currentSet    			set of rectangles selected for division
% fStar                 fmin - lgBalance    (used to calculate dist)
% j                     j = mod(sum(nDiv),nVar)  (used to calculate dist)
% k                     (sum(nDiv) - j)/nVar     (used to calculate dist)
% parent                currently selected rectangle for trisection, 
%                                       resulting in two "children"
% conViol               sum of constraint violations for a given rectangle
% oldfCount             temporary to store number of function calls before
%                                       current iteration
% childDist             Euclidean distance of child from parent 
% divDimen              dimension along which to divide a rectangle


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Step 1 - Evaluate center of design space (or load prior results)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(restartFile)
   %pick up from where DIRECT left off
   if ischar(restartFile)
      %reload prior DIRECT run results from file
      load(restartFile)
   elseif isstruct(restartFile)
      %the structure array was loaded directly
      results=restartFile;
   else
      %otherwise, the user screwed up
      error('restartFile must either be a structure array or a valid filename')
   end
   
   %pull out info from results structure array
   fileInfo=results.fileInfo;
   lb=results.lb;
   ub=results.ub;
   options=results.options;
   xBest=results.xBest;
   fBest=results.fBest;
   rect=results.rect;
   history=results.history;
   fCount=results.fCount;
   iter=results.iter+1;
   
   %determine some parameters from results
   nRect=size(rect.x,1);
   nVar=size(rect.x,2);
   nCon=size(rect.g,2);
   gBest=[];
   newBest=0;
   if ~isempty(history.fCount)
      %if history file isn't empty, a feasible pt has been found
      feasFlag=1;
   end
   if isempty(lb), lb=results.lb; end
   if isempty(ub), ub=results.ub; end
   
   %pull out more info from results for constrained problems
   if nCon > 0
      gBest=results.gBest;
      fChangeRate=results.fChangeRate;
      gChangeRate=results.gChangeRate;
   end
   
   %override any options that were reset in defining the problem
   for i=1:length(newOptionsNames)
      options=setfield(options,newOptionsNames{i}, ...
         eval(['newOptions.',newOptionsNames{i}]));
   end
   
else
   %no resultsFile, so start from scratch...
   
   %first point is the design center point
   rect.x = lb + 0.5*(ub-lb);
   nRect = 1;
   
   %check out fName file to see if there are multiple outputs
   if isempty(fileInfo.nArgOut)
      fileInfo.nArgOut = nargout(fileInfo.fName);
      if fileInfo.nArgOut < 1
         error('Error executing function-calling file! Check file and re-run.');
      end
   end
   
   %evaluate first point
   [rect.f, rect.g] = evalFiles(fileInfo,rect.x);
   
   %count number of constraint functions
   nCon=length(rect.g);
   
   %check for consistency in conTol assignment
   if length(options.conTol)~=nCon
      if length(options.conTol)==1
         options.conTol=options.conTol*ones(1,nCon);
      else
         if nCon ~= 0
            error('length of conTol must agree with number of constraints')
         end
      end
   end
   %force row orientation of conTol vector
   options.conTol=options.conTol(:)';
   
   %check if center point is feasible
   if nCon>0
      if sum(rect.g>options.conTol)==0
         feasFlag=1;
      end
   else
      %unconstrained, thus feasible
      feasFlag=1;
   end   
   
   %track info on rectangle division
   rect.dist=0.5*sqrt(nVar);
   rect.nDiv=zeros(1,nVar);
   rect.nDivAll=zeros(1,nVar);
   
   %initialize some other variables
   iter = 1;
   fCount = 1;
   fChangeRate = 0;
   gChangeRate = zeros(1,nCon);
   
   %start keeping track of history of best points
   if feasFlag
      xBest=rect.x;
      fBest=rect.f;
      gBest=rect.g;
      history.fCount=1;
      history.iter=1;
   else
      xBest=[];
      fBest=Inf; %give a value to show on display Step 5
      gBest=Inf; %give a value to show on display Step 5
      history.fCount=[];
      history.iter=[];
   end
   
end

%initialize results structure array (to force display in a certain order)
results.xBest=[];
results.fBest=[];
if nCon>0, results.gBest=[]; end
results.fCount=[];
results.iter=[];
results.history=[];
results.fileInfo=fileInfo;
results.lb=lb;
results.ub=ub;
results.options=options;
results.rect=[];
if nCon>0
   results.fChangeRate=[];
   results.gChangeRate=[];
   results.conWeight=[];
end

%start sending updates to screen if requested to
if options.display>0
   disp('Iter    Total     Evals        fBest       max(g)')
   disp('        Evals   this iter')
   disp('-----------------------------------------------------')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Begin main iteration loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%will terminate if EITHER fCount or iteration limit is reached
while stopFlag==0
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %  Step 2 - Select set of candidate rectangles for division
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if fCount == 1
      currentSet = [1];	% select only existing rect. for first iteration
   else
      if nCon == 0
         %unconstrained problem
         currentSet = rectSelection(rect,fBest-options.lgBalance);
      else
         %constrained problem
         currentSet = constrSelection(rect,feasFlag,fBest-options.lgBalance);      
      end
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %  Step 3 - Split candidates & evaluate centers of new rectangles
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %track fCount at beginning of iteration and reset newBest
   oldfCount = fCount;
   newBest = 0;
   
   for i=1:size(currentSet,1)
      
      %use this to break loop if searching for 1st feasible pt
      oldFeasFlag=feasFlag;
      
      % name "parent" rectangle index as current rectangle for trisection
      parent = currentSet(i);
      % find all longest sides of rectangle
      divDimen = find(rect.nDiv(parent,:)==min(rect.nDiv(parent,:)));
      if length(divDimen)>1
         % tied, so break by total number of divisions in each dim'n
         divDimen=divDimen(find(rect.nDivAll(divDimen)==min(rect.nDivAll(divDimen))));
         if length(divDimen)>1
            % still tied, so just use lowest dimension
            divDimen=min(divDimen);
         end
      end
      
      % split parent rectangle into two children
      delta = 3^(-(rect.nDiv(parent,divDimen)+1)); %fraction of range
      delta = delta*(ub(divDimen)-lb(divDimen)); 	%scale to true range
      
      xChild=ones(2,1)*rect.x(parent,:); % create children
      xChild(:,divDimen)=xChild(:,divDimen)+[1;-1]*delta; % offset children by delta
      
      % update rect.nDiv and rect.nDivAll
      rect.nDiv(parent, divDimen) = rect.nDiv(parent, divDimen) + 1;
      rect.nDivAll(divDimen) = rect.nDivAll(divDimen) + 1;
      
      % calculate new center-vertex distance (rect.dist) for parent rectangle
      j = mod(sum(rect.nDiv(parent,:)),nVar);
      k = (sum(rect.nDiv(parent,:)) - j)/nVar;
      rect.dist(parent) = 3^(-k)/2*(7/9+nVar-j)^0.5;
      
      for j=1:2
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % evaluate new rectangle center points
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         %call function files
         [fChild,gChild]=evalFiles(fileInfo,xChild(j,:));
         
         %record new evaluations, "child" rectangle inherits rect.dist, rect.nDiv from "parent"
         rect.x=[rect.x ; xChild(j,:)];
         rect.f=[rect.f ; fChild];
         rect.g=[rect.g ; gChild(:)'];
         rect.dist=[rect.dist ; rect.dist(parent)];
         rect.nDiv=[rect.nDiv ; rect.nDiv(parent,:)];
         
         % update rate of change variables here, if constrained
         if nCon > 0
            
            if 0
               %originally used this in updates below, but ChildDist == delta above
               %
               % find Euclidean distance of child from parent (midpoint to midpoint)
               childDist=0;
               for k=1:nVar
                  childDist=childDist+(xChild(j,k)-rect.x(parent,k))^2;
               end
               childDist=sqrt(childDist);         
            end
            
            % roc_obj update
            childChangeRate = abs(fChild-rect.f(parent))/delta;
            fChangeRate=fChangeRate+childChangeRate;
            
            % roc update
            for k=1:nCon
               childChangeRate = abs(gChild(k)-rect.g(parent,k))/delta;
               gChangeRate(k) = gChangeRate(k)+childChangeRate;
            end          
         end
         
         %update fCount and check if there's a new best point
         fCount=fCount+1;
         if fChild<fBest
            if nCon==0 | sum(gChild>options.conTol)==0
               %we'll account for multiple global optima at the end
               history.fCount=[history.fCount ; fCount];
               history.iter=[history.iter ; iter];
               fBest = fChild;
               gBest = gChild(:)';
               xBest = xChild(j,:);
               newBest = 1;
               feasFlag = 1;
            end
         end
         
      end %exit loop for evaluating children of current parent
      
      if oldFeasFlag==0 & feasFlag==1
         break %stop this iteration if just looking for 1st feasible point
      end
      
   end %exit loop for calling CurrentSet
   
   fStar=fBest-options.lgBalance;
   
   % update auxiliary function values if constrained version
   if nCon > 0
      
      % update constraint weights
      for i=1:nCon
         conWeight(i)=fChangeRate/max(gChangeRate(i),10e-10);
      end
      
      % first, update parent rectangles
      for i=1:(fCount-oldfCount)/2
         
         % name "parent" rectangle index as current rectangle to update
         parent = currentSet(i);
         
         % sum constraint violations
         conViol=0;
         for j=1:nCon
            conViol=conViol+conWeight(j)*max([rect.g(parent,j) 0]);
         end
         rect.auxCon(parent) = conViol/rect.dist(parent);
         rect.aux(parent) = (max(rect.f(parent)-fStar,0)/rect.dist(parent))+rect.auxCon(parent);
      end
      
      % next, update children, where there are 2 times as many children as parents
      for i=1:(fCount-oldfCount)
         
         % name "child" rectangle index as current rectangle to update
         child = oldfCount+i;
         
         % sum constraint violations
         conViol=0;
         for j=1:nCon
            conViol=conViol+conWeight(j)*max([rect.g(child,j) 0]);
         end
         rect.auxCon(child) = conViol/rect.dist(child);
         rect.aux(child) = (max(rect.f(child)-fStar,0)/rect.dist(child))+rect.auxCon(child);
      end
      
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %  Step 4 - Check for termination of DIRECT
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   if options.termType==1
      %stop if fCount limit reached
      if fCount>=options.maxfCount
         stopFlag=1;
      end  
      
   elseif options.termType==2
      %stop if iter limit reached
      if iter>=options.maxIter
         stopFlag=1;
      end
      
   elseif options.termType==3
      %stop if EITHER fCount or iter limit reached
      if (fCount>=options.maxfCount) | (iter>=options.maxIter)
         stopFlag=1;
      end
      
   elseif options.termType==4
      %stop when BOTH fCount and iter limit reached
      if (fCount>=options.maxfCount) & (iter>=options.maxIter)
         stopFlag=1;
      end
      
   elseif options.termType==5
      %stop if little improvement has been shown in last 100 fevals
      if ~isfield(options,'termParams'), options.termParams=[100 0.001]; end
      if length(options.termParams)<2, options.termParams(2)=0.001; end
      %start checking once a feasible point is found
      if ~isempty(history.fCount)
         %start checking once taken 100 fn calls after first feasible pt
         if fCount>=options.termParams(1)+history.fCount(1)
            %best point 100 fn calls ago
            oldf=rect.f(history.fCount(max(find( ...
               history.fCount<=(fCount-options.termParams(1))))));
            %latest best point
            %keyboard
            newf=rect.f(history.fCount(end));
            %make sure we've made 0.1% relative improvement in last 100 fn calls
            if newf>eps
               if (oldf-newf)/newf < options.termParams(2)
                  stopFlag=1;
               end
            else
               if (oldf-newf) < options.termParams(2)
                  stopFlag=1;
               end
            end
         end
      end
      
   elseif options.termType==6
      %stop if no improvement has been shown in last 10 iters
      if (iter-history.iter(end))>10
         stopFlag=1;
      end
      
   end
   
   %never exceed maximum allowable fCount
   if fCount>=options.maxfCount
      stopFlag=1;
   end  
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %  Step 5 - Update results to screen and file
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   
   if options.display > 0
      disp(sprintf('%3.0i    (%4.0i)     %3.0i         %7.4f      %7.4f', ...
         iter,fCount,(fCount-oldfCount),fBest,max(gBest)))
   end
   
   % the below saves at every iteration if save file specified, or will
   % save at end of run to default file name.
   if (~strcmpi(options.saveFile,'noSave') | stopFlag)
      %collect results structure array
      results.xBest=xBest;
      results.fBest=fBest;
      results.gBest=[];
      results.rect=rect;
      results.history=history;
      results.fCount=fCount;
      results.iter=iter;
      if nCon > 0
         results.gBest=gBest;
         results.fChangeRate=fChangeRate;
         results.gChangeRate=gChangeRate;
         results.conWeight=conWeight;
      end
      %save results to .mat file if requested
      if ~strcmpi(options.saveFile,'noSave')
         save(options.saveFile,'results');
      end
   end
   
   % update iteration counter
   iter=iter+1;
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %  Step 6 - Launch a local search if necessary
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   
   if options.localSearch & feasFlag & fCount>50 & (newBest|firstLocal)
      
      %set up SQP options using Matlab's optimset
      sqpOptions=optimset('Display','off','DiffMaxChange',1e-2*min(ub-lb),...
         'DiffMinChange',1e-5*min(ub-lb),'LargeScale','off','TypicalX',xBest(1,:),...
         'MaxFunEvals',nVar*100,'MaxIter',100);
      
      warning off
      if nCon==0
         %launch unconstrained local search from xBest
         if ~isempty(fileInfo.fParams)
            [xsqp,fsqp,sqpFlag,sqpOutput] = ...
               fminsearch(fileInfo.fName,xBest(1,:),sqpOptions,fileInfo.fParams{:});
         else
            [xsqp,fsqp,sqpFlag,sqpOutput] = ...
               fminsearch(fileInfo.fName,xBest(1,:),sqpOptions);
         end
      else
         %launch constrained local search from xBest
         if ~isempty(fileInfo.fParams)
            [xsqp,fsqp,sqpFlag,sqpOutput] = fmincon(fileInfo.fName,xBest(1,:), ...
               [],[],[],[],lb,ub,[],sqpOptions,fileInfo.fParams{:});
         else
            [xsqp,fsqp,sqpFlag,sqpOutput] = fmincon(fileInfo.fName,xBest(1,:), ...
               [],[],[],[],lb,ub,[],sqpOptions);
         end
      end
      warning backtrace
      
      %update fCount and iter
      fCount = fCount + sqpOutput.funcCount;
      iter = iter + 1;
      
      %check if there's a new best point
      if sqpFlag & fsqp<=fBest
         %history.fCount=[history.fCount ; fCount];
         %history.iter=[history.iter; iter-1];
         fBest = fsqp;
         if fsqp==fBest
            xBest=[xsqp; xBest];   
         else
            xBest = xsqp;
         end
         feasFlag = 1;
      end
      
      if options.display > 0
         disp(sprintf('%3.0iL   (%4.0i)     %3.0i         %7.4f      %7.4f', ...
            iter-1,fCount,(sqpOutput.funcCount),fBest,max(gBest)))
      end
      
      %update fmin and other variables
      firstLocal = 0;
      
      %REMAINING WORK FOR LOCAL SEARCH TO BE READY:
      %  how to report results every ten fn calls?
      %  how to "count" fn calls for sake of history array?
      %  how to store pts that were run considering rect already exists?
      %		do we use a redundant structure array to store every function call?
      %		do we store just local search pts and when local search started?
      %  fmincon and fminsearch can still go into infinite loops thanks to Matlab
      %  
      
   end %end local search if branch   
   
end %end Master Loop of UMDIRECT

%use the post-processing code "findlocals" to distinguish all
% local optima within a given fTol and xTol of each other
if 0
    [xLocal] = findlocals(rect.x,rect.f,rect.g,options.xTol, ...
        options.fTol,options.conTol,lb,ub);
else
    xLocal=[];
end
results.xLocal=xLocal;

%send warning if couldn't find a feasible point
if ~feasFlag
   %disp('Warning!  No feasible point was found!')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  QUESTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RECENT CHANGES:
%	- added options.fTol and options.xTol to identify unique global optima
%
% QUESTIONS:
%   - how are linear vs. nonlinear constraints handled differently?
%           maybe use linearity for computing exact rate of change info
%
% TO DO LIST:
%	- add options.termParams
%	- fix local search (Step 6)
%	- enable DIRECT-local search possibilities


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%								SUBROUTINE FOR EVALUATING FUNCTIONS						%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fout,gout]=evalFiles(fileInfo,x)

if fileInfo.nArgOut==1
   %run objective function file
   if ~isempty(fileInfo.fParams)
      %include fParams
      fout=feval(fileInfo.fName,x,fileInfo.fParams{:});
   else
      fout=feval(fileInfo.fName,x);
   end
   %run constraint file
   if ~isempty(fileInfo.gName)
      if ~isempty(fileInfo.gParams)
         %include gParams
         gout=feval(fileInfo.gName,x,fileInfo.gParams{:});
      else
         gout=feval(fileInfo.gName,x);            
      end
   else
      gout=[];
   end
else
   %two outputs from fName file
   if ~isempty(fileInfo.fParams)
      %include fParams
      [fout,gout]=feval(fileInfo.fName,x,fileInfo.fParams{:});
   else
      [fout,gout]=feval(fileInfo.fName,x);
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBROUTINE FOR UNCONSTRAINED RECTANGLE SELECTION                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [currentSet] = rectSelection(rect,fStar)
% last updated 1/26/02, MJS


% determine unique set of distances
distList = unique(rect.dist);

% determine best f value for any given distance
fBest=zeros(size(distList)); %initialize fBest to reduce overhead
for i = 1:length(distList)
   fBest(i) = min(rect.f(find(rect.dist == distList(i))));
end

% determine the lowest f value for each distance
setLowest = []; % initialize set of lowest values
multiplesSet = []; % initialize set to hold multiples
for i = 1:length(distList)
   j = find(rect.dist == distList(i) & rect.f == fBest(i));
   if length(j) > 1
      multiplesSet = [multiplesSet j'];
   end
   setLowest = [setLowest j(1)];
end
setLowest = [0 setLowest];

% compute the convex hull for the remaining rectangles
hull = grahamshull([0 distList'],[fStar fBest']);
currentSet = setLowest(hull);
currentSet = currentSet(2:end);
currentSet = union(currentSet,multiplesSet);
currentSet = currentSet(:);	 %make column vector


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grahamshull function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = grahamshull(x,y)

% function h = grahamshull(x,y)
%
% grahamshull returns all points on the convex hull, even redundant ones.
% The vector x must be in ascending order.
%
% Based on the algorithm GRAHAMSHULL in: 
% "Computational Geometry: An Introduction", Franco P. Preparata and
% Michael Ian Shamos, Springer-Verlag, New York, pp. 108-9, 1985.
%

x = x(:);
y = y(:);
m = length(x);
if length(y) ~= m
   error('Input dimensions must agree for convex hull')
end

% Index vector for points in convex hull
h = [1:m]';
if m < 3
   return;
end

v = 1;
flag = 0;


while (next(v,m) ~= 1) | (flag == 0)
   if next(v,m) == m
      flag = 1;
   end
   a = v;
   b = next(v,m);
   c = next(next(v,m),m);
   tol = 1e-6;
   if det([ x(a) y(a) 1 ; x(b) y(b) 1 ; x(c) y(c) 1 ]) >= -tol
      v = next(v,m);
   else
      j = next(v,m);
      x = [ x(1:j-1) ; x(j+1:m)];
      y = [ y(1:j-1) ; y(j+1:m)];
      h = [ h(1:j-1) ; h(j+1:m)];
      m = m - 1;
      v = pred(v,m);
   end
end

function i = next(v,m);
if v == m
   i = 1;
else
   i = v + 1;
end

function i = pred(v,m);
if v == 1
   i = m;
else
   i = v - 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBROUTINE FOR CONSTRAINED RECTANGLE SELECTION                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [currentSet] = constrSelection(rect, feasFlag, fStar)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  constrSelection.m - a function to select rectangles for the constrained DIRECT algorithm
%
%  Last updated: 2/6/02 by J. Whitehead
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERNALLY USED VARIABLES:
% --------------------------
%
% nVar                  number of design variables
% nCon                  number of constraint functions
% fCount                number of function calls made
% currentSet            set of rectangles selected for division
% feasFlag              has a feasible point been found
% fStar                 fmin - lgBalance    (used to calculate dist)
% j                     j = mod(nTri,nVar)  (used to calculate dist)
% k                     (nTri - j)/nVar     (used to calculate dist)
% stopFlag              stop DIRECT and exit main "do while" loop            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nRect = size(rect.x,1);
tieTol = abs(min([((1e-6)*fStar) 1e-6]));	% relative tolerance to determine intersection ties
currentSet = [];        		% initialize set of selected rectangles
domin = zeros(1, nRect);  		% initialize rectangle dominance check variable
ties = [];                 	% initialize ties vector

% if feasible, select rectangles by lower envelope to create 'currentSet'
if (feasFlag == 1)
   
   % select rectangle(s) with lowest rect.aux value as first rectangle
   % sorting takes too long. Find min value and all ties.
   minVal = min(rect.aux);
   for i=1:nRect
      if (rect.aux(i) <= (minVal + tieTol))
         currentSet = [currentSet ; i];
         ties = [ties i];
      end
   end
   
   breakFlag = 0;
   maxIntersect = fStar;
   
   while breakFlag == 0
      
      % 1. if necessary, do tie-breaking to select which rectangle to follow, put that index in currRect
      if (length(ties) > 1)
         % select first rectangle in 'ties' vector
         incumbent = ties(1);
         ties(1) = [];
         
         for i=1:length(ties)
            % check if right side intersect of incumbent
            if ( maxIntersect > rect.f(incumbent) )
               
               % check if also right side intersect of other rect, otherwise follow incumbent
               if ( maxIntersect > rect.f(ties(i)) )                    
                  % check if other rectangle should be followed, instead of incumbent
                  
                  % if other is further left, follow it
                  if ( rect.f(ties(i)) < rect.f(incumbent) )
                     incumbent = ties(i);
                     
                     % special case, if both intersect on left and have same f-value
                  elseif ( rect.f(ties(i)) == rect.f(incumbent) )
                     if ( rect.dist(ties(i)) > rect.dist(incumbent) )
                        incumbent = ties(i);
                     end
                  end 
               end
               
            else % is left side intersect of incumbent
               
               % if on right side of other rect, follow it
               if ( maxIntersect > rect.f(ties(i)) )
                  incumbent = ties(i);
                  
               else % is on left side of other rect, so check slopes
                  
                  if ( rect.dist(ties(i)) > rect.dist(incumbent) )
                     incumbent = ties(i);
                  end
               end
               
            end
         end
      else
         incumbent = ties(1);	% if none tied, pick only element in 'ties' vector
      end
      %keyboard
      ties = [];
      currRect = incumbent;
      
      % 2. once selected one rectangle to follow, do dominance check to weed out others
      for i = 1:nRect
         
         % if already dominated, ignore
         if ( (domin(i) ~= 1) & (i ~= currRect) )
            
            if (rect.dist(currRect) >= rect.dist(i))
               if ( (rect.auxCon(currRect) <= rect.auxCon(i)) & (rect.f(currRect) < rect.f(i)) )
                  domin(i) = 1;	%rectangle i is dominated
                  
               elseif ( (max([(rect.f(currRect) - rect.f(i)) 0]) / rect.dist(currRect)) ...
                     < (rect.auxCon(i) - rect.auxCon(currRect)) )
                  domin(i) = 1; % rectangle 'i' is dominated, keep track of it
               end
            end
         end
      end
      
      % 3. now check where all other rectangles intersect with sloped (left) side of current rectangle
      breakFlag = 1; 
      intersect = [];
      for i = 1:nRect
         % if dominated, ignore 
         if (domin(i) ~= 1) & (i ~= currRect)
            % calculate the f-star value for right and left side intersection
            rSideInter = rect.f(currRect) +(rect.auxCon(currRect) - rect.auxCon(i))*rect.dist(currRect);
            
            % for left side intersect, make sure not to get divide by zero error!
            if (rect.dist(currRect) ~= rect.dist(i))
               lSideInter = (rect.dist(currRect)*rect.dist(i)*(rect.auxCon(i)-rect.auxCon(currRect))+rect.f(i)*rect.dist(currRect) - ...
                  rect.f(currRect)*rect.dist(i))/(rect.dist(currRect) - rect.dist(i));
            else
               lSideInter = -inf;
            end
            
            % check for false intersections
            if (rSideInter < rect.f(i)) | (rSideInter > rect.f(currRect)) | (rSideInter >= (maxIntersect-tieTol))
               rSideInter = -inf;
            end
            if (lSideInter > rect.f(currRect)) | (lSideInter > rect.f(i)) | (lSideInter >= (maxIntersect-tieTol))
               lSideInter = -inf;
            end
            
            % select closest intersection point for that rectangle, i.e. the maximum one
            intersect(i) = max(rSideInter, lSideInter);
            
            % if intersection is a 'true' intersection, make sure break flag is false
            if (intersect(i) > -inf)
               breakFlag = 0;
            end
         else
            intersect(i) = -inf;	% if dominated or current rectangle, set intersection to -inf
         end
      end
      if breakFlag == 0
         
         % 4. select intersection point "closest" to current point, find all ties, go to step 1.
         maxIntersect = max(intersect);
         for i = 1:nRect
            % ignore current rectangle
            if (i ~= currRect)
               
               % select as part of lower envelope all rectangles within the tie tolerance
               if ( intersect(i) >= (maxIntersect - tieTol) )
                  currentSet = [currentSet; i];
                  
                  % select only those rectangles which have exact same intersect point for following along lower envelope
                  if ( intersect(i) == maxIntersect )
                     ties = [ties i];
                  end   
               end
            end
         end 
      end
      
   end % 'while breakFlag'
   
   
   % if not feasible, select the rectangle that minimizes the rate of
   % change necessary to bring the constraint violations down to zero
else
   % sorting takes too long. Find min value and all ties.
   minRect = min(rect.auxCon);
   for i=1:nRect
      if (rect.auxCon(i) <= (minRect + tieTol))
         currentSet = [currentSet; i];
      end
   end
end

currentSet = unique(currentSet);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot envelope of rectangles (for debugging)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ploton=0;
s = {'m-'; 'c-'; 'r-'; 'g-'; 'b-'; 'm:' ;'c:'; 'r:'; 'g:' ;'b:'};
if ploton==1
   clf
   hold on
   for i=1:nRect
      x=[(rect.f(i)-10) rect.f(i) (rect.f(i)+10)];
      y(1) = (10/rect.dist(i))+rect.auxCon(i);
      y(2) = rect.auxCon(i);
      y(3) = rect.auxCon(i);
      plot(x,y,s{(1+mod(i,9))})
      keyboard
   end
   hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SUBROUTINE FOR LOCATING ALL LOCAL OPTIMA FROM A GIVEN DATA SET AND TOL %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xLocal,sortID]=findlocals(x,f,g,xTol,fTol,gTol,lb,ub)

if nargin<8
   ub=[];
   if nargin<7
      lb=[];
      if nargin<6
         gTol=[];
         if nargin<5
            fTol=[];
            if nargin<4
               xTol=[];
               if nargin<3
                  g=[];
                  if nargin<2
                     error('Error!  Findlocals requires at least two arguments!')
                  end
               end
            end
         end
      end
   end
end

if isempty(xTol), xTol=0.01; end
if isempty(fTol), fTol=0.01; end
if isempty(gTol), gTol=1e-6; end
if isempty(lb), lb=min(x); end
if isempty(ub), ub=max(x); end

%find all feasible points
if isempty(g)
   feasID=[1:length(f)];
else
   feasID = find(all(g<=gTol,2));
end

%find best feasible point
fBest = min(f(feasID));

%find all feasible points within fTol of fBest
candidates = feasID(find(f(feasID)<=fBest+fTol*abs(fBest)));
[fSort,sortID] = sort(f(candidates));

%make ID number of original sample, not of candidates
sortID = candidates(sortID);

%find which of these candidates are xTol away from the nearest best candidate
% by sweeping through their objective function values and checking distances
xLocal=x(sortID(1),:);

for i=2:length(candidates)
   %update current point to check
   xCurrent=x(sortID(i),:);
   xBetter=x(sortID(1:i-1),:);
   
   %find distances to nearest known better solution
   clear dist
   for ii=1:size(xBetter,1)
      dist(ii)=norm((xCurrent-xBetter(ii,:))./[ub-lb]);
   end
   
   %if it's the best point within xTol, call it a local optimum
   % also accounts for situations where two points have identical f
   if min(dist)>xTol | f(sortID(i))==f(sortID(i-1))
      xLocal = [xLocal ; xCurrent];
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
