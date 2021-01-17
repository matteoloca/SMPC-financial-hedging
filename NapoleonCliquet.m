function f = NapoleonCliquiet(numTimes,numPaths,fixingDates,previousPrices,startingT)
%BARRIEROPTION Price European barrier call option
%
% Syntax:
%
%   f = BarrierOption(numTimes,numPaths)
%
% Description:
%
%   End-of-period processing function to price European up-and-out (barrier)
%   call options by Monte Carlo simulation, assuming a constant risk-free
%   rate. The function also illustrates how to save the maximum and final
%   prices on each sample path, and how to update, access, and share
%   information among several nested functions.
%
% Input Arguments:
%
%   numTimes - Number of simulation times steps (numPeriods*numSteps).
%   numPaths - Number of independent sample paths (simulation trials).
%
% Output Argument:
%
%   f - Structure of nested function handles to save the maximum and final 
%       prices on a path-by-path basis, and to compute prices and standard 
%       errors of European up-and-in (barrier) call options.

% Copyright 1999-2010 The MathWorks, Inc.
% $Revision: 1.1.4.1 $   $Date: 2013/10/09 06:17:34 $

prices = zeros(numPaths,numTimes);  % Pre-allocate terminal price vector
iTime = startingT;                   % Counter for the period index
iPath = 1;                   % Counter for the trial index
T = 0;                       % Initialize time-to-expiry
tStart = 0;
simTimes=numTimes-startingT;

fixingDates=[fixingDates];
times = size(fixingDates);
times=times(2);

prices(:,1:startingT)=repmat(previousPrices',numPaths,1);

variance=zeros(numPaths,numTimes);

f.Sim = @sim;     % End-of-period processing function
f.OptionPrice = @getOptionPrice; % Barrier option pricing utility
f.FinalOptionPrices=@getFinalOptionPrices;
f.StandardError = @getStandardError;    % Standard error calculator utility
f.AllOptionPrices=@getAllOptionPrices;
f.Prices = @getPrices;

f.Variance= @getVariance;

function X = sim(t,X)
  if t==0
      %prices(:,startingT) = X(1);
      X(1);
      iTime = startingT + 1;
  else
     if iTime < numTimes
        prices(iPath,iTime) = X(1);
        if size(X,1)>1
        variance(iPath,iTime) = X(2);
        end
        iTime = iTime + 1;
        %t% Update the period counter
     else % The last period of the current path
        prices(iPath,numTimes) = X(1); % Save the terminal price for this path
        if size(X,1)>1
        variance(iPath,numTimes) = X(2);
        end
        iTime = startingT+1;         % Re-set the time period counter
        iPath = iPath + 1; % A new path will begin next time
        T = t - tStart;    % Accumulate the time-to-expiration
     end
  end
  %end
  
end

function value = getOptionPrice(strike,rate,C)
  po=zeros(numPaths,times-1);
  
  po=((prices(:,fixingDates(2:times))-prices(:,fixingDates(1:times-1)))./(prices(:,fixingDates(1:times-1))));
  values = zeros(numPaths,1);
  
  values = exp(-rate*T)*max(C+min(po'),0);
    
  %values = max(C+((prices(:,2)-prices(:,))/prices(:,1:Tfix-1)),0);
  value = mean(values);
  
end



function value = getAllOptionPrices(strike,rate,dt)
  time=T:-dt:0;
  indexes=numTimes+1:-1:1;
  value=zeros(numTimes+1);
      value = exp(-rate*time).*max(prices(:,numTimes)-strike,0);
      value=mean(value');
  
end

function value = getFinalOptionPrices(strike,rate)
value = exp(-rate*T)*max(prices(:,numTimes)-strike,0);
end

function value = getPrices()
  value=prices(1,startingT:numTimes);
end

function value =getVariance()
    value=variance(1,startingT:numTimes);
end

function value = getStandardError(strike,rate,Antithetic)
  values = exp(-rate*T)*max(prices-strike,0);
  if (nargin >= 4) && Antithetic
     value = std((values(1:2:end-1)+values(2:2:end))/2)/sqrt(numPaths/2);
  else
     value = std(values)/sqrt(numPaths);
  end
  
end

end % End of outer/primary function
