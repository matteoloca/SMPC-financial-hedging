function f = BarrierOption(numTimes,numPaths)
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
maximum = zeros(numPaths,1); % Pre-allocate maximum  price vector
iTime = 1;                   % Counter for the period index
iPath = 1;                   % Counter for the trial index
T = 0;                       % Initialize time-to-expiry
tStart = 0;                  % Initialize sample time
variance=zeros(numPaths,numTimes);

f.Sim = @sim;     % End-of-period processing function
f.OptionPrice = @getBarrierOptionPrice;% Barrier option pricing utility
f.FinalOptionPrices=@getFinalBarrierOptionPrice;
f.StandardError = @getStandardError;    % Standard error calculator utility
f.AllOptionPrices=@getAllOptionPrices;
f.Prices = @getPrices;
f.Variance= @getVariance;

function X = sim(t,X)
    
   maximum(iPath) = max(maximum(iPath),X(1));
   if t==0
      prices(:,1) = X(1);
      if size(X,1)>1
        variance(:,1) = X(2);
      end
      iTime = iTime + 1;
   else
     if iTime < numTimes
        prices(iPath,iTime) = X(1);
        if size(X,1)>1
        variance(iPath,iTime) = X(2);
        end
        iTime = iTime + 1; % Update the period counter
     else % The last period of the current path
        prices(iPath,numTimes) = X(1); % Save the terminal price for this path
        if size(X,1)>1
        variance(iPath,numTimes) = X(2);
        end
        iTime = 2;         % Re-set the time period counter
        iPath = iPath + 1; % A new path will begin next time
        T = t - tStart;    % Accumulate the time-to-expiration
     end
  end
  
end

function value = getBarrierOptionPrice(strike,rate,barrier)
    
  values = zeros(numPaths,1);
  i = maximum < barrier; % Paths on which barrier was exceeded
  values(i) = exp(-rate*T)*max(real(prices(i,numTimes))-strike,0);
  value = mean(values);
  
end

function value = getFinalBarrierOptionPrice(strike,rate,barrier)
    
  values = zeros(numPaths,1);
  i = maximum < barrier; % Paths on which barrier was exceeded
  values(i) = exp(-rate*T)*max(prices(i,numTimes)-strike,0);
  %value = mean(values);
  
end

function value = getAllOptionPrices(strike,rate,dt)
  time=T:-dt:0;
  i = maximum < barrier;
  values = exp(-rate*time).*max(prices(i,numTimes)-strike,0);
  value = mean(values);
  
end

function value = getPrices()
  value=prices(1,:);
end

function value = getStandardError(strike,rate,barrier,Antithetic)
    
  values = zeros(numPaths,1);
  i = maximum < barrier; % Paths on which barrier was exceeded
  values(i) = exp(-rate*T)*max(prices(i)-strike,0);
  if (nargin >= 4) && Antithetic
     value = std((values(1:2:end-1)+values(2:2:end))/2)/sqrt(numPaths/2);
  else
     value = std(values)/sqrt(numPaths);
  end
  
end

function value =getVariance()
    value=variance;
end
end % End of outer/primary function
