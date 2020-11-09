classdef SlugDerivativeFilter < Filter
    %SLUGTESTFILTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    ProcessParameters
    Slugs
    Nozzle

    end
    
    methods
        function this = SlugDerivativeFilter()
% Set slugtest filter properties

 prompt = {'Enter nozzle diameter (mm)'};
            title = 'SlugDerivativeFilter properties';
            dims = [1 35];
            definput = {'20'} ;
            answer = inputdlg(prompt,title,dims,definput) ;
            this.Nozzle.d = str2num(answer{1});
            
            % Initialize slug fields.
            this.Slugs = struct('SlugwiseTime',NaN,'SlugwiseMass',NaN,...
            'AverageSlugMass',NaN,'SlugwiseFlowrate',NaN,...
            'SlugwiseTotalMass',NaN,'SlugwiseTotalTime',NaN,...
            'Index',1,'TimeBeforeImpact',NaN,'MassBeforeImpact',NaN);
            % Initialize filtered data
            this.FilteredData = struct('Mass',NaN,'Time',NaN);
            % Initialize process parameters
            this.ProcessParameters = struct('ExpYieldStress',NaN,'ExpYieldStressError',NaN,'Flowrate',NaN);
       
        end
        
        function this = FilterData(this,rawData,rawDataRate)

    % Averaging kernel filter, to be convoluted with input signal.
    filter = ones( 1, this.KernelSize ) / this.KernelSize ;
    
    % Construct a filter doing averaging (as the precedent) and
    % derivation.
    if mod(this.KernelSize,2) == 0
            dfilter = [ ones( 1, this.KernelSize ), -ones( 1, this.KernelSize ) ] / ( 2*this.KernelSize ) ;
        else
            dfilter = [ ones( 1, floor(this.KernelSize) ), 0, -ones( 1, floor(this.KernelSize) ) ] / (2*this.KernelSize) ;
    end
        
     % Apply averaging filter to input data.
    m_conv = conv( rawData.Mass , filter, 'valid' ) ;
    t_conv = conv( rawData.Time, filter, 'valid' ) ;
    
    % Sample averaged data
    m_samp = m_conv( 1 : floor( rawDataRate/this.Rate ) : end ) ;    
    t_samp = t_conv( 1 : floor( rawDataRate/this.Rate ) : end ) ;
    
    % Build index list out of sampled-averaged mass data.
    id_samp = 1:length( m_samp ) ;
                
    % Apply averaging-derivation filter to raw input data.
    dm_conv = conv( rawData.Mass , dfilter, 'valid' ) ;
    % Sample averaging-derivation filtered data.
    dm_samp = dm_conv( 1 : floor( rawDataRate/this.Rate ) : end ) ;
                
    % Find threshold.
    threshold = var(dm_samp)^.5 * 2 ;
    
    % Detect drips
    id_impact = id_samp ( and( dm_samp(2:end) > threshold, dm_samp(1:end-1) < threshold ) ) ;
    
    
    if ~isempty(id_impact)
    % get timestamp at the index of detected drip
    t_impact = t_samp(id_impact+1) ;
                
    m_impact = zeros( 1,length(id_impact) )  ;
    f = this.Rate;
    
    % Criterion to avoid something about under-sampling ?  
    if id_impact(1)-floor(0.2*f )<=0 
        minT = 1 ;
    else 
        minT = id_impact(1)-floor(0.2*f ) ;
    end
                
    % Get some mass data points before impact. How much is
    % based on sampling frequency and an arbitrary coefficient
    m_av_impact = [m_samp( minT:(id_impact(1) - floor(0.05*f )) ),NaN] ;
    t_av_impact = [t_samp( minT:(id_impact(1) - floor(0.05*f )) ),NaN] ;
                
    % Average masses of previous weights data points before
    % impact
    m_impact(1) = mean( m_samp( minT:(id_impact(1) - floor(0.05*f )) ) ) ;
    for k = 2:length(id_impact)
     pond = 0.5  ;
     id_mass = [ floor( ( id_impact(k-1)*pond + id_impact(k)*(1-pond) ) ), id_impact(k)- floor(0.05*f )] ;
     m_impact(k) = mean( m_samp( id_mass(1):id_mass(2) ) ) ;
    end
    
    % Compute results into short-named variables so its readable
    
    slugwiseTime = diff( t_impact ) ;
    slugwiseMass = diff( m_impact ) ;
    
    averageSlugMass = mean(slugwiseMass(~isnan(slugwiseMass)));

    slugwiseFlowrate = slugwiseMass ./ slugwiseTime ;
        
    slugwiseTotalMass = m_impact;
    slugwiseTotalTime = t_impact;
    index = id_impact;
    timeBeforeImpact = t_av_impact;
    massBeforeImpact = m_av_impact;

    % Return filtered data.
    filteredMass = m_samp';
    filteredTime = t_samp';
    %m_derivative = dm_samp';
    
    % Apply short-named variables to long-named output variables
    
    this.Slugs.SlugwiseTime = slugwiseTime;
    this.Slugs.SlugwiseMass = slugwiseMass;
    this.Slugs.AverageSlugMass = averageSlugMass;
    this.Slugs.SlugwiseFlowrate = slugwiseFlowrate;
    this.Slugs.SlugwiseTotalMass = slugwiseTotalMass;
    this.Slugs.SlugwiseTotalTime = slugwiseTotalTime;
    this.Slugs.Index = index;
    this.Slugs.TimeBeforeImpact = timeBeforeImpact;
    this.Slugs.MassBeforeImpact = massBeforeImpact;
    this.FilteredData.Mass = filteredMass;
    this.FilteredData.Time = filteredTime;
    
    else
        % return the expected structs but with NaN fields
        this.Slugs = struct('SlugwiseTime',NaN,'SlugwiseMass',NaN,...
            'AverageSlugMass',NaN,'SlugwiseFlowrate',NaN,...
            'SlugwiseTotalMass',NaN,'SlugwiseTotalTime',NaN,...
            'Index',1,'TimeBeforeImpact',NaN,'MassBeforeImpact',NaN);

        this.FilteredData = struct('Mass',NaN,'Time',NaN);
    %m_derivative = NaN;
        
    end
    
    
    
            
        end
        
        function processParameters = get.ProcessParameters(this)
        processParameters.ExpYieldStress = ExpYieldStress(this);
        processParameters.ExpYieldStressError = ExpYieldStressError(this);
        processParameters.Flowrate = DripFlowrate(this);
        end
        
                
        % Calculation subfunctions
        
        function [coefficients, rSquare] = LinearRegression(this, x,y)

% Flip matrices
x = x';
y = y';
% Detect NaNs.
xNan = isnan(x);
yNan = isnan(y);

% Remove NaNs.
try
xCleaned = x(xNan == 0,:);
yCleaned = y(yNan == 0,:);
catch
xCleaned = x(xNan == 0);
yCleaned = y(yNan == 0);    
end

% Remap to fit matrix dimensions
if length(xCleaned) > length(yCleaned)
xCleaned = xCleaned(1:length(yCleaned));
end
if length(yCleaned) > length(xCleaned)
yCleaned = yCleaned(1:length(xCleaned));
end

% Calculate b by padding x with a column of ones and using the \ operator. 
X = [ones(length(xCleaned),1) xCleaned];

% Do the linear regression
coefficients = X\yCleaned;

% Resulting slope.
slope = coefficients(2);

% Calculate R2 of fit.
yCalc = X*coefficients;
rSquare = 1 - sum((yCleaned - yCalc).^2)/sum((yCleaned - mean(yCleaned)).^2);
end

        
    end
    methods (Access = private)
        
        % Process parameters subfunctions
        function expYieldStress = ExpYieldStress(this)
        expYieldStress = mean(this.Slugs.AverageSlugMass(~isnan(this.Slugs.AverageSlugMass))) * 9.81 / (sqrt(3)*pi*this.Nozzle.d^2/4) * 1000 ;
        
        end
        function expYieldStressError = ExpYieldStressError(this)
        expYieldStressError = var( this.Slugs.AverageSlugMass(~isnan(this.Slugs.AverageSlugMass)) * 9.81 / (sqrt(3)*pi*this.Nozzle.d^2/4) * 1000 )^.5 *...
                    2 / sqrt(length(this.Slugs.AverageSlugMass) ) ;    
        end
        function dripFlowrate = DripFlowrate(this)
            
            y = this.Slugs.SlugwiseTotalMass;
            x = this.Slugs.SlugwiseTotalTime;
            
            if y < 0
                y = zeros(1:length(y));
            end
            
            if sum(~isnan(x))>0 || sum(~isnan(y))>0
            [coefficients,~] = this.LinearRegression(x,y);
            dripFlowrate = coefficients(2);
            else
            dripFlowrate = NaN;
            end
        end

        
    end
end

