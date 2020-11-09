
% Load example sensor data.
load('slugtest_example/DataExample.mat');

% Instantiate slug detection filter.
filter = SlugDerivativeFilter();

% Filter raw sensor data to detect slugs.
filter = filter.FilterData(rawData,rawData.Rate);

% Get filtering results from filter object.
slugs = filter.Slugs; 
filteredData = filter.FilteredData;
processParameters = filter.ProcessParameters;

% Plot raw data
plot(rawData.Time,rawData.Mass);
% Plot filtered data
plot(filteredData.Time,filteredData.Mass);
% Plot detected slugs.
scatter(slugs.SlugwiseTotalTime,slugs.SlugwiseTotalMass);

% Display process parameters.
disp(strcat('Apparent Yield Stress: ',num2str(processParameters.ExpYieldStress),' Pa'));
disp(strcat('Apparent Yield Stress Error: ',num2str(processParameters.ExpYieldStressError),' Pa'));
disp(strcat('Apparent Flowrate (based on detected slugs): ',num2str(processParameters.Flowrate),' g/s'));