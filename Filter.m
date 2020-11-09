classdef Filter
    %FILTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    FilteredData
    Rate
    KernelSize
    end
    
    methods
        function this = Filter(rawData)
% Link filter to the raw data.

            prompt = {'Enter sampling rate (Hz)','Enter kernel size'};
            title = 'Filter properties';
            dims = [1 35];
            definput = {'25','800'} ;
            answer = inputdlg(prompt,title,dims,definput) ;
            this.Rate = str2double( answer{1} );
            this.KernelSize = str2double( answer{2} ) ;
        end
    end
end

