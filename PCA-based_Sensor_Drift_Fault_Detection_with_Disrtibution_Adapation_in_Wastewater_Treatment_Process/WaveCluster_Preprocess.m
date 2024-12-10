%  The computational process function for wavelet clustering 

function [sigcells, datacellindices, counts, wdata] = WaveCluster_Preprocess(data, weights, num_cells, densitythreshold, level, wavename, useSWT)

    if (~exist('num_cells', 'var') || isempty(num_cells))
        num_cells = max(data) - min(data) + 1; 
    end
    
    if (~exist('densitythreshold', 'var') || isempty(densitythreshold))
        densitythreshold = '10%';
    end
    
    if (~exist('level', 'var') || isempty(level))
        level = 1;
    end

    if (~exist('wavename', 'var') || isempty(wavename))
        wavename = 'bior2.2';
    end

    if (~exist('useSWT', 'var') || isempty(useSWT))
        useSWT = 0; 
    end
    

    [counts, datacellindices] = data2grid(data, weights, num_cells);

   
    oldmode = dwtmode('status', 'nodisp');
    dwtmode('zpd', 'nodisp'); 
    
    %[decwaves, bookkeep] = wavedec2(counts, level, wavename);
    try
        if (useSWT)
            wdata = swtN(counts, level, wavename);
        else
            wdata = dwtN(counts, level, wavename);

           
            datacellindices = ceil(datacellindices ./ 2 ^ level);
        end
    catch me
        dwtmode(oldmode, 'nodisp');   
        rethrow(me);
    end
    
    dwtmode(oldmode, 'nodisp');   
    
   
    if (ischar(densitythreshold) && densitythreshold(end) == '%')
       
        pctthresh = str2double(densitythreshold(1:end-1));
        densitythreshold = prctile(wdata(wdata > 0), pctthresh);
%         disp(['Automatic density threshold @ ' num2str(pctthresh) '%: ' num2str(densitythreshold)])
    end

   
    sigcells = (wdata >= densitythreshold);
end
