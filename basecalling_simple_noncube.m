

load(fullfile(params.transcriptResultsDir,sprintf('%s_puncta_noncubepixels.mat',params.FILE_BASENAME)));
% 
if ~exist('gtlabels','var')
    load('groundtruth_dictionary_slice.mat');
end
%load('groundtruth_dictionary_neurons.mat')
%% normalize all puncta intensities by their Z value


%Create a new variable, puncta_set, which is the cropped puncta_set_median
%for only the bases that we want to call. Given the first three bases are
%magenta, this could screw things up.
ROUNDS_TO_CALL = [4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 20];
NUMCALLEDROUNDS = length(ROUNDS_TO_CALL);

puncta_set_median_znormalized = zeros(size(puncta_set_median));

max_z = max(floor(puncta_centroids(:,3)));
% Here we 
for rnd_idx = 1:20
    
    means_z = zeros(max_z,4);
    stds_z = zeros(max_z,4);
    
    for z = 1:max_z
        indices = find(floor(puncta_centroids(:,3))==z);
        if length(indices)<=1
            puncta_set_median_znormalized(rnd_idx,:,indices)=...
                puncta_set_median(rnd_idx,:,indices);
            continue 
        end
        %Get the statistics for all the puncta with the center at each Z
        for c = 1:4
            means_z(z,c) = mean(puncta_set_median(rnd_idx,c,indices));
            stds_z(z,c) = std(puncta_set_median(rnd_idx,c,indices));
        end
        
        %Aply the normalization
        for c = 1:4
            p_nonnorm = squeeze(puncta_set_median(rnd_idx,c,indices));
            
            if stds_z(z,c) ~= 0
                p_norm = (p_nonnorm-means_z(z,c))/stds_z(z,c);
            else %If all the puncta for this chan at this Z are the same
                p_norm = 0; %this is a rare case but worth catching
            end
            if any(isnan(p_norm(:)))
               barf() 
            end
            puncta_set_median_znormalized(rnd_idx,c,indices)=p_norm;
        end
    end
    
    

   
    fprintf('Z-normalized Round%i\n',rnd_idx);
end

%%  Convert all the data into zscores (very cheap base calling)
num_puncta = size(puncta_centroids,1);

%Pre-initialize the cell arrray and determine the basecalls

base_calls_quickzscore = zeros(num_puncta,params.NUM_ROUNDS);
base_calls_quickzscore_confidence = zeros(num_puncta,params.NUM_ROUNDS);
base_calls_pixel_intensity = zeros(num_puncta,params.NUM_ROUNDS,params.NUM_CHANNELS);
base_calls_zscores = zeros(num_puncta,params.NUM_ROUNDS,params.NUM_CHANNELS);

puncta_set_roundnormed = zeros(size(puncta_set_median_znormalized));

for rnd_idx = 1:NUMCALLEDROUNDS
    
    actual_rnd_idx = ROUNDS_TO_CALL(rnd_idx);
    clear chan_col; %Just in case, otherwise the for loop can error.
    for c = params.COLOR_VEC
        chan_col(:,c) = reshape(puncta_set_median_znormalized(actual_rnd_idx,c,:),[],1);
    end
    
    % cols_normed = quantilenorm(chan_col);
    cols_normed = zscore(chan_col);
    
    for c = params.COLOR_VEC
        puncta_set_roundnormed(actual_rnd_idx,c,:) = reshape(cols_normed(:,c),size(squeeze(puncta_set_median_znormalized(actual_rnd_idx,c,:))));
    end
    
    fprintf('Looping through puncta from round %i\n',rnd_idx);
    for p_idx= 1:num_puncta
        
        
        scores = squeeze(puncta_set_roundnormed(actual_rnd_idx,:,p_idx));
        
        %add the minimum score to shift everything non-zero
        scores = scores - min(scores);
        
        %and the new baseguess
        [~, newbaseguess] = max(scores);
        base_calls_quickzscore(p_idx,actual_rnd_idx) = newbaseguess;
        
        [scores_sorted,~] = sort(scores,'descend');
        base_calls_quickzscore_confidence(p_idx,actual_rnd_idx) = scores_sorted(1)/(scores_sorted(1)+ scores_sorted(2));
        
        %use the baseguess to get the absolute brightness of the puncta
        base_calls_pixel_intensity(p_idx,actual_rnd_idx,:) = puncta_set_median(rnd_idx,:,p_idx);
        
        base_calls_zscores(p_idx,actual_rnd_idx,:) = scores;
    
    end
    
end

[unique_transcipts,~,~] = unique(base_calls_quickzscore,'rows');
fprintf('Found %i transcripts, %i of which are unique\n',size(base_calls_quickzscore,1),size(unique_transcipts,1));

%% Do a quick visualization of the base calls across rounds
for base_idx = 1:params.NUM_CHANNELS
    perc_base(:,base_idx) = sum(base_calls_quickzscore==base_idx,1)/size(base_calls_quickzscore,1);
end
figure;
% Chan 1 = Blue
% Chan 2 = Green
% Chan 3 = Magenta
% Chan 4 = Red

plot(perc_base(:,1)*100,'b','LineWidth',2); hold on;
plot(perc_base(:,2)*100,'g','LineWidth',2)
plot(perc_base(:,3)*100,'m','LineWidth',2)
plot(perc_base(:,4)*100,'r','LineWidth',2); hold off;
legend('Chan1 - FITC','Chan2 - CY3', 'Chan3 - Texas Red', 'Chan4 - Cy5');
title(sprintf('Percentage of each base across rounds for %i puncta',size(base_calls_quickzscore,1)));


%% Make sets of transcripts and create a new transcript object
%The ground truth starts on Round4, so subtract 3 to get in alignment
gt_mask = ROUNDS_TO_CALL-3;
transcript_objects = cell(size(base_calls_quickzscore,1),1);
output_cell = cell(size(base_calls_quickzscore,1),1);
parfor p_idx = 1:size(base_calls_quickzscore,1)
    
    transcript = struct;
    %Search for a perfect match in the ground truth codes
    img_transcript = base_calls_quickzscore(p_idx,ROUNDS_TO_CALL);
    
    %Sanity check: randomize the img_transcript
    %     img_transcript = img_transcript(randperm(length(img_transcript)));
    
    %Search for a perfect match in the ground truth codes
    hits = (groundtruth_codes(:,gt_mask)==img_transcript);
    
    %Calculate the hamming distance (now subtracting primer length)
    scores = length(img_transcript)- sum(hits,2);
    [values, indices] = sort(scores,'ascend');
    
    best_score = values(1);
    %Get the first index that is great than the best score
    %If this is 1, then there is a unique fit
    %     idx_last_tie = find(values>best_score,1)-1;
    idx_last_tie = 1;
    while values(idx_last_tie)<=best_score %Trying a faster method than find()
        idx_last_tie = idx_last_tie +1;
    end
    idx_last_tie = idx_last_tie-1;
    
    intensities = squeeze(base_calls_pixel_intensity(p_idx,:,:));

    %Assuming the groundtruth options are de-duplicated
    %Is there a perfect match to the (unique ) best-fit
    transcript.img_transcript=base_calls_quickzscore(p_idx,:);
    transcript.img_transcript_confidence=base_calls_quickzscore_confidence(p_idx,:);
    transcript.img_transcript_absValuePixel=intensities;
    transcript.img_transcript_zscores=squeeze(base_calls_zscores(p_idx,:,:));
    transcript.pos = puncta_centroids(p_idx,:);
    transcript.hamming_score = best_score;
    transcript.numMatches = idx_last_tie;
    
%     transcript.nn_distance = spacings(p_idx);
    
    %If the alignment to known sequences is one error or less
    if best_score <=1
        row_string = sprintf('%i,%s,',p_idx,mat2str(img_transcript));
        
        %If it's a unique match, use the name
        if idx_last_tie==1
            transcript.known_sequence_matched = groundtruth_codes(indices(idx_last_tie),:);
            transcript.name = gtlabels{indices(idx_last_tie)};
        end
        
        for tiedindex = 1:idx_last_tie
            row_string = strcat(row_string,sprintf('%s,',gtlabels{indices(tiedindex)}));
        end
        row_string = sprintf('%s\n',row_string);
        %fprintf('%s',row_string);
        output_cell{p_idx} = row_string; 
    end
    
    transcript_objects{p_idx} = transcript;
    
    if mod(p_idx,1000) ==0
        fprintf('%i/%i matched\n',p_idx,size(base_calls_quickzscore,1));
    end
end
fprintf('Done!\n');

output_cell = output_cell(~cellfun('isempty',output_cell));  
output_csv = strjoin(output_cell,'');

output_file = fullfile(params.transcriptResultsDir,sprintf('%s_simpleextractedcodes_noncube.csv',params.FILE_BASENAME));

fileID = fopen(output_file,'w');
fprintf(fileID,output_csv);
fclose(fileID);

%%
%remove the transcripts that had an empty round, likely because of the registration

save(fullfile(params.transcriptResultsDir,sprintf('%s_transcriptmatches_objects_noncube.mat',params.FILE_BASENAME)),'transcript_objects','ROUNDS_TO_CALL','-v7.3');
fprintf('Saved transcript_matches_objects!\n');

