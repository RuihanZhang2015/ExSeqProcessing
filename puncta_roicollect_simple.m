loadParameters;

filename_centroids = fullfile(params.punctaSubvolumeDir,sprintf('%s_centroids+pixels.mat',params.FILE_BASENAME));
load(filename_centroids,'puncta_centroids','puncta_voxels')
%% Collect the subvolumes we started this with, but now only with the pixels from the puncta!

num_insitu_transcripts = size(puncta_voxels,1);

%Define a puncta_set object that can be parallelized
puncta_set_cell = cell(params.NUM_ROUNDS,1);


%the puncta indices are here in linear form for a specific round

parpool(5); %arbitrary but this parallel loop is memory intensive
parfor exp_idx = 4:params.NUM_ROUNDS
    disp(['round=',num2str(exp_idx)])
    pixels_per_rnd = cell(num_insitu_transcripts,params.NUM_CHANNELS);
    
    for c_idx = params.COLOR_VEC
        filename_in = fullfile(params.registeredImagesDir,sprintf('%s_round%.03i_%s_%s.tif',params.FILE_BASENAME,exp_idx,params.CHAN_STRS{c_idx},regparams.REGISTRATION_TYPE));
        img =  load3DTif_uint16(filename_in);
        
        for puncta_idx = 1:num_insitu_transcripts
            
            indices_for_puncta = puncta_voxels{puncta_idx};
            
            %This makes a volume that is all zeros except for the punctafeinder_indices_for_puncta
            pixels_for_puncta_set = img(indices_for_puncta);
            
            %Then we take the PUNCTA_SIZE region around those pixels only
            pixels_per_rnd{puncta_idx,c_idx} = pixels_for_puncta_set;
            
            if mod(puncta_idx,5000)==0
                fprintf('Rnd %i, Chan %i, Puncta %i processed\n',exp_idx,c_idx,puncta_idx);
            end
            
        end
        
    end
    puncta_set_cell{exp_idx} = pixels_per_rnd;

end

disp('reducing processed puncta')
puncta_set_median = zeros(params.NUM_ROUNDS,params.NUM_CHANNELS,num_insitu_transcripts);
% reduction of parfor
for puncta_idx = 1:num_insitu_transcripts
    for exp_idx = 1:params.NUM_ROUNDS
        for c_idx = params.COLOR_VEC
            % Each puncta_set_cell per exp is
            % pixels_per_rnd = cell(num_insitu_transcripts,params.NUM_CHANNELS);
            puncta_set_median(exp_idx,c_idx,puncta_idx) = median(puncta_set_cell{exp_idx}{puncta_idx,c_idx});
            %             puncta_set_mean(exp_idx,c_idx,puncta_idx) = mean(puncta_set_cell{exp_idx}{puncta_idx,c_idx});
            %             puncta_set_max(exp_idx,c_idx,puncta_idx) = max(puncta_set_cell{exp_idx}{puncta_idx,c_idx});
        end
    end
end

% Using median values to ignore bad points
channels_present = squeeze(sum(puncta_set_median==0,2));
num_roundspresent_per_puncta = squeeze(sum(channels_present==0,1));
%Create a mask of all puncta that have a non-zero signal in all four
%channels for all rounds
signal_complete = num_roundspresent_per_puncta==params.NUM_ROUNDS;

fprintf('Number of complete puncta: %i \n',sum(signal_complete));

puncta_set_median = puncta_set_median(:,:,signal_complete);
puncta_centroids = puncta_centroids(signal_complete,:);

puncta_size = zeros(sum(signal_complete),1);
indices_of_good_puncta = find(signal_complete);
for puncta_idx = 1:size(puncta_centroids,1)
    puncta_size(puncta_idx) = length(puncta_voxels{indices_of_good_puncta(puncta_idx)});
end
%just save puncta_set
% outputfile = '/mp/nas0/ExSeq/AutoSeq2/xy10/6_transcripts/exseqauto-xy10_puncta_noncubepixels';
outputfile = fullfile(params.transcriptResultsDir,sprintf('%s_puncta_noncubepixels.mat',params.FILE_BASENAME));
save(outputfile,'puncta_set_median','puncta_centroids','puncta_size','-v7.3');

%% Analyze the output
% Using median values to ignore bad points

% figure;
% plot3(puncta_centroids(~signal_complete,1),puncta_centroids(~signal_complete,2),...
%     puncta_centroids(~signal_complete,3),'b.');
% hold on;
% plot3(puncta_centroids(signal_complete,1),puncta_centroids(signal_complete,2),...
%     puncta_centroids(signal_complete,3),'rx','MarkerSize',5);

