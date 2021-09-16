%% General
datadir = '/host/fladgate/local_raid/HCP_data/structuralData/mpc/myelinMap/';
subjList = []; % enter your subject list here
numSurf = 12; % or 14
nVertices = 20484; % fsa5 surface

% Load parcellation
parc = 'schaefer-200';
[vert_lh, label_lh, ctb_lh] = read_annotation(['/data_/mica1/03_projects/jessica/micapipe/parcellations/lh.', parc,'_mics.annot']);
[vert_rh, label_rh, ctb_rh] = read_annotation(['/data_/mica1/03_projects/jessica/micapipe/parcellations/rh.', parc,'_mics.annot']);
parcFS = zeros(1,length(label_lh));
for ii = 1:length(label_lh)
    parcFS(ii) = find(ctb_lh.table(:,5) == label_lh(ii));
end
for ii = 1:length(label_rh)
    parcFS(ii+length(label_lh)) = ...
        find(ctb_rh.table(:,5) == label_rh(ii)) + length(ctb_rh.table(:,5));
end

% Mask
maskIdx_lh = 1; maskIdx_rh = 102;
mask = ones(1,nVertices);
mask(parcFS == maskIdx_lh) = 0; mask(parcFS == maskIdx_rh) = 0; 
mask = logical(mask);

parcFS_mask = parcFS(mask);


%% Build subject intensity profiles
profiles = zeros(numSurf, nVertices, length(subjList));
for subj = 1:length(subjList)
    this_subj = subjList(subj);
    
    for surf = 1:numSurf
        lh = SurfStatReadData([dataDir, this_subj, '.lh.', surf, '.fsaverage5.mgh']);
        rh = SurfStatReadData([dataDir, this_subj, '.rh.', surf, '.fsaverage5.mgh']);
        concat = [lh';rh'];

        % If the surfaces are not smoothed, you can use SurfStatSmooth to
        % smooth them. FWHM = 2 or 3 is usually good. I would check with
        % casey though as she might have smoothed them already when mapping
        % the native surface data to fsa5
        % concat_smooth = SurfStatSmooth(concat', surf, FWHM)
        
        profiles(surf,:,subj) = concat; % or concat_smooth
    end
end

% Mean profile of all subjects
profilesM = mean(profiles,3);
save('/data/mica1/03_projects/yezhou/project1/data_HCP/mpc/HCP_MPCprofiles.mat','profiles');
%% Build MPC

% MPC function path:
% /micaopen/MPC/scripts/build_mpc.m 

MPCall = [];
for subj = 1:length(subjList)
    this_profile = profiles(:,:,subj);
    
    % Remove mask from profiles and parcellation
    profileMask = this_profile(:,mask);
    
    % Build MPC; I = parcellated intensity profiles
    [MPC, I] = build_mpc(profileMask, parcFS_mask);
    
    MPCall(:,:,subj) = MPC;
    
end
MPC_group=mean(MPCall,3);
save('/data/mica1/03_projects/yezhou/project1/data_HCP/mpc/HCP_MPCall.mat','MPCall');

%% Build MPC gradient
myelinG1_REF    = GradientMaps('n_components', 13);
myelinG1_REF    = myelinG1_REF.fit(MPC_group);
Gref            = myelinG1_REF.gradients{1};

MPCGradient_all_aligned=zeros(200,13,length(subjList));
for subj = 1:length(subjList)
    myelinG1    = GradientMaps('n_components', 13,'alignment','pa');
    myelinG1    = myelinG1.fit(MPCall(:,:,subj), 'reference',Gref);
    MPCGradient_all_aligned(:,:,subj)    = myelinG1.aligned{1}(:, 1:13);
    
end
MPCGradient_group=mean(MPCGradient_all_aligned,3);
MPCG1_group=MPCGradient_group(:,1);

Gs=zeros(202,1); Gs([1 102],:)=-inf;
Gs(2:101,:)=MPCG1_group(1:100,1); Gs(103:202,:)=MPCG1_group(101:200,1);
A=plot_hemispheres(Gs,{surf_l, surf_r}, ...
    'parcellation', parcFSNew, 'labeltext',{'HCP_MPCG1'});
% A.colorlimits([5.71,13.5271]);
A.colormaps([0.8 0.8 0.8; parula]);

save('/data/mica1/03_projects/yezhou/project1/data_HCP/mpc/HCP_MPCGradient_all_aligned.mat','MPCGradient_all_aligned');

