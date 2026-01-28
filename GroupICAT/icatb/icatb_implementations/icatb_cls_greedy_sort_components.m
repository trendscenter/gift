classdef icatb_cls_greedy_sort_components < handle
    % Date 9/30/25
    % Cyrus Eierud
    % Code that greedily sorts the best correlated matches between your
    % components and a template
    % 1/27/26 Added function to sort between two arbitrary files.
    % Example to calcuate and save state guided ICA:
    %    oc_sort = icatb_cls_greedy_sort_components(sesInfo); %sesInfo structure holds param file info
    %    oc_sort.m_dialog % engages greedy sort    
    
    properties
        stru_sesInfo % parameter file
        s_mask % mask used to remove non brain voxels
        ron_mask_ix  % mask used to remove non brain voxels
        s_components2compare % your blind components
        s_template %template to compare against
        cob_components_filter_inclusion % components you want to keep as 1 and filter out as 0
        ard_corrs_table
        ari_ordered_pairs_table
        ari_ordered_pairs_orig %The order number for both s_components2compare and the template
    end
    
    methods
        function n_err = m_calc_greed_match(o)
            % sorts components between your components and a template
            
            n_err = 1;

            disp(['icatb_info [' char(datetime) '] cls_greedy_sort_components.m: Starting Greedy Sort'])

            if isempty(o.s_mask)
                %Get the mask variable
                ixMas = o.ron_mask_ix;
            else
                % Mask file exists - read it in and Flatten mask
                [VV, HInfo2] = icatb_returnHInfo(o.s_mask); 
                icasig2 = icatb_spm_read_vols(VV);
                structDIM2 = HInfo2.DIM;
                icasig2 = reshape(icasig2, 1, prod(structDIM2));
                icasig2(isnan(icasig2))=0; %clean nans
                icasigMask=icasig2;
                clear icasig2;
                tol = 1e-2;
                ixMas=find(abs(icasigMask - 1) < tol);
            end

            [o.ard_corrs_table, o.ari_ordered_pairs_table, o.ari_ordered_pairs_orig] = mpr_correlation(o, ixMas, o.s_components2compare, o.cob_components_filter_inclusion, o.s_template);

            save('o');

            disp(['icatb_info [' char(datetime) '] cls_greedy_sort_components.m: Completed Greedy Sort'])
            disp(['icatb_info [' char(datetime) '] cls_greedy_sort_components.m: Components sorted from file ' o.s_components2compare])
            disp(['icatb_info [' char(datetime) '] cls_greedy_sort_components.m: Template sorted against ' o.s_template])
            mkdir(fullfile(fileparts(o.s_components2compare),'utilities'));
            save(fullfile(fileparts(o.s_components2compare),['utilities' filesep 'greedsort' datestr(datetime('now'),'yyyymmddHHMMSS') '.mat']), 'o');
            disp(['icatb_info [' char(datetime) '] cls_greedy_sort_components.m: cls_greedy_sort_components.m: Greedy sort component pairs across files, saved in var o.ari_ordered_pairs_table under ' fullfile(fileparts(o.s_components2compare),['utilities' filesep 'greedsort' datestr(datetime('now'),'yyyymmddHHMMSS') '.mat'])])
            
            n_err = 0;
        end      

        function o = icatb_cls_greedy_sort_components(stru_sesInfo)
            % Constructor
            if isempty(stru_sesInfo)
                o.s_mask = ''; 
                o.ron_mask_ix = [];
                o.s_components2compare = [];
            else
                % needs to be the sesinfo variable in stru_sesInfo
                o.stru_sesInfo = stru_sesInfo;
                o.s_mask = fullfile(o.stru_sesInfo.outputDir, [o.stru_sesInfo.userInput.prefix, 'Mask.nii']); 
                o.ron_mask_ix = [];
                o.s_components2compare = fullfile(o.stru_sesInfo.outputDir, o.stru_sesInfo.icaOutputFiles(1).ses.name);                 
            end

        end

        function val = get_s_mask(o)
            % Method to get Value
            val = o.s_mask;
        end

        function val = get_s_components2compare(o)
            % Method to get Value
            val = o.s_components2compare;
        end  

        function o = set_cob_components_filter_inclusion(o, cob_new)
            % Method to set Value

            %check that the incoming cob_web has same length as the
            %components
            o.cob_components_filter_inclusion = cob_new;
        end

        function val = get_cob_components_filter_inclusion(o)
            % Method to get Value
            val = o.cob_components_filter_inclusion;
        end

        function o = set_s_template(o, s_new_file)
            % Method to set Value
            o.s_template = s_new_file;
        end

        function val = get_s_template(o)
            % Method to get Value
            val = o.s_template;
        end

        function h_popup = m_dialog(o)

            % set up the defaults
            icatb_defaults;
            global BG_COLOR;
            global BG2_COLOR;
            global BUTTON_COLOR;
            global FG_COLOR;
            global AXES_COLOR;
            global BUTTON_FONT_COLOR;
            
            global UI_FONTNAME;
            global UI_FONTUNITS;
            global FONT_COLOR;
            global UI_FS;
            
            % set up fonts
            titleFont = 14;
            textFont = 12;
            
            % Position of the controls
            axisPos = [0 0 1 1];
            titlePos = [0.5 0.95];
            textPos = [0.1 0.65];
            xOffSet = 0.01; yOffSet = 0.04;
            
            % title color
            titleColor = [0 0.9 0.9];

            % set the defaults for the figure window
            figHandle = figure('Resize','off', ...
                'menubar', 'none', ...
                'DefaultTextColor', FONT_COLOR,...
                'DefaultTextInterpreter', 'none',...
                'DefaultAxesColor', AXES_COLOR,...
                'DefaultAxesXColor', 'k',...
                'DefaultAxesYColor', 'k',...
                'DefaultAxesZColor', 'k',...
                'DefaultPatchFaceColor', 'k',...
                'DefaultPatchEdgeColor', 'k',...
                'DefaultSurfaceEdgeColor', 'k',...
                'DefaultLineColor', 'k',...
                'DefaultUicontrolInterruptible', 'on',...
                'PaperType', 'usletter',...
                'PaperUnits', 'normalized',...
                'PaperPositionMode', 'auto', ...
                'InvertHardcopy', 'off',...
                'Renderer', 'zbuffer',...
                'color', BG_COLOR, 'resize', 'off', 'name', 'Greedy Search', ...
                'visible', 'off');
            
            % change the size of the figure window here
            pos = get(figHandle, 'position');
            screenSize = get(0, 'screensize');
            figurePos(1) = 0.5*(screenSize(1) + screenSize(3)) - 0.5*pos(3);
            figurePos(2) = 0.5*(screenSize(2) + screenSize(4)) - 0.5*pos(4);
            figurePos(3) = 0.95*pos(3);
            figurePos(4) =  0.7*pos(4);
            
            set(figHandle, 'position', figurePos);
            
            % set axis handle to off
            axisHandle = axes('Parent', figHandle, 'Position', axisPos, 'Visible', 'off');
            
            % Name of the toolbox
            text('units', 'normalized', 'string', 'Greedy Sort of Components vs a Template', 'position', titlePos, 'fontsize', titleFont, 'HorizontalAlignment', 'center', ...
                'fontweight', 'bold', 'FontName', UI_FONTNAME, 'color', titleColor, 'FontAngle', 'italic');
            
            textPos(2) = titlePos(2) - 0.15-0.32;

            % label
            ui_template_pos(3) = 0.8; ui_template_pos(4) = 0.40;
            ui_template_pos(1) = 0.1; ui_template_pos(2) = textPos(2);
            icatb_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'text', 'position', ui_template_pos, ...
                'string', 'After selecting parameter file pointing to your custom blind components, you may follow steps 1,2 and 3. Step 1 selects template. Optional step 2, inserting comma separated booleans, 1 for signal and 0 for noisy components (Autolabeller may discriminate). After step 3 (run), check the saved results file.', 'tag', 'tag_lab', 'fontsize', UI_FS - 1, 'HorizontalAlignment', 'left');

            textPos(2) = textPos(2) - 0.11;

            % label
            ui_template_pos(3) = 0.5; ui_template_pos(4) = 0.08;
            ui_template_pos(1) = 0.1; ui_template_pos(2) = textPos(2);
            icatb_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'text', 'position', ui_template_pos, ...
                'string', '1. NMark Template', 'tag', 'tag_lab', 'fontsize', UI_FS - 1, 'HorizontalAlignment', 'left');

            % template combobox
            stru_files_templates = dir([fullfile(fileparts(which('gift.m')),'icatb_templates/Neuromark*.nii')]);
            ui_template_pos(3) = 0.5; ui_template_pos(4) = 0.08;
            ui_template_pos(1) = 0.4; ui_template_pos(2) = textPos(2);
            icatb_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'popup', 'position', ui_template_pos, ...
                'string', str2mat(string({stru_files_templates(:).name})), ... 
                'tag', 'tag_ui_pop_templ', 'fontsize', UI_FS - 1);


            textPos(2) = textPos(2) - 0.15;

            % label
            ui_template_pos(3) = 0.5; ui_template_pos(4) = 0.08;
            ui_template_pos(1) = 0.1; ui_template_pos(2) = textPos(2);
            icatb_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'text', 'position', ui_template_pos, ...
                'string', '2. (Signal[1]/Noise[0] Vec.)', 'tag', 'tag_lab', 'fontsize', UI_FS - 1, 'HorizontalAlignment', 'left');

            % text autolabeller position
            ed_auto_pos(3) = 0.4; ed_auto_pos(4) = 0.08;
            ed_auto_pos(1) = 0.5; ed_auto_pos(2) = textPos(2);
            icatb_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'edit', 'position', ed_auto_pos, 'string', '(optionally replace with bool vector)', ... 
                'tag', 'tag_ui_ed_b_autolab_filt', 'fontsize', UI_FS - 1);

            textPos(2) = textPos(2) - 0.15;
            % label
            ui_template_pos(3) = 0.5; ui_template_pos(4) = 0.08;
            ui_template_pos(1) = 0.1; ui_template_pos(2) = textPos(2);
            icatb_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'text', 'position', ui_template_pos, ...
                'string', '3. Action:', 'tag', 'tag_lab', 'fontsize', UI_FS - 1, 'HorizontalAlignment', 'left');

            % run position
            rod_run_pos(3) = 0.12; rod_run_pos(4) = 0.08;
            rod_run_pos(1) = 0.5 - 0.5*rod_run_pos(3); rod_run_pos(2) = textPos(2);
            icatb_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'pushbutton', 'string', 'Run', ...
                'callback', @(src,evt)o.mpr_run_callback(src,evt,figHandle), 'position', rod_run_pos, 'fontunits', UI_FONTUNITS, 'fontname', UI_FONTNAME, ...
                'FontSize', UI_FS, 'ForegroundColor', BUTTON_FONT_COLOR, 'backgroundcolor', BUTTON_COLOR);

            % Initialise cancel position
            rod_cancel_pos = rod_run_pos;
            rod_cancel_pos(1) = 0.75 - 0.5*rod_run_pos(3);        
            % plot button
            icatb_uicontrol('parent', figHandle, ...
                'units', 'normalized', 'style', 'pushbutton', ...
                'string', 'Cancel', ...
                'callback', @(src,evt)o.mpr_cancel_callback(src, evt, figHandle), ...
                'position', rod_cancel_pos, ...
                'fontunits', UI_FONTUNITS, 'fontname', UI_FONTNAME, ...
                'FontSize', UI_FS, 'ForegroundColor', BUTTON_FONT_COLOR, ...
                'backgroundcolor', BUTTON_COLOR, 'tag', 'tag_cancel');

            set(figHandle, 'visible', 'on');
        end

        function arAllComp = m_greedy_simple(o, s_nii_1, s_nii_2)

            o.s_components2compare=s_nii_1;
            o.cob_components_filter_inclusion=1;
            o.s_template=s_nii_2;

            [ard_corrs_table, ari_ordered_pairs_table, ari_ordered_pairs_orig] = mpr_correlation(o,[], o.s_components2compare, 1, o.s_template);
            n_err = o.m_calc_greed_match();

        end

    end

    methods (Access = private)

        function mpr_cancel_callback(o, handleObj, evd, figHandle)
            delete(figHandle);
        end   

        function mpr_run_callback(o, handleObj, evd, figHandle)            

            % get template
            h_ui_pop_templ = findobj(figHandle, 'tag', 'tag_ui_pop_templ'); 
            i_val = get(h_ui_pop_templ,'Value');              % numeric index of selection
            ars_items = get(h_ui_pop_templ,'String');           % cell array of all items
            s_templ = strtrim(ars_items(i_val,:));             % the actual selected string
            o.s_template = [fullfile(fileparts(which('gift.m')),'icatb_templates/'), s_templ];

            % Data to match noise to 
            v_tmp = icatb_spm_vol(o.s_components2compare);
            n_tot_comps =  size(v_tmp,1);

            h_ui_ed_b_autolab_filt = findobj(figHandle, 'tag', 'tag_ui_ed_b_autolab_filt'); 
            s_rob_autolab_filt = get(h_ui_ed_b_autolab_filt, 'string');
            nums = split(s_rob_autolab_filt, {',',' ','\t','[',']'});
            nums = str2double(nums);
            if all(ismember(nums,[0 1]) & ~isnan(nums))
                vec = logical(nums);
                if n_tot_comps == size(vec,1)
                    o.cob_components_filter_inclusion = vec;
                    disp(['icatb_info [' char(datetime) '] cls_greedy_sort_components.m: componrent model order and bool vector matches for all ' num2str(n_tot_comps) ' components.']);
                else
                    o.cob_components_filter_inclusion = ones(n_tot_comps,1);
                    disp(['icatb_info [' char(datetime) '] cls_greedy_sort_components.m: componrent model order and bool vector not matching, having ' num2str(size(vec,1)) ' booleans for your model order of ' num2str(n_tot_comps) '.']);
                end
            else
                vec = [];
                o.cob_components_filter_inclusion = ones(n_tot_comps,1);
                disp(['icatb_info [' char(datetime) '] cls_greedy_sort_components.m: No booleans found for non noise. Using all components']);
            end

            n_err = o.m_calc_greed_match();
        end     

        function arAllComp = mpr_get_comps(~, ixMas, cob_ic_keep, sIcFile)
        
            v_tmp = icatb_spm_vol(sIcFile);
            n_tot_comps =  size(v_tmp,1);

            % A single 1 takes all the components within a template
            if isscalar(cob_ic_keep)
                if cob_ic_keep == 1
                    cob_ic_keep = ones(n_tot_comps,1);
                end
            else
                if ( size(cob_ic_keep,1) ~= n_tot_comps )
                    error(['icatb_err [' char(datetime) '] icatb_cls_greedy_sort_components.m: Filter vector of boleans carrying useful components vs noise does not match the number of components in component file/template'])
                end
            end

            %get ic loadings
            arAllComp=zeros(sum(cob_ic_keep),length(ixMas));
            iOrder=1;
            ix_loop = find(cob_ic_keep);
            for ixIc = ix_loop'
	            [VV, HInfo2] = icatb_returnHInfo([sIcFile ',' num2str(ixIc)]); 
	            icasig2 = icatb_spm_read_vols(VV);
	            structDIM2 = HInfo2.DIM;
	            icasig2 = reshape(icasig2, 1, prod(structDIM2));
	            icasig2(isnan(icasig2))=0; %clean nans
	            arAllComp(iOrder,:)=squeeze(icasig2(ixMas));
	            iOrder=iOrder+1;
	            clear icasig2;
            end
            arAllComp=arAllComp';
        end   


        function [ard_corrs_table, ari_ordered_pairs_table, ari_ordered_pairs_orig] = mpr_correlation(o, ixMas, s_file_icas, rob_filter, s_file_template)
            

            if isempty(ixMas)
                %create mask assuming a zero voxel is a outside of brain
                [VV1, HInfo1] = icatb_returnHInfo(s_file_icas); 
                ard_icas = icatb_spm_read_vols(VV1);
            
                [VV2, HInfo2] = icatb_returnHInfo(s_file_template); 
                ard_template = icatb_spm_read_vols(VV2);
            
                % zeros determines the voxels that will be masked out
                rod_mask = ard_icas .* ard_template;
            
                structDIM1 = HInfo1.DIM;
                rod_mask = reshape(rod_mask, 1, prod(structDIM1));
            
                ix_rod_mask=find(rod_mask ~= 0); 
                rod_mask(ix_rod_mask)=1; 
            
                ixMas=find(abs(rod_mask - 1) == 0);
                o.ron_mask_ix = ixMas;
            end

            mat_flat_vox_by_coms_compare = o.mpr_get_comps(ixMas, rob_filter, s_file_icas);
            mat_flat_vox_by_coms_template = o.mpr_get_comps(ixMas, 1, s_file_template); 

            if rob_filter == 1
                rob_filter = ones(size(mat_flat_vox_by_coms_compare,2),1);
                o.cob_components_filter_inclusion = rob_filter;
            end

            sim = corr(mat_flat_vox_by_coms_compare,mat_flat_vox_by_coms_template);
            nC = min(size(sim));
            for ii=1:nC
                [x,y] = find(sim == max(sim(:)));
                x = x(1); y = y(1);
                IDX1(ii,1) = x;
                IDX2(ii,1) = y;
                codCorr(ii,1) = sim(x,y);
                sim(x,:) = nan;
                sim(:,y) = nan;
            end
            ard_corrs_table = corr(mat_flat_vox_by_coms_compare,mat_flat_vox_by_coms_template);
            ari_ordered_pairs_table = [IDX1 IDX2];
            pos = find(rob_filter);
            for ii = 1:length(IDX1)
                IDX1orig(ii,1) = pos(IDX1(ii));
            end
            ari_ordered_pairs_orig = [IDX1orig IDX2];

        end

        function [ard_corrs_table, ari_ordered_pairs_table, ari_ordered_pairs_orig] = mpr_sort_greed(o, ixMas, s_file_icas, rob_filter, s_file_template)
            % sorts components between your components and a template

            disp(['icatb_info [' char(datetime) '] cls_greedy_sort_components.m: Starting Greedy Sort'])

            % Flatten mask
            [VV, HInfo2] = icatb_returnHInfo(o.s_mask); 
            icasig2 = icatb_spm_read_vols(VV);
            structDIM2 = HInfo2.DIM;
            icasig2 = reshape(icasig2, 1, prod(structDIM2));
            icasig2(isnan(icasig2))=0; %clean nans
            icasigMask=icasig2;
            clear icasig2;
            tol = 1e-2;
            ixMas=find(abs(icasigMask - 1) < tol);

            [o.ard_corrs_table, o.ari_ordered_pairs_table, o.ari_ordered_pairs_orig] = mpr_correlation(o, ixMas, o.s_components2compare, o.cob_components_filter_inclusion, o.s_template);

            save('o');

            disp(['icatb_info [' char(datetime) '] cls_greedy_sort_components.m: Completed Greedy Sort'])
            disp(['icatb_info [' char(datetime) '] cls_greedy_sort_components.m: Components sorted from file ' o.s_components2compare])
            disp(['icatb_info [' char(datetime) '] cls_greedy_sort_components.m: Template sorted against ' o.s_template])
            mkdir(fullfile(fileparts(o.s_components2compare),'utilities'));
            save(fullfile(fileparts(o.s_components2compare),['utilities' filesep 'greedsort' datestr(datetime('now'),'yyyymmddHHMMSS') '.mat']), 'o');
            disp(['icatb_info [' char(datetime) '] cls_greedy_sort_components.m: cls_greedy_sort_components.m: Greedy sort component pairs across files, saved in var o.ari_ordered_pairs_table under ' fullfile(fileparts(o.s_components2compare),['utilities' filesep 'greedsort' datestr(datetime('now'),'yyyymmddHHMMSS') '.mat'])])
            
        end

    end
end

