classdef icatb_cls_greedy_sort_components < handle
    % Date 9/26/25
    % Cyrus Eierud
    % Code that greedily finds best correlated matches between your
    % components and a template
    % Example to calcuate and save state guided ICA:
    %ce092925 %       oc_sgica = icatb_cls_sgica(dfncInfo);
    %       oc_sgica = oc_sgica.set_s_outputDir(outputDir);
    %       oc_sgica = oc_sgica.set_input_FNCdynflat(FNCdynflat);
    %       n_ret = oc_sgica.calc_save;
    % Example to return the report data in structure variable:
    %       oc_sgica = icatb_cls_sgica(dfncInfo);
    %       oc_sgica = oc_sgica.set_s_outputDir(outputDir); 
    
    properties
        stru_sesInfo % parameter file
        s_mask % mask used to remove non brain voxels
        s_components2compare % your blind components
        i_components_order % model order number for your blind components
        s_template %template to compare against
        cob_components_filter_inclusion % components you want to keep as 1 and filter out as 0
        ard_corrs_table
        ari_ordered_pairs_table
        ari_ordered_pairs_orig
        %ce092925 coi_components_sorted %the blind source components you compare with template
        %ce092925 coi_template_sorted % the template you compared your components to
        %ce092925 cod_corrs_sorted % the correlation strength between your components and the template
    end

    
    methods
        function n_err = m_calc_greed_match(o)
            % sorts components between your components and a template

            %ce092625 create an error if mask, template or components have
            %incompatible resolutions 
            
            n_err = 1;

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

            mat_flat_vox_by_coms_compare = o.mpr_get_comps(ixMas, o.cob_components_filter_inclusion, o.s_components2compare);
            mat_flat_vox_by_coms_template = o.mpr_get_comps(ixMas, 1, o.s_template); 

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
            o.ard_corrs_table = corr(mat_flat_vox_by_coms_compare,mat_flat_vox_by_coms_template);
            o.ari_ordered_pairs_table = [IDX1 IDX2];
            pos = find(o.cob_components_filter_inclusion);
            for ii = 1:length(IDX1)
                IDX1orig(ii,1) = pos(IDX1(ii));
            end
            o.ari_ordered_pairs_orig = [IDX1orig IDX2];
            save('o');
            n_err = 0;
        end      

        function o = icatb_cls_greedy_sort_components(stru_sesInfo)
            % Constructor
            if nargin > 0
                o.stru_sesInfo = stru_sesInfo;
            end
            o.s_mask = fullfile(o.stru_sesInfo.outputDir, [o.stru_sesInfo.userInput.prefix, 'Mask.nii']); 
            o.s_components2compare = fullfile(o.stru_sesInfo.outputDir, o.stru_sesInfo.icaOutputFiles(1).ses.name); 
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

        function h_popup = mpr_dialog(o)
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
            
            textPos(2) = titlePos(2) - 0.15;
            
%             text('units', 'normalized', 'string', 'ce093025b', 'position', textPos, 'fontsize', textFont, 'HorizontalAlignment', 'left', ...
%                 'fontweight', 'normal', 'FontName', UI_FONTNAME);

            textPos(2) = textPos(2) - 0.15;

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

            disp(['icatb_info [' char(datetime) '] cls_greedy_sort_components.m: Starting Greedy Sort'])
            n_err = o.m_calc_greed_match();
            disp(['icatb_info [' char(datetime) '] cls_greedy_sort_components.m: Completed Greedy Sort'])
            disp(['icatb_info [' char(datetime) '] cls_greedy_sort_components.m: Components sorted from file ' o.s_components2compare])
            disp(['icatb_info [' char(datetime) '] cls_greedy_sort_components.m: Template sorted against ' o.s_template])
            mkdir(fullfile(fileparts(o.s_components2compare),'utilities'));
            save(fullfile(fileparts(o.s_components2compare),['utilities' filesep 'greedsort' datestr(datetime('now'),'yyyymmddHHMMSS') '.mat']), 'o');
            disp(['icatb_info [' char(datetime) '] cls_greedy_sort_components.m: Greedy sort data saved under ' fullfile(fileparts(o.s_components2compare),['utilities' filesep 'greedsort' datestr(datetime('now'),'yyyymmddHHMMSS') '.mat'])])

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
   
    end
end

