% -------------------------------------------------------
%
%    map  - class for storing geometric map data
%
%    Ver. 1.0.0
%
%    Created:           Mark Nothstein (25.02.2020)
%    Last modified:     Mark Nothstein (25.02.2020)
%
%    Institute of Biomedical Engineering
%    Karlsruhe Institute of Technology
%
%    http://www.ibt.kit.edu
%
%    Copyright 2000-2020 - All rights reserved.
%
% ------------------------------------------------------
%
% map(xyzin,facesin,namein,normalsin,varargin)
% class storing geometric map data.
%
% used to store map data
%
%
%
% Inputs:
% xyz       (points x xyz) 
% faces     (faces x xyz_id)
% normals   (normals x xyz_id)
% elecs     (legacy)
% ellocs    (legacy)
%
% Outputs:
%    map class
%
%
% Example Usage:
%   Has the same format as matlab "patch"
%   patch('Vertices',xyz,'Faces',faces,'EdgeColor','k','FaceColor','none','EdgeAlpha',0.5);
%
% Revision history:

classdef map < handle
    
    properties
        xyz
        faces
        normals
        elecs
        ellocs
        
        name
        status
        origin
        colormap
        comment
    end
    
    methods
        function obj = map(xyzin,facesin,namein,normalsin,varargin)
            %map(xyzin,facesin,namein,colormapin,normalsin,elecsin,ellocsin,statusin,originin,commentin)
            if nargin == 0
                obj.xyz = [];
                obj.faces  = [];
                obj.normals  = [];
                obj.elecs  = [];
                obj.ellocs  = [];
                obj.status  = [];
                obj.origin  = [];
                obj.comment  = [];
            else
                p = inputParser;
                addRequired(p,'xyzin',@(x) isnumeric(x) );
                addRequired(p,'facesin',@(x) isnumeric(x) );
                addRequired(p,'namein',@(x) ischar(x) && ~isempty(x) );
                addRequired(p,'normalsin',@(x) isnumeric(x) );
                addParameter(p,'colormapin',[],@(x) ~isempty(x) && isnumeric(x) );
                addParameter(p,'elecsin',[],@(x) ~isempty(x) && isnumeric(x) );
                addParameter(p,'statusin','',@(x) ischar(x) && ~isempty(x) );
                addParameter(p,'originin','',@(x) ischar(x) && ~isempty(x) );
                addParameter(p,'commentin','',@(x) ischar(x) && ~isempty(x) );
                addParameter(p,'ellocsin',[],@(x) isnumeric(x) && ~isempty(x) );
                try
                    p.parse( xyzin,facesin,namein,normalsin, varargin{:} );
                catch MyError
                    rethrow(MyError);
                end
                
                %needed values
                obj.name = p.Results.namein;
                obj.xyz = p.Results.xyzin;
                obj.faces = p.Results.facesin;
                obj.normals = p.Results.normalsin;
                
                %optional values
                obj.colormap = p.Results.colormapin;
                obj.elecs = p.Results.elecsin;
                obj.status = p.Results.statusin;
                obj.origin = p.Results.originin;
                obj.comment = p.Results.commentin;
                obj.ellocs = p.Results.normalsin;
                
                
                
            end
            
        end
        
        function a = plot(obj,id)
            if nargin<1; id=1; end
            name = 'surface mesh';
            
            a = figure;
            %patch('Vertices',xyz,'Faces',faces,'EdgeColor','k','FaceColor','none','EdgeAlpha',0.5);
            %patch('Vertices',xyz,'Faces',faces,'EdgeColor','k','EdgeAlpha',0.1,'FaceColor','k','FaceAlpha',0.0)
            if isempty(obj.colormap)
                patch('Vertices',obj.xyz,'Faces',obj.faces,'FaceColor',[0.4,0.4,0.4],'EdgeColor','none','FaceAlpha',0.1);
            else
                %NICE TO HAVE FOR VISUALISTION
                %mypatient.invasive.map.plot(3)
                %hold on
                %plot3(squeeze(mypatient.invasive.posterior.spi1.elec_loc_xyz_bi(1,1,:)),squeeze(mypatient.invasive.posterior.spi1.elec_loc_xyz_bi(1,2,:)),squeeze(mypatient.invasive.posterior.spi1.elec_loc_xyz_bi(1,3,:)),'k*','LineWidth',10)
                colr = obj.colormap(id).data;
                %colr(colr==0) = nan;
                %interpolate?
                TR = triangulation(double(obj.faces), double(obj.xyz));
                geoInd = getElecSurfaceVertex(TR.Points, obj.xyz);
                [~,interp_matrix] = generateLaplacianInterpolationMatrix(getElecSurfaceVertex,geoInd)
                
                patch('Vertices',obj.xyz,'Faces',obj.faces,'EdgeColor','k','FaceVertexCData',colr,'FaceColor','flat');
                colorbar
                
                 f1 = figure;
                ax1 = gca;
                data_struct_id = 1;
                data_struct.vertices = obj.xyz;
                data_struct.faces = obj.faces;
                data_struct.data = colr;
                plotGeometryWithData(ax1 ,data_struct, [], [quantile(data_struct(data_struct_id).data, .05), quantile(data_struct(data_struct_id).data, .95)],...
                    'hsv', 'interp');
    
                %Colormap
                pink = [255/255 100/255 175/255];
                num_cmap_vals = 500;
                fraction_colorbar_pink_vals = 0.5; % 50% pink
                max_voltage = 1;
                min_voltage = 0.1;
                max_lim = 4;
                min_lim = 0;
                num_colorscale_vals = num_cmap_vals*(1-fraction_colorbar_pink_vals);
                cmap = flipud(jet(num_colorscale_vals));
                cmap_size = size(cmap,1);
                for i=1:num_cmap_vals-num_colorscale_vals
                    cmap(cmap_size+i, :) = pink; %[0.85 0.4 0.8]; %some sort of pink
                end
                difference = max_lim-min_lim;
                if min_voltage==0 || min_lim==0
                    c_ax_low = 0;
                else
                    c_ax_low = min_voltage/min_lim;
                end
                c_ax_high = max_voltage/fraction_colorbar_pink_vals;
                c_ax_low = min_voltage/fraction_colorbar_pink_vals/(max_voltage/min_voltage);
                
                caxis([c_ax_low c_ax_high]);
                colorbar('Ticks',[1]);
                colorbar('Limits',[min_lim max_lim]);
                colormap(cmap);
                
                new_maxlim = fraction_colorbar_pink_vals*max_lim;
                
                cmap = flipud(jet(num_colorscale_vals));
                cmap(end+1,:)=pink;
                caxis([0.1 1.0]);
                
                axis equal
                axis off
                shading interp;
                alpha 0.7
            end
        end
        
        function import_colorgeo(obj,cdata)
            %adds colordata into colorgeo
            %input is struct with field "description" for name of data
            % and field data with as  many rows ans mesh elemnts 
            csiz = length(obj.colormap);
            cdatsiz = length(cdata);
            
            if isempty(obj.colormap)
                desc_lst = [];
            else
                desc_tbl = struct2table(obj.colormap);
                desc_lst = desc_tbl.description;
            end
            
            for inp = 1:cdatsiz
                if sum(strcmpi(desc_lst,cdata(inp).description))
                    fprintf('Description already exists. Skipping.\n')
                else
                    obj.colormap(csiz+1).description = cdata(inp).description;
                    obj.colormap(csiz+1).data = cdata(inp).data;
                    csiz = csiz+1;
                end
            end
            fprintf('Import of colordata finished.\n')
        end
        
        
    end
        
end

