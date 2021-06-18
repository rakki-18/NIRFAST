function data = new_femdata_spectral(mesh,frequency,wv)

    %
    % Calculates data (phase and amplitude) for a given
    % spectral me sh at a given frequency (MHz).
    % outputs phase and amplitude in structure data
    % and mesh information in mesh
    
    %
    % mesh is the input mesh (variable or filename)
    % frequency is the modulation frequency (MHz)
    % wv is optional wavelength array
    
    % This spectral mesh already has the mua, mus and kappa stored for each wavelength.
    parallel = parallel_init();
    
    %Number of nodes
    nodes = size(mesh.mua,1)/size(mesh.wv,1);
    % Store the mua, mus and kappa values for each wavelength for convenience
    MUA = []; MUS = []; KAPPA = [];
    for i = 1: size(mesh.wv,1)
        MUA = [MUA mesh.mua((i-1)*nodes+1:i*nodes)];
        MUS =  [MUS mesh.mus((i-1)*nodes+1:i*nodes)];
        KAPPA = [KAPPA mesh.kappa((i-1)*nodes+1:i*nodes)];
    end
    
    % error checking
    if frequency < 0
        errordlg('Frequency must be nonnegative','NIRFAST Error');
        error('Frequency must be nonnegative');
    end
    
    % If not a workspace variable, load mesh
    if ischar(mesh)== 1
      mesh = load_mesh(mesh);
    end
    
    if nargin == 2
      wv_array = sort(mesh.wv);
      linki = 1:length(wv_array);
    elseif nargin == 3
      wv_array = sort(wv);
      for i = 1:length(wv_array)
        flag = find(mesh.wv == wv_array(i));
        linki(i) = flag; % keeps correct link column with wv_array
        if isempty(flag) == 1
          disp('ERROR: wv_array contains wavelengths not present in the extinction coefficients')
          data = [];
          return
        end
      end
    end
    
    nwv = length(wv_array);
    
    %**************************************************************
    % Initiate log file
    % fid_log = fopen([mesh.name,'.log'],'w');
    % fprintf(fid_log,'**************************\n\nForward Model Log\n\n*************************\n\n');
    % fprintf(fid_log,'Forward Mesh  = %s\n',mesh.name);
    % fprintf(fid_log,'Wavelength Array  =  %dnm\n',wv_array);
    
    %****************************************************************
    % Run femdata for each wavelength and save data as .paa
    
    data.paa = []; data.wv=[]; data.phi = [];

      % disp(sprintf('Calculating data for: %g nm',(wv_array(1))))
      
      %****************************************************************
      % calculate absorption and scattering coefficients from concetrations and
      % scattering parameters a and b
      % This is already calculated in this case.
      % [mesh.mua, mesh.mus, mesh.kappa] = calc_mua_mus(mesh,wv_array(1));

      mesh = update(mesh, 1, MUA,MUS,KAPPA);
      % if sources are not fixed, move sources depending on mus
      if mesh.source.fixed == 0
        mus_eff = mesh.mus;
        [mesh]=move_source(mesh,mus_eff,3);
        clear mus_eff
      end
      
      mesh_temp = mesh;
      % build link file for first wv from first wv of function call.
      mesh_temp.link = [mesh.link(:,1:2) mesh.link(:,linki(1)+2)];
      % disp(sprintf('Calculating with mua start = %d and mus = %d and kappa = %d',mesh_temp.mua(1), mesh_temp.mus(1),mesh_temp.kappa(1)));
      % disp(size(mesh_temp.kappa));
      [data_single_wv,junk] = femdata_stnd(mesh_temp,frequency);
      data.paa = [data.paa, data_single_wv.paa];
      data.wv = [data.wv wv_array(1)];
      data.link = mesh_temp.link;
      % Convert from a sparse matrix to a full matrix
      data_single_wv.phi = full(data_single_wv.phi); 
      % Take average along rows as there are more than one source
      data_single_wv.phi = sum(data_single_wv.phi,2)/size(data_single_wv.phi,2); 
      data.phi = [data.phi; data_single_wv.phi];
      clear data_single_wv junk
      
      % PARALLEL
      if parallel
        parfor i = 2 : nwv
          data_single_wv=[];
          mus_eff=[];
          % disp(sprintf('Calculating data for: %g nm in parallel',(wv_array(i))))
    
          %****************************************************************
          % calculate absorption and scattering coefficients from concetrations and
          % scattering parameters a and b
          mesh_temp=mesh;
          mesh_temp.link = [mesh.link(:,1:2) mesh.link(:,linki(i)+2)];
          [mesh_temp.mua, mesh_temp.mus, mesh_temp.kappa] = calc_mua_mus(mesh_temp,wv_array(i));
    
          % if sources are not fixed, move sources depending on mus
          if mesh_temp.source.fixed == 0
            mus_eff = mesh_temp.mus;
            [mesh_temp]=move_source(mesh_temp,mus_eff,3);
           % clear mus_eff
          end
    
          [data_single_wv,mesh_temp] = femdata_stnd(mesh_temp,frequency);
          data_amp(:,i) = data_single_wv.paa(:,1);
          data_phase(:,i)= data_single_wv.paa(:,2);
          data_wv(i) = wv_array(i);
          data_link(:,i-1) = mesh.link(:,linki(i)+2);
          % clear data_single_wv
        end
      else
      % SERIAL
        for i = 2 : nwv
          data_single_wv=[];
          mus_eff=[];
          % disp(sprintf('Calculating data for: %g nm in serial',(wv_array(i))))
    
          %****************************************************************
          % calculate absorption and scattering coefficients from concetrations and
          % scattering parameters a and b
          mesh_temp=mesh;
          mesh_temp.link = [mesh.link(:,1:2) mesh.link(:,linki(i)+2)];
          
          % Already calculated in this case.
          % [mesh_temp.mua, mesh_temp.mus, mesh_temp.kappa] = calc_mua_mus(mesh_temp,wv_array(i));
          mesh_temp = update(mesh_temp, i, MUA,MUS,KAPPA);

          % if sources are not fixed, move sources depending on mus
          if mesh_temp.source.fixed == 0
            mus_eff = mesh_temp.mus;
            [mesh_temp]=move_source(mesh_temp,mus_eff,3);
           % clear mus_eff
          end
          
          % disp(sprintf('Calculating with mua start = %d and mus = %d and kappa = %d',mesh_temp.mua(1), mesh_temp.mus(1),mesh_temp.kappa(1)));
          % disp(size(mesh_temp.kappa));
          [data_single_wv,mesh_temp] = femdata_stnd(mesh_temp,frequency);
          data_amp(:,i) = data_single_wv.paa(:,1);
          data_phase(:,i)= data_single_wv.paa(:,2);
          data_wv(i) = wv_array(i);
          data_link(:,i-1) = mesh.link(:,linki(i)+2);
          % Convert from a sparse matrix to a full matrix
          data_single_wv.phi = full(data_single_wv.phi); 
          % Take average along rows as there are more than one source
          data_single_wv.phi = sum(data_single_wv.phi,2)/size(data_single_wv.phi,2); 
          data.phi = [data.phi; data_single_wv.phi];
          % clear data_single_wv
        end
      end
    
    for i=2:nwv
        data.paa=[data.paa,data_amp(:,i), data_phase(:,i)];
        data.wv(i)=data_wv(i);
        data.link(:,i+2) = data_link(:,i-1);
    end
        
    if isfield(mesh,'R')
        mesh=rmfield(mesh,'R');
    end
end

function mesh_temp = update(Mesh, i, MUA, MUS, KAPPA)
    mesh_temp  = Mesh;
    mesh_temp.mua = MUA(:,i);
    mesh_temp.mus = MUS(:,i);
    mesh_temp.kappa = KAPPA(:,i);
end 