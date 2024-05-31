%%% pTX3DPulseDesigner.m --- a class that designs 3D spoke pTX pulses. The
%spokes are placed symmetrically on a circle in kspace. the placement is
%optimized over different radiuses and rotation angles.
%
% 
%%
%% 
%% Filename: pTX3DPulseDesigner.m
%% Description: 
%% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
%% Maintainer: 
%% Copyright (C) 2010 CMRR at UMN
%% Created: Thu Aug  5 20:03:38 2010 (CDT)
%% Version: 
%% Last-Updated: Tue Mar 12 22:31:43 2013 (CDT)
%%           By: Xiaoping Wu
%%     Update #: 239
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%% Commentary: 
%% 
%% 
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%% Change log:
%% 
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%% Code:

classdef pTX3DPulseDesigner < handle
    properties 
      B1Map = []
      Mask = []
      Target = []
      FieldOfExcit = [] % cm
      NumOfSpokes = 3
      SubRFType = 'gauss'
      Thickness = 5 % mm
      SubRFDuration = 0.5 % ms
      NominalFOX = 10:2:36 % cm
      ConvergenceTolerance = 1e-5
      Lambda = [1e-8 1e-7 5e-7 1e-6 2e-6 5e-6 8e-6 1e-5 2e-5 5e-5]
      OptimalLambda = 5e-7
      NominalFlipAngle = 10 % deg
      B0Map = []
      MaxGradSlewRate = 160 % T/m/s
      MaxGradAmplitude = 50 % mT/m
      DwellTime = 4 % us
      ReadOutOffset = [0 0 0] % mm
      Kappa = 1 % a scaling factor used for VERSE
      OptimalFOX = 10
    end % properties
    
    properties (SetAccess = protected)
      NominalTheta = 0:10:170 % deg      
    end % properties

    properties (SetAccess = private)
      IsOptimal = false      
      OptimalTheta = 0
      ResidualError = 0
    end % properties

 
    methods 
        function obj = pTX3DPulseDesigner (b1map, mask, fox)
        % Constructor
        % Usage: obj = classname (b1map, mask, fox)
        % b1map:
        % mask:
        % fox: in cm
        
        % Pre init.
        % anything not using output (obj)

          if nargin == 0
            disp(['-> The object is constructed. Please specify necessary ' ...
                  'properties...'])
            return
          end

        % compvalue = classname.staticMethod();

        
        % Object init.
        % call super class, if applicable, before accessing object
        % obj = obj@superClass(args{:});

        
        % Post init.
        % anything including accessing object
        % obj.classMethod();
        % obj.Property = compvalue;
          obj.B1Map = b1map;
          obj.Mask = mask;
          obj.FieldOfExcit = fox;
          
          obj.Target = mask;
          
          disp('-> Object constructed. Please specify other needed properties...')
        end

        function obj = set.Kappa (obj, kappa)
           if kappa<=0 || kappa>1
             kappa = 1;
             disp('-> Kappa has to be 0< kappa<= 1)'); 
           end
           
           obj.Kappa = kappa;
        end
        
        function obj = set.NumOfSpokes (obj, nspokes)
        % Purpose:
        % Usage: output = obj.funcname (a,b,c)
        % a:
        % b:
        % c:

          if nspokes < 1
            nspokes = 1;
            disp('-> NumOfSpokes must be >= 1 and is set to 1.')
          end
                         
          if ~isequal(obj.NumOfSpokes,nspokes)          
            obj.NumOfSpokes = nspokes;          
            obj.IsOptimal = false;
          end        
        
        end
        
        function obj = set.Target (obj, targ)
          if ~isequal(obj.Target, targ) 
            obj.Target = targ;
            obj.IsOptimal = false;
          end
        end

        function obj = set.NominalFOX (obj, myfox)
        % Purpose:
        % Usage: output = obj.funcname (a,b,c)
        % a:
        % b:
        % c:

        if ~isequal(obj.NominalFOX, myfox)
          obj.NominalFOX = myfox;        
          obj.IsOptimal = false;
        end
        
        end
        
        function obj = set.SubRFType (obj, rftype)
        % Purpose:
        % Usage: output = obj.funcname (a,b,c)
        % a:
        % b:
        % c:
        if ~(strcmpi(rftype,'gauss')||strcmpi(rftype,'sinc'))
          disp('SubRFType must be gauss or sinc')
          rftype = 'gauss';
        end
        
        obj.SubRFType = rftype;

        end
        
        function obj = set.B0Map (obj, b0map)
        % Purpose:
        % Usage: output = obj.funcname (a,b,c)
        % a:
        % b:
        % c:
        if isempty(b0map)
          b0map = zeros(size(obj.Mask));
        end  
        
        if ~isequal(obj.B0Map, b0map)
          obj.B0Map = b0map;
          obj.IsOptimal = false;
        end

        end


        function stat = plotLCurve (obj)
        % Purpose:
        % Usage: stat = obj.plotLCurve ()
        % 
        % 
        % 
        lambda = obj.Lambda;
        myfox = obj.NominalFOX;
        mytheta = obj.NominalTheta;
        m = obj.constructTargetVector;
        
        [~,rho,eta] = obj.calcWeights(mean(myfox),mytheta(1),lambda,m);
        figure, plot_lc(rho,eta,'o',1,lambda)
        
        stat = 1;
        end

        function [rfs,grad,kp,m,sysmat,wt1,mxypat,rfi,gm] = design (obj)
        % Purpose:
        % Usage: [rfs,grad,kp,myrho] = obj.design ()
        % rfs:
        % grad:
        % kp:
          m = obj.constructTargetVector;
          fox0 = obj.OptimalFOX;
          %fox0= obj.NominalFOX;
          theta0 = obj.OptimalTheta;
          %lambda0 = obj.OptimalLambda;
          myrho = obj.ResidualError;
          
          if ~obj.IsOptimal && (obj.NumOfSpokes>1)
            [fox0,theta0,myrho,mykp,mykp1] = obj.findOptimals(m);end

          [wt1,~,~,rfsub1,grad,kp,sysmat] = obj.calcWeights(fox0,theta0,m);

          dt = 1e-6*obj.DwellTime;
          if obj.Kappa< 1
            zz = find(grad(3,:)==0);
            gspoke = grad(3,1:zz(2));
            gamp = max(abs(gspoke));
            b1max = max(abs(rfsub1));
            
            opt.KEEPDUR=1;
            opt.KAPPA=obj.Kappa;
            opt.GRT = dt;
            opt.VGRT = dt;
            opt.RFOSFACTOR = 1;
            [rfsubv,gspokev] = verseItAdv(rfsub1(gspoke==gamp).',gamp,opt);
            rfsubv = b1max.*rfsubv(:).';
            gspokev = gspokev(:).';
            
            srmax = obj.MaxGradSlewRate;
            gmax = obj.MaxGradAmplitude;
            
            [gspoke1,np1] = add_ramp(gspokev,0,srmax,dt); % rampup
            [gspoke1,np2] = add_ramp(fliplr(gspoke1),0,srmax,dt); % rampdown
            gspoke1 = fliplr(gspoke1);
            grad = assemble_grad_spoke3d(gspoke1,kp,dt,gmax,srmax);
            
            rfsub1 = [zeros(1,np1) rfsubv zeros(1,np2)];
          end
          
          b1map = obj.B1Map;
          b0map = obj.B0Map;
          mask = obj.Mask;
          fox = obj.FieldOfExcit;
          %fa0 = obj.NominalFlipAngle;

          nchs = size(b1map,4);
          wt = reshape(wt1,[],nchs);
          rfs = assemble_rf_spoke3d(rfsub1,wt);

          ns = size(mask,3);
          soi = round(ns/2);
          roOffset= obj.ReadOutOffset;
          mxypat = run_bloch_sim (rfs,grad(1:2,:),...
                                  b1map(:,:,soi,:),...
                                  mask(:,:,soi),...
                                  1e-2*fox, ...
                                  b0map(:,:,soi),0,[],dt,roOffset);

           
           m= abs(m).*exp(1i*angle(sysmat*wt1));
           
       gamma = 2.675e8;    
       kp1=[kp [0;0]]; % ends at 0
       gm= diff(kp1/gamma,1,2);

       rfi= 1e3*dt*wt.';
           
        end

    end % methods

    methods (Access = protected)
      
      function subrf = designSubRF (obj)
    % Purpose:
    % Usage: output = obj.funcname (a,b,c)
    % a:
    % b:
    % c:
      rftype = obj.SubRFType;
      tp = obj.SubRFDuration;
      dt = 1e-6* obj.DwellTime; % s
      switch rftype
       case 'gauss'
        subrf = design_rf_gauss(tp, 10, 0.01, dt);
       case 'sinc'
        subrf = design_rf_sinc(tp, 10,1, dt);
       otherwise
        subrf = [];
      end
      
      subrf= subrf./sum(subrf);
      end
      
      function [fox0,theta0,myrho,mykp,mykp1] = findOptimals (obj, m)
      % Purpose:
      % Usage: [fox0,theta0,myrho] = findOptimals (obj)
      % 
      % 
      % 
      myfox = obj.NominalFOX;
      mytheta = obj.NominalTheta;
      %lambda0 = obj.OptimalLambda;
      
      myrho = zeros(length(myfox),length(mytheta));
      mykp = cell(length(myfox),length(mytheta));
      mykp1 = mykp;
      for ifox= 1:length(myfox),  
        myfox1 = myfox(ifox);
        %parfor jtheta= 1:length(mytheta)
        for jtheta= 1:length(mytheta)
          [~,rho,~,~,~,kp,~,kp1] = obj.calcWeights(myfox1,mytheta(jtheta),m);
          myrho(ifox,jtheta) = rho;
          mykp{ifox,jtheta} = kp;
          mykp1{ifox,jtheta} = kp1;
        end
      end

      [ifox,jtheta] = find(myrho==min(myrho(:)));
      fox0= myfox(ifox(1));
      theta0= mytheta(jtheta(1));
      
      obj.OptimalFOX = fox0;
      obj.OptimalTheta = theta0;
      obj.IsOptimal = true;
      obj.ResidualError = myrho;
      
      end
      
      function [wt,rho,eta,rfsub1,grad,kp,sysmat,kp1] = calcWeights (obj,nfox,theta,m)
      % Purpose:
      % Usage: output = obj.funcname (a,b,c)
      % a:
      % b:
      % c:
        
        nspokes = obj.NumOfSpokes;        
        %     
        dk=2*pi/(0.01*nfox);
        theta= deg2rad(theta);
        switch nspokes
            case 1
                kp=[0;0];
            case 2
                kp=[0 0;dk*sin(theta) dk*cos(theta)]';
            case 3
                kp=[0 0;dk*sin(theta) dk*cos(theta); -dk*sin(theta) -dk*cos(theta)]';
            case 4
                kp=[0 0; dk*sin(theta) dk*cos(theta); -dk*sin(theta) -dk*cos(theta);...
                    dk*sin(theta+pi/2) dk*cos(theta+pi/2)]';
            case 5
                kp=[0 0; dk*sin(theta) dk*cos(theta); -dk*sin(theta) -dk*cos(theta);...
                    dk*sin(theta+pi/2) dk*cos(theta+pi/2); -dk*sin(theta+pi/2) -dk*cos(theta+pi/2)]';
            otherwise
        end
        
        %kp = place_spoke_symmetric(nfox,nspokes,theta);
        
        rfsub = obj.designSubRF;
        thk = 1e-3* obj.Thickness; % m
        maxsr = obj.MaxGradSlewRate;
        maxamp = obj.MaxGradAmplitude;
        dt = 1e-6* obj.DwellTime; % s
        
        [grad,rfsub1] = design_spoke3d (rfsub,thk,kp,false,maxsr,maxamp,dt);
                
        b1map = obj.B1Map;
        mask = obj.Mask;
        fox = 1e-2* obj.FieldOfExcit; % m
        b0map = obj.B0Map;

        sysmat = complex(zeros(length(find(mask(:))),nspokes*size(b1map,4)));
        nslices = size(mask,3);
        idx0= 0;
        for islice = 1:nslices,
          imask = mask(:,:,islice);
          idxn = idx0 + length(find(imask(:)));
          [sysmat(idx0+1:idxn,:),kp1] = construct_sysmat_spoke3d (grad, ...
                                                            b1map(:,:, islice,:), ...
                                                            mask(:,:,islice), ...
                                                            fox, b0map(:,:,islice), ...
                                                            dt,obj.ReadOutOffset); ...
              

          idx0= idxn;
        end
        
        tol = obj.ConvergenceTolerance;
        [U,s,~]=csvd(sysmat);
        lambda= l_curve(U,s,m);
          
        if ~obj.IsOptimal %|| length(lambda)>1|| nspokes>1 % called by plotLCurve or multi spoke design
          [wt,rho,eta]= solve_mlstr(sysmat,m,lambda,tol);
        else  % called by single spoke design with a single specified lambda
          n= size(sysmat,2);
          ntrials= 600;
%           if nspokes>1
%               ntrials=24*4;
%           end
          wts= complex(zeros(n,ntrials));
          x0s= rand(n,ntrials)+ 1i*rand(n,ntrials);
          maxmins= zeros(1,ntrials);
          disp('-> Diff trials of initial points started...')
          parfor itrial=1:ntrials
            [iwt,irho,~]= solve_mlstr(sysmat,m,lambda,tol,x0s(:,itrial));
            im= sysmat* iwt;
            maxmins(itrial)= (max(abs(im))-min(abs(im)))./ mean(abs(im));
            wts(:,itrial)= iwt;
            rhos(itrial)= irho;
          end
        
          myidx= find(maxmins== min(maxmins));
          wt= wts(:,myidx(1));
          rho=rhos(:,myidx(1));
          eta=[];
        end
        
      end
      
      function m = constructTargetVector (obj)
      % Purpose:
      % Usage: output = obj.funcname (a,b,c)
      % a:
      % b:
      % c:
      mask = obj.Mask;
      targ = obj.Target;
      m = zeros(length(find(mask)),1);
      fa0= obj.NominalFlipAngle;
      
      phs= angle(sum(obj.B1Map,4));
      
      nslices = size(mask,3);
      idx0= 0;
      for islice = 1:nslices,
        imask = mask(:,:,islice);
        itarg = targ(:,:,islice);
        iphs= phs(:,:,islice);
        idxn = idx0 + length(find(imask(:)));
        m(idx0+1:idxn) = exp(1i*iphs(imask)).*construct_target_vector(itarg,imask,fa0,'l');        
        idx0= idxn;
      end

      end


    end % methods

end % classdef




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% pTX3DPulseDesigner.m ends here
