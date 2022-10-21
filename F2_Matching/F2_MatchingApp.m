classdef F2_MatchingApp < handle & F2_common
  properties
    guihan
    QuadScanData
    TwissFitSource string {mustBeMember(TwissFitSource,["Model","Analytic"])} = "Model"
    ProfFitMethod string {mustBeMember(ProfFitMethod,["Gaussian","Asymmetric"])} = "Asymmetric"
    LiveModel
    Optimizer string {mustBeMember(Optimizer,["fminsearch","lsqnonlin"])} = "lsqnonlin"
    DimSelect string {mustBeMember(DimSelect,["X" "Y" "XY"])} = "XY"
    LM
    ShowPlotLegend logical = false
    UseMatchQuad logical % Which matching quads to use
    LoadingScanData logical = false % Sets to true when loading raw emittance scan data (overrides loading in externally calculated Twiss parameters)
    MatchDir string {mustBeMember(MatchDir,["BKW" "FWD"])} = "BKW"
  end
  properties(SetAccess=private)
    MatchQuadNames string = ["QUAD:IN10:425" "QUAD:IN10:441" "QUAD:IN10:511" "QUAD:IN10:525"] % Read in when select Prof device and/or when change # matching quads
    TwissFitModel(1,8) = [1 1 0 0 1 1 5 5] % beta_x,beta_y,alpha_x,alpha_y,bmag_x,bmag_y,nemit_x,nemit_y fitted twiss parameters at profile device (nemit in um-rad)
    TwissFitAnalytic(1,8) = [1 1 0 0 1 1 5 5] % beta_x,beta_y,alpha_x,alpha_y,bmag_x,bmag_y,nemit_x,nemit_y fitted twiss parameters at profile device (nemit in um-rad)
    MultiWireData(2,4) % [x|y , wire # 1-4] [um]
    MultiWireDataErr(2,4) % [x|y , wire # 1-4] [um]
  end
  properties(SetAccess=private,SetObservable)
    goodmatch logical = false
  end
  properties(SetAccess=private,Hidden)
    InitMatch
    TwissMatch
    TwissPreMatch
    MatchQuadInitVals
    MatchQuadID % ID into LiveModel.Mags
    ModelQuadScan % [2 x Nscan]
    InitRestore
    ProfModelInd
    MatchQuadModelInd
    LEMQuadID
    ScanDataSelect
  end
  properties(SetObservable,AbortSet)
    ModelSource string {mustBeMember(ModelSource,["Live" "Archive" "Design"])} = "Live"
    ModelDate(1,6) = [2021,7,1,12,1,1] % [yr,mnth,day,hr,min,sec]
    UseFudge logical = false
  end
  properties(SetObservable)
    ProfName string = "PROF:IN10:571"
    NumMatchQuads uint8 {mustBeGreaterThan(NumMatchQuads,3)} = 4 % Number of quadrupoles upstream of profile monitor to use in match
    UndoAvailable logical = false
  end
  properties(Dependent)
    quadscan_k
    TwissFit(1,8)
  end
  properties(Constant)
    plotcols = 'bgrm'
    EmitDataProfs string = [] % profile devices for which there are NO emittance PVs
%     EmitDataProfs = ["PROF:IN10:571" "PROF:LI11:375" "CAMR:LI20:103" "WIRE:IN10:561" "WIRE:LI11:444" "WIRE:LI19:144"] % profile devices for which there are NO emittance PVs
    InitMatchProf = ["WIRE:IN10:561","PROF:IN10:571"] % Profile devices to associate with initial match conditions
    MWNames_L2 = ["WIRE:LI11:444","WIRE:LI11:614","WIRE:LI11:744","WIRE:LI12:214"]
    MWNames_L3 = ["WIRE:LI18:944","WIRE:LI19:144","WIRE:LI19:244","WIRE:LI19:344"]
  end
 
  methods
    function obj = F2_MatchingApp(ghan)
      if exist('ghan','var') && ~isempty(ghan)
        obj.guihan=ghan;
      end
      obj.LiveModel = F2_LiveModelApp ; % Initially update live model (default)
      obj.LM=copy(obj.LiveModel.LM);
      obj.LM.ModelClasses=["PROF" "WIRE"];
      obj.ProfModelInd = obj.LM.ModelID(obj.LM.ControlNames==obj.ProfName) ;
      obj.ProfModelInd = obj.ProfModelInd(1) ;
      obj.UseArchive = true ; % default to getting data from archive from now on
      obj.LiveModel.ModelSource = "Live" ;
      obj.LEMQuadID = obj.LiveModel.LEM.Mags.LM.ModelClassList == "QUAD" ; % Tag quadrupoles in LEM magnet list
%       obj.LiveModel.LEM.Mags.WriteEnable=false; % disable writing to magnets to test
    end
    function DoMatch(obj,cmd)
      %DOMATCH Perform matching to profile monitor device based on fitted Twiss parameters and live model
      %DoMatch() Match initial Twiss parameters and quad settings to restore Model
      %DoMatch("init") Match only initial Twiss parameters
      global BEAMLINE PS
      obj.goodmatch = false ;
      obj.MatchQuadInitVals = []; obj.TwissPreMatch=[]; obj.TwissMatch=[];
      LM_mags=obj.LiveModel.LEM.Mags.LM;
      pele=obj.ProfModelInd;
      betax_fit = obj.TwissFit(1) ;
      betay_fit = obj.TwissFit(2) ;
      alphax_fit = obj.TwissFit(3) ;
      alphay_fit = obj.TwissFit(4) ;
      
      % Form list of matching quads
      if obj.MatchDir=="BKW"
        iquads = LM_mags.ModelID(:) < pele & obj.LEMQuadID(:) ;
      else
        iquads = LM_mags.ModelID(:) > pele & obj.LEMQuadID(:) ;
      end
      if sum(iquads)<double(obj.NumMatchQuads)
        error('Insufficient matching quads available');
      end
      if obj.MatchDir=="BKW"
        id1=find(iquads,double(obj.NumMatchQuads),'last') ; id1=id1(obj.UseMatchQuad);
      else
        id1=find(iquads,double(obj.NumMatchQuads),'first') ; id1=id1(obj.UseMatchQuad);
      end
      quadps = arrayfun(@(x) BEAMLINE{x}.PS,LM_mags.ModelID(id1)) ;
      idm=ismember(obj.LiveModel.LEM.Mags.LM.ModelID,LM_mags.ModelID(id1)) ;
      bmin = obj.LiveModel.LEM.Mags.BMIN(idm)./10;
      bmax = obj.LiveModel.LEM.Mags.BMAX(idm)./10;
      
      % Get initial PS settings
      ps_init = arrayfun(@(x) PS(x).Ampl,quadps) ;
      obj.MatchQuadInitVals = ps_init.*10 ;
      
      % If matching on PR10571 or WS10561 then get initial Twiss parameters and record them in PVs
%       if ismember(obj.ProfName,obj.InitMatchProf)
%         i1=1;
%       else % match from entrance of first matching quad or match point
        if obj.MatchDir=="BKW"
          i1=PS(quadps(1)).Element(1);
          i2=pele;
        else
          i1=pele;
          i2=PS(quadps(end)).Element(end);
        end
%       end
      
      % Design twiss parameters at destination location
      betax_design = obj.LiveModel.DesignTwiss.betax(i2) ;
      betay_design = obj.LiveModel.DesignTwiss.betay(i2) ;
      alphax_design = obj.LiveModel.DesignTwiss.alphax(i2) ;
      alphay_design = obj.LiveModel.DesignTwiss.alphay(i2) ;

      % Get Twiss parameters at initial match point by back-propogating using flipped Lattice
      I = TwissToInitial(obj.LiveModel.DesignTwiss,i1,obj.LiveModel.Initial);
      if obj.MatchDir=="BKW"
        isol=findcells(BEAMLINE,'Class','SOLENOID');
        for ind=isol; BEAMLINE{ind}.B=0; end
        [~,R]=RmatAtoB(i1,i2);
        S1=[betax_fit   -alphax_fit                0          0 ;
            -alphax_fit (1+alphax_fit^2)/betax_fit 0          0 ;
            0           0                          betay_fit   -alphay_fit ;
            0           0                          -alphay_fit (1+alphay_fit^2)/betay_fit ] ;
        R=[R(1,1) R(1,2) 0      0      ;
           R(2,1) R(2,2) 0      0      ;
           0      0      R(3,3) R(3,4) ;
           0      0      R(4,3) R(4,4) ] ;
        S0=inv(R)*S1*inv(R'); %#ok<MINV>
        if any(isinf(S0(:))) || any(isnan(S0(:)))
          error('Match of Initial Twiss Parameters Failed');
        end
        I.x.Twiss.beta=S0(1,1); I.x.Twiss.alpha=-S0(1,2); I.y.Twiss.beta=S0(3,3); I.y.Twiss.alpha=-S0(3,4);
      else
        I.x.Twiss.beta=betax_fit; I.x.Twiss.alpha=alphax_fit; I.y.Twiss.beta=betay_fit; I.y.Twiss.alpha=alphay_fit;
      end
      I.Momentum = BEAMLINE{i1}.P;
      
      % - Fine tune initial Twiss parameters using Matching
      if obj.MatchDir=="BKW"
        M=Match;
        M.beam=MakeBeam6DGauss(I,1e3,3,1);
        M.iInitial=i1;
        M.initStruc=I;
        M.verbose=false; % see optimizer output or not
        M.optim='fminsearch';
        M.optimDisplay='off';
        M.useParallel=false;
        M.addVariable('INITIAL',1,'betax',0.01,1000);
        M.addVariable('INITIAL',1,'alphax',-100,100);
        M.addVariable('INITIAL',1,'betay',0.01,1000);
        M.addVariable('INITIAL',1,'alphay',-100,100);
        M.addMatch(pele,'beta_x',betax_fit,betax_fit/1000);
        M.addMatch(pele,'alpha_x',alphax_fit,1e-3);
        M.addMatch(pele,'beta_y',betay_fit,betay_fit/1000);
        M.addMatch(pele,'alpha_y',alphay_fit,1e-3);
        M.doMatch;
        display(M);
        if any(abs(M.optimVals-[betax_fit alphax_fit betay_fit alphay_fit])>[betax_fit/100 0.01 betay_fit/100 0.01])
          error('Match of Initial Twiss Parameters Failed');
        end
        obj.InitMatch=M.initStruc;
      else
        obj.InitMatch=I;
      end
      obj.InitMatch.x.NEmit = obj.TwissFit(7)*1e-6 ;
      obj.InitMatch.y.NEmit = obj.TwissFit(8)*1e-6 ;
      [~,T]=GetTwiss(i1,i2,obj.InitMatch.x.Twiss,obj.InitMatch.y.Twiss);
      obj.TwissPreMatch=T;
      obj.TwissPreMatch.z = [arrayfun(@(x) BEAMLINE{x}.Coordi(3),i1:i2) BEAMLINE{i2}.Coordf(3)] ;
      % If this is all is required, exit now
      if exist('cmd','var') && cmd=="init"
        return
      end
      
      % Calculate re-match to fitted Twiss parameters at profile monitor device
      M=Match;
      M.beam=MakeBeam6DGauss(obj.InitMatch,1e3,3,1);
      M.iInitial=i1;
      M.initStruc=obj.InitMatch;
      M.verbose=false; % see optimizer output or not
      M.optim=char(obj.Optimizer);
      M.optimDisplay='iter';
      M.useParallel=false;
      for ips=1:length(quadps)
        M.addVariable('PS', quadps(ips),'Ampl',bmin(ips),bmax(ips));
      end
      M.addMatch(i2,'beta_x',betax_design,betax_design/1e3);
      M.addMatch(i2,'beta_y',betay_design,betay_design/1e3);
      M.addMatch(i2,'alpha_x',alphax_design,1e-3);
      M.addMatch(i2,'alpha_y',alphay_design,1e-3);
      M.doMatch;
      display(M);
      [~,T]=GetTwiss(i1,i2,obj.InitMatch.x.Twiss,obj.InitMatch.y.Twiss);
      obj.TwissMatch=T;
      obj.TwissMatch.z = [arrayfun(@(x) BEAMLINE{x}.Coordi(3),i1:i2) BEAMLINE{i2}.Coordf(3)] ;
      
      % Check match in range and store in Mags BDES field
      if any(M.varVals(:)>bmax(:) | M.varVals(:)<bmin(:))
        for ips=1:length(quadps)
          PS(quadps(ips)).Ampl = ps_init(ips) ;
        end
        error('Required quad match values outside limits');
      end
      if any(abs(M.optimVals-[betax_design betay_design alphax_design alphay_design])>[betax_design/100 betay_design/100 0.01 0.01])
        for ips=1:length(quadps)
          PS(quadps(ips)).Ampl = ps_init(ips) ;
        end
        error('Match failed to converge to required accuracy');
      end
      
      obj.MatchQuadID = idm ;
      obj.LiveModel.LEM.Mags.BDES(idm) = M.varVals * 10 ;
      obj.LiveModel.LEM.Mags.SetBDES_err(false,~idm) ; % set just matching quads to be written
      obj.QuadScanData.twissmatch = M.matchVals ;
      obj.QuadScanData.quadmatch = M.varVals * 10 ;
      
      % Restore pre-matched magnet strengths in model
      for ips=1:length(quadps)
        PS(quadps(ips)).Ampl = ps_init(ips) ;
      end
      
      obj.goodmatch = true ;
      obj.UndoAvailable = false ;
      
    end
    function msg = RestoreMatchingQuads(obj)
      %RESTOREMATCHINGQUADS Undo last operation of WriteMatchQuads
      
      if isempty(obj.MatchQuadInitVals) || isempty(obj.MatchQuadID)
        error('No stored match quad settings to restore');
      end
      obj.LiveModel.LEM.Mags.BDES(obj.MatchQuadID) = obj.MatchQuadInitVals ;
      obj.LiveModel.LEM.Mags.SetBDES_err(false,~obj.MatchQuadID) ;
      obj.LiveModel.LEM.Mags.SetBDES_err(true,obj.MatchQuadID) ;
      msg = obj.LiveModel.LEM.Mags.WriteBDES ;
      if ~isempty(obj.InitRestore)
        obj.LiveModel.Initial = obj.InitRestore;
      end
      obj.goodmatch = false ;
      obj.UndoAvailable = false ;
      
      if obj.ModelSource ~= "Design"
        obj.LiveModel.UpdateModel ;
      end
      
    end
    function msg = WriteMatchingQuads(obj)
      %WRITEMATCHQUADS Write matched quadrupole fields to control system and update any Twiss PVs
      
      % Require DoMatch to have successfully run
      if ~obj.goodmatch
        error('No good match solution found, aborting quad writing');
      end
      
      % Write matched values to control system
      msg = obj.LiveModel.LEM.Mags.WriteBDES ;
      
      obj.UndoAvailable = true ;
      
      if obj.ModelSource ~= "Design"
        obj.LiveModel.UpdateModel ;
      end
      
    end
    function didload=LoadQuadScanData(obj)
      %LOADSCANDATA Load corr plot data for quad scan or multiwire
      %didload=1: quad scan didload=2: multiwire else error
      
      didload=0;
%       
%       obj.TwissMatch=[];
%       obj.TwissPreMatch=[];
%       obj.InitMatch=[];
      obj.goodmatch = false ;
%       obj.QuadScanData = [] ;
      LM_mags=obj.LiveModel.LEM.Mags.LM;
      LM_all=obj.LiveModel.LM;
      
      % If there is a QUAD scan file in the current day's directory, default to loading that
      files = dir(obj.datadir);
      qscanfiles = startsWith(string(arrayfun(@(x) x.name,files,'UniformOutput',false)),"CorrelationPlot-QUAD") | ...
       cellfun(@(xx) ~isempty(xx),regexp(string(arrayfun(@(x) x.name,files,'UniformOutput',false)),"CorrelationPlot-LI%d%d:QUAD", 'once')) | ...
       startsWith(string(arrayfun(@(x) x.name,files,'UniformOutput',false)),"Emittance-scan") | ...
       startsWith(string(arrayfun(@(x) x.name,files,'UniformOutput',false)),"Emittance-multi-WIRE");
      if any(qscanfiles)
        [~,latestfile]=max(datenum({files(qscanfiles).date}));
        ifile=find(qscanfiles); ifile=ifile(latestfile);
        [fn,pn]=uigetfile('*.mat','Select Cor Plot File',fullfile(obj.datadir,files(ifile).name));
      else % just offer up current data dir
        [fn,pn]=uigetfile('*.mat','Select Cor Plot File',obj.datadir);
      end
      if ~ischar(fn)
        return
      end
      % Update live model to match data file date
      dd = dir(fullfile(pn,fn));
      obj.ModelDate = datevec(dd.datenum) ;
      dat = load(fullfile(pn,fn));
      
      % If multi-wire, just need to upload data into multi-wire GUI boxes and done
      if startsWith(string(fn),"Emittance-multi-WIRE")
        obj.ProfName = string(dat.data.name{1}) ;
        if ~isempty(obj.guihan)
          obj.guihan.TabGroup.SelectedTab = obj.guihan.MultiWireEmittanceTab;
        end
        if isfield(dat.data.beam,'xStat')
          dim=1;
          obj.MultiWireData(1,:) = arrayfun(@(x) dat.data.beam(x,2).xStat(3),1:4) ;
          obj.MultiWireDataErr(1,:) = arrayfun(@(x) dat.data.beam(x,2).xStatStd(3),1:4) ;
        else
          dim=2;
          obj.MultiWireData(2,:) = arrayfun(@(x) dat.data.beam(x,2).yStat(3),1:4) ;
          obj.MultiWireDataErr(2,:) = arrayfun(@(x) dat.data.beam(x,2).yStatStd(3),1:4) ;
        end
        if ~isempty(obj.guihan)
          if dim==1
            obj.guihan.MeasurementPlaneDropDown.Value = "Horizontal" ;
          else
            obj.guihan.MeasurementPlaneDropDown.Value = "Vertical" ;
          end
          if obj.ProfName == "WIRE:LI18:944"
            obj.guihan.LinacDropDown.Value = "L3" ;
          else
            obj.guihan.LinacDropDown.Value = "L2" ;
          end
          for iwire=1:4
            obj.guihan.("Wire"+iwire+"_Sigma").Value = obj.MultiWireData(dim,iwire) ;
            obj.guihan.("Wire"+iwire+"_SigmaErr").Value = obj.MultiWireDataErr(dim,iwire) ;
          end
        end
        didload=2;
        return
      end
      
      % The rest is dealing with an actual quad scan...
      if startsWith(string(fn),'Emittance') % data from emittance_gui?
        iscorplot=false;
      else % data from corr plot
        iscorplot=true;
      end
      if iscorplot
        xpv = find(endsWith(string(arrayfun(@(x) x.name,dat.data.profPV(:,1),'UniformOutput',false)),"XRMS"),1) ;
        obj.ProfName = string(regexprep(dat.data.profPV(xpv).name,':XRMS$','')) ; % also updates model
        obj.QuadScanData.QuadName = string(regexprep(dat.data.ctrlPV(1).name,'(:BCTRL)|(:BDES)$','')) ;
      else
        obj.ProfName = string(dat.data.name{1}) ;
        obj.QuadScanData.QuadName = string(dat.data.quadName) ;
      end
      obj.QuadScanData.ProfInd = find(LM_all.ControlNames==obj.ProfName,1) ;
      if iscorplot
        stat = logical([dat.data.status]) ;
      else
        stat = logical([dat.data.use]) ;
      end
      if iscorplot
        obj.QuadScanData.QuadVal = [dat.data.ctrlPV(stat).val] ;
      else
        obj.QuadScanData.QuadVal = dat.data.quadVal(stat) ;
      end
      obj.QuadScanData.nscan = sum(stat) ;
      fitm=["Gaussian" "Asymmetric"];
      if iscorplot
        sz=size(dat.data.beam);
      else
        sz=size(dat.data.beamList) ;
      end
      % Check for existence of X/Y data
      if ~isfield(dat.data.beam(1),'xStat')
        obj.DimSelect="Y";
        if ~isempty(obj.guihan)
          obj.guihan.UseYDataButton.Value=true;
          obj.guihan.UseXDataButton.Value=false;
          drawnow
        end
      end
      if ~isfield(dat.data.beam(1),'yStat')
        if obj.DimSelect=="Y"
          error('No data to load');
        end
        obj.DimSelect="X";
        if ~isempty(obj.guihan)
          obj.guihan.UseXDataButton.Value=true;
          obj.guihan.UseYDataButton.Value=false;
          drawnow
        end
      end
      for ii=1:2
        if iscorplot
          if contains(obj.DimSelect,"X")
            xdat = reshape([dat.data.beam(stat,:,ii).xStat],5,sum(stat),sz(2)) ;
            xerr = reshape([dat.data.beam(stat,:,ii).xStatStd],5,sum(stat),sz(2)) ;
          end
          if contains(obj.DimSelect,"Y")
            ydat = reshape([dat.data.beam(stat,:,ii).yStat],5,sum(stat),sz(2)) ;
            yerr = reshape([dat.data.beam(stat,:,ii).yStatStd],5,sum(stat),sz(2)) ;
          end
        else
          if contains(obj.DimSelect,"X")
            xdat = reshape([dat.data.beamList(stat,:,ii).xStat],5,sum(stat),sz(2)) ;
            xerr = reshape([dat.data.beamList(stat,:,ii).xStatStd],5,sum(stat),sz(2)) ;
          end
          if contains(obj.DimSelect,"Y")
            ydat = reshape([dat.data.beamList(stat,:,ii).yStat],5,sum(stat),sz(2)) ;
            yerr = reshape([dat.data.beamList(stat,:,ii).yStatStd],5,sum(stat),sz(2)) ;
          end
        end
        if sz(2)>1
          if contains(obj.DimSelect,"X")
            xdat_err=std(squeeze(xdat(3,:,:)),[],2);
            xdat=mean(squeeze(xdat(3,:,:)),2);
          end
          if contains(obj.DimSelect,"Y")
            ydat_err=std(squeeze(ydat(3,:,:)),[],2);
            ydat=mean(squeeze(ydat(3,:,:)),2);
          end
        else
          if contains(obj.DimSelect,"X")
            xdat=squeeze(xdat(3,:,:));
            xdat_err=squeeze(xerr(3,:,:));
          end
          if contains(obj.DimSelect,"Y")
            ydat=squeeze(ydat(3,:,:));
            ydat_err=squeeze(yerr(3,:,:));
          end
        end
        if contains(obj.DimSelect,"X")
          obj.QuadScanData.x.(fitm(ii)) = xdat ;
          obj.QuadScanData.xerr.(fitm(ii)) = xdat_err ;
        end
        if contains(obj.DimSelect,"Y")
          obj.QuadScanData.y.(fitm(ii)) = ydat ;
          obj.QuadScanData.yerr.(fitm(ii)) = ydat_err ;
        end
      end
      obj.QuadScanData.QuadInd = find(LM_mags.ControlNames==obj.QuadScanData.QuadName,1) ;
      if isempty(obj.QuadScanData.QuadInd)
        errordlg(sprintf('Scan Quad (%s) Not found in model',obj.QuadScanData.QuadName),'QUAD not found');
        obj.QuadScanData=[];
        return
      end
      didload=1;
    end
    function WriteEmitData(obj,dim)
      %WRITEENMITDATA Write twiss and emittance data to PVs
      %WriteEmitData([dim])
      % dim = "x" | "y" | "xy" (default="xy")
      if ismember(obj.ProfName,obj.MWNames_L2)
        pn = obj.MWNames_L2 ;
      elseif ismember(obj.ProfName,obj.MWNames_L3)
        pn = obj.MWNames_L3 ;
      else
        pn= obj.ProfName ;
      end
      if ~exist('dim','var')
        dim="xy";
      else
        dim=lower(string(dim));
      end
      for ipn=1:length(pn)
        if contains(dim,"x")
          lcaPutNoWait(sprintf('%s:BETA_X',pn(ipn)),obj.TwissFit(1));
          lcaPutNoWait(sprintf('%s:ALPHA_X',pn(ipn)),obj.TwissFit(3));
          lcaPutNoWait(sprintf('%s:BMAG_X',pn(ipn)),obj.TwissFit(5));
          lcaPutNoWait(sprintf('%s:EMITN_X',pn(ipn)),obj.TwissFit(7));
        end
        if contains(dim,"y")
          lcaPutNoWait(sprintf('%s:BETA_Y',pn(ipn)),obj.TwissFit(2));
          lcaPutNoWait(sprintf('%s:ALPHA_Y',pn(ipn)),obj.TwissFit(4));
          lcaPutNoWait(sprintf('%s:BMAG_Y',pn(ipn)),obj.TwissFit(6));
          lcaPutNoWait(sprintf('%s:EMITN_Y',pn(ipn)),obj.TwissFit(8));
        end
        % Write initial match conditions to PVs
%         if ~isempty(obj.InitMatch) && ismember(obj.ProfName,obj.InitMatchProf) && obj.goodmatch
%           lcaPutNoWait(char(obj.LiveModel.Initial_betaxPV),obj.InitMatch.x.Twiss.beta) ;
%           lcaPutNoWait(char(obj.LiveModel.Initial_betaxPV),obj.InitMatch.x.Twiss.alpha) ;
%           lcaPutNoWait(char(obj.LiveModel.Initial_emitxPV),obj.InitMatch.x.NEmit*1e6) ;
%           lcaPutNoWait(char(obj.LiveModel.Initial_betaxPV),obj.InitMatch.y.Twiss.beta) ;
%           lcaPutNoWait(char(obj.LiveModel.Initial_betaxPV),obj.InitMatch.y.Twiss.alpha) ;
%           lcaPutNoWait(char(obj.LiveModel.Initial_emityPV),obj.InitMatch.y.NEmit*1e6) ;
%         end
      end
      obj.InitRestore = obj.LiveModel.Initial ;
      if ~isempty(obj.InitMatch)
        obj.LiveModel.Initial=obj.InitMatch ;
      end
    end
    function ReadEmitData(obj)
      if ~ismember(obj.ProfName,obj.EmitDataProfs) && ~obj.LoadingScanData
        [betax,pvtime] = lcaGet(sprintf('%s:BETA_X',obj.ProfName));
        betay = lcaGet(sprintf('%s:BETA_Y',obj.ProfName));
        alphax = lcaGet(sprintf('%s:ALPHA_X',obj.ProfName));
        alphay = lcaGet(sprintf('%s:ALPHA_Y',obj.ProfName));
        bmagx = lcaGet(sprintf('%s:BMAG_X',obj.ProfName));
        bmagy = lcaGet(sprintf('%s:BMAG_Y',obj.ProfName));
        emitnx = lcaGet(sprintf('%s:EMITN_X',obj.ProfName));
        emitny = lcaGet(sprintf('%s:EMITN_Y',obj.ProfName));
        obj.TwissFitModel=[betax betay alphax alphay bmagx bmagy emitnx emitny];
        obj.TwissFitAnalytic=[betax betay alphax alphay bmagx bmagy emitnx emitny];
        % Set achive date to match PV write date if not already set by reading quad values
        if ~isfield(obj.QuadScanData,'QuadName')
          obj.ModelDate = datevec(F2_common.epics2mltime(pvtime)) ;
          obj.LiveModel.ArchiveDate = obj.ModelDate ;
        end
      end
    end
    function FitQuadScanData(obj)
      %FITQUADSCANDATA Fit twiss parameter to live model using quad scan data
      
      global BEAMLINE PS
      % Find scan quad ID, record PS setting and set to zero for initial calculations
      DesignTwiss = obj.LiveModel.DesignTwiss;
      qid = obj.QuadScanData.QuadInd ;
      LM_mags=obj.LiveModel.LEM.Mags.LM;
      pele=obj.ProfModelInd;
      qele = LM_mags.ModelID(qid);
      ips = BEAMLINE{qele}.PS;
      ps0 = PS(ips).Ampl;
      k = obj.quadscan_k ;
      if contains(obj.DimSelect,"X")
        x = obj.QuadScanData.x.(obj.ProfFitMethod).*1e-6 ; 
        xerr = obj.QuadScanData.xerr.(obj.ProfFitMethod).*1e-6 ; 
      end
      if contains(obj.DimSelect,"Y")
        y = obj.QuadScanData.y.(obj.ProfFitMethod).*1e-6 ; 
        yerr = obj.QuadScanData.yerr.(obj.ProfFitMethod).*1e-6 ; 
      end
      rgamma = LM_mags.ModelP(qid)/0.511e-3;
      
      % Get Twiss at profile device using Model (RmatAtoB)
      qscanvals=obj.QuadScanData.QuadVal./10; % Quad scan values / T
      Rscan=cell(1,length(qscanvals));
      for iscan=1:length(qscanvals)
        PS(ips).Ampl=qscanvals(iscan);
        [~,R]=RmatAtoB(qele,pele);
        Rscan{iscan}=R;
      end
      PS(ips).Ampl=ps0;
      [~,R]=RmatAtoB(qele,pele);
      obj.ModelQuadScan=zeros(2,length(qscanvals));
      if contains(obj.DimSelect,"X")
        xopt = lsqnonlin(@(xx) obj.ModelTwissFitFn(xx,1:2,Rscan,x,xerr),...
          [1 0 obj.LiveModel.Initial.x.NEmit/rgamma],[0.01 -100 0.1e-6/rgamma],...
          [1e4 100 1e-3/rgamma],optimset('Display','iter'));
        emitx = xopt(3) ;
        Sx = emitx .* [xopt(1) -xopt(2);-xopt(2) (1+xopt(2)^2)/xopt(1)] ;
        Sx_prof = R(1:2,1:2) * Sx * R(1:2,1:2)' ; Sx_prof=Sx_prof./emitx ;
        obj.TwissFitModel([1 3 7]) = [Sx_prof(1,1) -Sx_prof(1,2) emitx*rgamma*1e6] ;
        obj.TwissFitModel(5) = bmag(obj.LiveModel.DesignTwiss.betax(pele),obj.LiveModel.DesignTwiss.alphax(pele),...
          Sx_prof(1,1),-Sx_prof(1,2));
        for iscan=1:length(qscanvals)
          S_prof = Rscan{iscan}(1:2,1:2) * Sx * Rscan{iscan}(1:2,1:2)' ; % Sigma matrix at profile monitor location
          obj.ModelQuadScan(1,iscan) = sqrt(S_prof(1,1)).*1e6 ;
        end
      end
      if contains(obj.DimSelect,"Y")
        yopt = lsqnonlin(@(xx) obj.ModelTwissFitFn(xx,3:4,Rscan,y,yerr),...
          [1 0 obj.LiveModel.Initial.x.NEmit/rgamma],[0.01 -100 0.1e-6/rgamma],...
          [1e4 100 1e-3/rgamma],optimset('Display','iter'));
        emity = yopt(3) ;
        Sy = emity .* [yopt(1) -yopt(2);-yopt(2) (1+yopt(2)^2)/yopt(1)] ;
        Sy_prof = R(3:4,3:4) * Sy * R(3:4,3:4)' ; Sy_prof = Sy_prof ./ emity ;
        obj.TwissFitModel([2 4 8]) = [Sy_prof(1,1) -Sy_prof(1,2) emity*rgamma*1e6] ;
        obj.TwissFitModel(6) = bmag(obj.LiveModel.DesignTwiss.betay(pele),obj.LiveModel.DesignTwiss.alphay(pele),...
          Sy_prof(1,1),-Sy_prof(1,2));
        for iscan=1:length(qscanvals)
          S_prof = Rscan{iscan}(3:4,3:4) * Sy * Rscan{iscan}(3:4,3:4)' ; % Sigma matrix at profile monitor location
          obj.ModelQuadScan(2,iscan) = sqrt(S_prof(1,1)).*1e6 ;
        end
      end
      
      % Twiss parameters using analytic approach
      %https://indico.cern.ch/event/703517/contributions/2886144/attachments/1602810/2541739/2018-02-19_Emittance_measurement_with_quadrupole_scan.pdf
      PS(ips).Ampl=0;
      if contains(obj.DimSelect,"X")
        qx=noplot_polyfit(-k,x.^2,2.*x.*xerr,2);
      end
      if contains(obj.DimSelect,"Y")
        qy=noplot_polyfit(k,y.^2,2.*y.*yerr,2);
      end
      % sigma matrix at (thin) quad entrance
      [~,R]=RmatAtoB(qele+1,pele); % from center of quad to profile monitor
      % -- x
      if contains(obj.DimSelect,"X")
        d=R(1,2);
        Lq=sum(arrayfun(@(x) BEAMLINE{x}.L,PS(ips).Element));
        A=qx(3); B=qx(2); C=qx(1);
        sig11 = A / (d^2*Lq^2) ;
        sig12 = (B-2*d*Lq*sig11)/(2*d^2*Lq) ;
        sig22 = (C-sig11-2*d*sig12) / d^2 ;
        emit = sqrt(sig11*sig22 - sig12^2) ;
        alpha = - sig12/emit ;
        beta = sig11/emit ;
        % Propogate Twiss back to start of Model quadrupole and then forward to profile monitor location
        S = emit .* [beta -alpha;-alpha (1+alpha^2)/beta] ;
        Rd = [1 Lq/2;0 1] ; S(1,2)=-S(1,2); S(2,1)=-S(2,1);
        S0 = Rd*S*Rd';
        PS(ips).Ampl=ps0;
        [~,R2]=RmatAtoB(qele,pele);
        S0(1,2)=-S0(1,2); S0(2,1)=-S0(2,1);
        S_prof = R2(1:2,1:2)*S0*R2(1:2,1:2)';
        T=S_prof./emit ;
        if ~isreal(T)
          obj.TwissFitAnalytic(2:2:end)=nan;
        else
          obj.TwissFitAnalytic(7) = rgamma * emit * 1e6 ;
          obj.TwissFitAnalytic(1) = T(1,1) ;
          obj.TwissFitAnalytic(3) = -T(1,2) ;
          obj.TwissFitAnalytic(5) = bmag(DesignTwiss.betax(pele),DesignTwiss.alphax(pele),...
            obj.TwissFitAnalytic(1),obj.TwissFitAnalytic(3));
        end
      end
      % -- y
      if contains(obj.DimSelect,"Y")
        d=R(3,4);
        Lq=sum(arrayfun(@(x) BEAMLINE{x}.L,PS(ips).Element));
        A=qy(3); B=qy(2); C=qy(1);
        sig11 = A / (d^2*Lq^2) ;
        sig12 = (B-2*d*Lq*sig11)/(2*d^2*Lq) ;
        sig22 = (C-sig11-2*d*sig12) / d^2 ;
        emit = sqrt(sig11*sig22 - sig12^2) ;

        alpha = - sig12/emit ;
        beta = sig11/emit ;
        % Propogate Twiss back to start of Model quadrupole and then forward to profile monitor location
        S = emit .* [beta -alpha;-alpha (1+alpha^2)/beta] ; S(1,2)=-S(1,2); S(2,1)=-S(2,1);
        Rd = [1 Lq/2;0 1] ;
        S0 = Rd*S*Rd'; S0(1,2)=-S0(1,2); S0(2,1)=-S0(2,1);
        PS(ips).Ampl=ps0;
        [~,R2]=RmatAtoB(qele,pele);
        S_prof = R2(3:4,3:4)*S0*R2(3:4,3:4)';
        T=S_prof./emit ;
        if ~isreal(T)
          obj.TwissFitAnalytic(2:2:end)=nan;
        else
          obj.TwissFitAnalytic(8) = rgamma * emit * 1e6 ;
          obj.TwissFitAnalytic(2) = T(1,1) ;
          obj.TwissFitAnalytic(4) = -T(1,2) ;
          obj.TwissFitAnalytic(6) = bmag(obj.LiveModel.DesignTwiss.betay(pele),obj.LiveModel.DesignTwiss.alphay(pele),...
            obj.TwissFitAnalytic(2),obj.TwissFitAnalytic(4));
        end
      end
      
      % Update initial Twiss match
      try
        obj.DoMatch("init");
      catch
        fprintf(2,"Failed to match initial Twiss Parameters!\n");
      end
    end
    function fhan = PlotQuadScanData(obj,newfig)
      fhan=[];
      if ~exist('newfig','var')
        newfig=false;
      end
      if ~isempty(obj.guihan) && ~newfig
        obj.guihan.UIAxes.reset;
        obj.guihan.UIAxes_2.reset;
        cla(obj.guihan.UIAxes);
        cla(obj.guihan.UIAxes_2);
      end
      if ~isfield(obj.QuadScanData,'x') && ~isfield(obj.QuadScanData,'y')
        return
      end
      k = obj.quadscan_k ;
      if contains(obj.DimSelect,"X")
        x = obj.QuadScanData.x.(obj.ProfFitMethod) ; 
        xerr = obj.QuadScanData.xerr.(obj.ProfFitMethod) ; 
      end
      if contains(obj.DimSelect,"Y")
        y = obj.QuadScanData.y.(obj.ProfFitMethod) ; 
        yerr = obj.QuadScanData.yerr.(obj.ProfFitMethod) ; 
      end
      if isempty(obj.guihan) || newfig
        fhan=figure;
        ah1=subplot(2,1,1);
        ah2=subplot(2,1,2);
      else
        ah1=obj.guihan.UIAxes;
        ah2=obj.guihan.UIAxes_2;
      end
      ah1.reset; ah2.reset;
      kabs=abs(k);
      if contains(obj.DimSelect,"X")
        [q,dq]=noplot_polyfit(k,x.^2,2.*x.*xerr,2);
        errorbar(ah1,kabs,x.^2,2.*x.*xerr,'ko'); ax=axis(ah1);
        hold(ah1,'on');
        k_fit=linspace(k(1),k(end),1000);
        ypl1=(q(1)-dq(1))+(q(2)-dq(2)).*k_fit+(q(3)-dq(3)).*k_fit.^2; 
        ypl2=(q(1)+dq(1))+(q(2)+dq(2)).*k_fit+(q(3)+dq(3)).*k_fit.^2;
        ypl=[ypl1(:) ypl2(:)-ypl1(:)];
        apl=area(ah1,abs(k_fit),ypl); apl(1).FaceColor='none'; apl(1).LineStyle='none'; apl(2).FaceColor=[0.3010 0.7450 0.9330]; apl(2).LineStyle='none'; apl(2).FaceAlpha=0.5;
        if length(obj.ModelQuadScan)==length(k)
          plot(ah1,abs(k_fit),interp1(k,obj.ModelQuadScan(1,:).^2,k_fit,'spline'),'r','LineWidth',2);
          if obj.ShowPlotLegend; legend(ah1,["Data" "" "Polynomial Fit" "Model Fit"]); else; legend(ah1,'off'); end
        else
          if obj.ShowPlotLegend; legend(ah1,["Data" "" "Polynomial Fit"]); else; legend(ah1,'off'); end
        end
        xlabel(ah1,sprintf('|k| (%s)',obj.QuadScanData.QuadName)); ylabel(ah1,'\sigma_x^2 [\mum^2]');
        axis(ah1,ax);
        grid(ah1,'on');
        hold(ah1,'off');
      end
      if contains(obj.DimSelect,"Y")
        [q,dq]=noplot_polyfit(k,y.^2,2.*y.*yerr,2);
        errorbar(ah2,kabs,y.^2,2.*y.*yerr,'ko'); ax=axis(ah2);
        hold(ah2,'on');
        k_fit=linspace(k(1),k(end),1000);
        ypl1=(q(1)-dq(1))+(q(2)-dq(2)).*k_fit+(q(3)-dq(3)).*k_fit.^2; 
        ypl2=(q(1)+dq(1))+(q(2)+dq(2)).*k_fit+(q(3)+dq(3)).*k_fit.^2;
        ypl=[ypl1(:) ypl2(:)-ypl1(:)];
        apl=area(ah2,abs(k_fit),ypl); apl(1).FaceColor='none'; apl(1).LineStyle='none'; apl(2).FaceColor=[0.3010 0.7450 0.9330]; apl(2).LineStyle='none'; apl(2).FaceAlpha=0.5;
        if length(obj.ModelQuadScan)==length(k)
          plot(ah2,abs(k_fit),interp1(k,obj.ModelQuadScan(2,:).^2,k_fit,'spline'),'r','LineWidth',2);
          if obj.ShowPlotLegend; legend(ah2,["Data" "" "Polynomial Fit" "Model Fit"]); else; legend(ah2,'off'); end
        else
          if obj.ShowPlotLegend; legend(ah2,["Data" "" "Polynomial Fit"]); else; legend(ah2,'off'); end
        end
        xlabel(ah2,sprintf('|k| (%s)',obj.QuadScanData.QuadName)); ylabel(ah2,'\sigma_y^2 [\mum^2]');
        axis(ah2,ax);
        grid(ah2,'on');
        hold(ah2,'off');
      end
      if ~isempty(obj.guihan) % Fill emit / bmag values for Model and Analytic fit results
        TA=obj.TwissFitAnalytic; TA(isnan(TA))=0;
        obj.guihan.EditField_2.Value = TA(7) ;
        obj.guihan.EditField.Value = TA(5) ;
        obj.guihan.EditField_4.Value = obj.TwissFitModel(7) ;
        obj.guihan.EditField_3.Value = obj.TwissFitModel(5) ;
        obj.guihan.EditField_6.Value = TA(8) ;
        obj.guihan.EditField_5.Value = TA(6) ;
        obj.guihan.EditField_8.Value = obj.TwissFitModel(8) ;
        obj.guihan.EditField_7.Value = obj.TwissFitModel(6) ;
      end
    end
    function fhan=PlotTwiss(obj,newfig)
      global BEAMLINE PS
      fhan=[];
      if ~exist('newfig','var')
        newfig=false;
      end
      if isempty(obj.guihan) || newfig
        fhan=figure;
        ah=axes;
      else
        ah=obj.guihan.UIAxes2;
      end
      ah.reset;
      LM_mags=obj.LiveModel.LEM.Mags.LM;
      if obj.MatchDir=="BKW"
        iquads = LM_mags.ModelID(:) < obj.ProfModelInd & obj.LEMQuadID(:) ;
      else
        iquads = LM_mags.ModelID(:) > obj.ProfModelInd & obj.LEMQuadID(:) ;
      end
      if obj.MatchDir=="BKW"
        id1=find(iquads,double(obj.NumMatchQuads),'last') ; id1=id1(obj.UseMatchQuad);
      else
        id1=find(iquads,double(obj.NumMatchQuads),'first') ; id1=id1(obj.UseMatchQuad);
      end
      quadps = arrayfun(@(x) BEAMLINE{x}.PS,LM_mags.ModelID(id1)) ;
      if obj.MatchDir=="BKW"
        i1=PS(quadps(1)).Element(1);
        i2=obj.ProfModelInd;
      else
        i1=obj.ProfModelInd;
        i2=PS(quadps(end)).Element(end);
      end
      z=arrayfun(@(x) BEAMLINE{x}.Coordi(3),i1:i2);
      txt=string.empty;
      if ~isempty(obj.TwissPreMatch)
        plot(ah,obj.TwissPreMatch.z(1:end-1),obj.TwissPreMatch.betax(1:end-1),'Color',[0 0.4470 0.7410],'LineStyle','-.'); hold(ah,'on');
        txt(end+1)="\beta_x (current)";
      end
      if ~isempty(obj.TwissMatch)
        plot(ah,obj.TwissMatch.z(1:end-1),obj.TwissMatch.betax(1:end-1),'Color',[0 0.4470 0.7410],'LineStyle','--'); hold(ah,'on');
        txt(end+1)="\beta_x (matched)";
      end
      plot(ah,z,obj.LiveModel.DesignTwiss.betax(i1:i2),'Color',[0 0.4470 0.7410]); hold(ah,'on');
      grid(ah,'on');
      txt(end+1)="\beta_x (design)";
      xlabel(ah,'Z [m]'); ylabel(ah,'\beta [m]');
      if ~isempty(obj.TwissPreMatch)
        plot(ah,obj.TwissPreMatch.z(1:end-1),obj.TwissPreMatch.betay(1:end-1),'Color',[0.8500 0.3250 0.0980],'LineStyle','-.'); hold(ah,'on');
        txt(end+1)="\beta_y (current)";
      end
      if ~isempty(obj.TwissMatch)
        plot(ah,obj.TwissMatch.z(1:end-1),obj.TwissMatch.betay(1:end-1),'Color',[0.8500 0.3250 0.0980],'LineStyle','--'); hold(ah,'on');
        txt(end+1)="\beta_y (matched)";
      end
      plot(ah,z,obj.LiveModel.DesignTwiss.betay(i1:i2),'Color',[0.8500 0.3250 0.0980]); hold(ah,'off');
      axis(ah,'tight'); 
      txt(end+1)="\beta_y (design)";
      if obj.ShowPlotLegend; legend(ah,txt) ; else; legend(ah,'off'); end
      if ~isempty(obj.guihan) && ~newfig
        ax=axis(ah);
        axis(obj.guihan.UIAxes2_2,ax);
        AddMagnetPlotZ(i1,i2,obj.guihan.UIAxes2_2,'replace');
      else
        AddMagnetPlotZ(i1,i2);
      end
    end
    function tab = TwissTable(obj)
      param = ["beta_x";"beta_y";"alpha_x"; "alpha_y";"bmag_x";"bmag_y";"nemit_x";"nemit_y"] ;
      meas = obj.TwissFit ;
      if isfield(obj.QuadScanData,'twissmatch')
        match = [obj.QuadScanData.twissmatch nan nan nan nan];
      else
        match = nan(1,8) ;
      end
      dt=obj.LiveModel.DesignTwiss;
      pi=obj.ProfModelInd;
      I=obj.LiveModel.DesignInitial ;
      design = [dt.betax(pi) dt.betay(pi) dt.alphax(pi) dt.alphay(pi) 1 1 I.x.NEmit*1e6 I.y.NEmit*1e6] ;
      tab=table(param(:),meas(:),match(:),design(:));
      tab.Properties.VariableNames=["Param";"Meas.";"Match";"Design"];
    end
    function tab = MagnetTable(obj)
      global BEAMLINE
      Z = arrayfun(@(x) BEAMLINE{x}.Coordi(3),obj.MatchQuadModelInd) ;
      E = arrayfun(@(x) BEAMLINE{x}.P,obj.MatchQuadModelInd) ;
      obj.LM.ModelClasses="QUAD";
      bdes = obj.LiveModel.LEM.Mags.LM.ModelBDES(ismember(obj.LiveModel.LEM.Mags.LM.ModelID,obj.MatchQuadModelInd)) ;
      if obj.ModelSource=="Design"
        bact = bdes ;
      else
        bact = obj.LiveModel.LEM.Mags.BDES_cntrl(ismember(obj.LiveModel.LEM.Mags.LM.ModelID,obj.MatchQuadModelInd)) ;
      end
      bmatch = nan(size(bdes)) ;
      if isfield(obj.QuadScanData,'quadmatch')
        bmatch(obj.UseMatchQuad) = obj.QuadScanData.quadmatch ;
      end
      b_design = arrayfun(@(x) obj.LiveModel.LM.DesignBeamline{x}.B*20,obj.MatchQuadModelInd) ;
      usequad = obj.UseMatchQuad ;
      tab=table(obj.MatchQuadNames(:),string(num2str(Z(:),6)),E(:),bdes(:),bact(:),bmatch(:),b_design(:),usequad(:));
      tab.Properties.VariableNames=["Name";"Z";"E (GeV)";"BDES";"BACT";"BMATCH";"B_DESIGN";"USE"];
    end
    function [emitData,txt_results] = emitMW(obj,dim,data,data_err,section)
      %EMITMW Multi-Wire emittance calculation
      %emitData = emitMW(dim,data,section)
      %
      % Compute beam sigma matrix, covariance matrix, and chisquare from measured
      % beam sizes and their errors at 3 or more wires
      %
      % INPUTs:
      %
      %   dim      = "X" or "Y"
      %   data     = 1x4 vector of measured wire scan rms sizes [m] (zeros or nans indicates not to use data for this wire)
      %   data_err = 1x4 vector of measured wire scan rms size errors [m]
      %   section  = "L2" or "L3"
      %
      % OUTPUTs:
      %
      %   emitData    = return data structure
      %   txt_results = text of results
      
      global BEAMLINE
      
      data=data(:); data_err=data_err(:);
      
      % Treat all zero errors as no error info
      if isempty(data_err)
        data_err=zeros(1,4);
      end
      data_err(data_err==0)=1e-6;
      
      % Which wires to include?
      wsel=true(4,1);
      wsel(data==0 | isnan(data))=false;
      if sum(wsel)<3
        error('Need min 3 wires');
      end
      nw=sum(wsel);
      
      % Get wire names
      switch section
        case "L2"
          wname={'WS11444' 'WS11614' 'WS11744' 'WS12214'};
          desemit=[lcaGet('PROF:IN10:571:EMITN_X') lcaGet('PROF:IN10:571:EMITN_Y')].*1e-6 ;
        case "L3"
          wname={'WS18944' 'WS19144' 'WS19244' 'WS19344'};
          desemit=[lcaGet('WIRE:LI11:614:EMITN_X') lcaGet('WIRE:LI11:614:EMITN_Y')].*1e-6 ;
        otherwise
          error('Incorrect section selection');
      end
      wname=wname(wsel);
      id1=findcells(BEAMLINE,'Name',wname{1});
      
      % get pointers
      idw=zeros(1,nw);
      for n=1:nw
        idw(n)=findcells(BEAMLINE,'Name',wname{n});
      end
      
      % get the model
      R=cell(1,nw);
      Rx=zeros(nw,2);Ry=zeros(nw,2);
      for n=1:nw
        [stat,Rab]=RmatAtoB(id1,idw(n));
        if (stat{1}~=1),error(stat{2}),end
        R{n}=Rab(1:4,1:4);
        Rx(n,:)=[Rab(1,1),Rab(1,2)];
        Ry(n,:)=[Rab(3,3),Rab(3,4)];
      end
      
      % get design Twiss at first wire
      energy=BEAMLINE{id1}.P;
      egamma = energy/0.51099906e-3;
      if dim=="X"
        bx0 = obj.LiveModel.DesignTwiss.betax(id1) ;
        ax0 = obj.LiveModel.DesignTwiss.alphax(id1) ;
      else
        bx0 = obj.LiveModel.DesignTwiss.betay(id1);
        ax0 = obj.LiveModel.DesignTwiss.alphay(id1);
      end
      
      % load analysis variables
      x=data(wsel).^2 ; % m^2
      dx=2.*x.*data_err(wsel) ;
      
      % compute least squares solution
      M=zeros(nw,3);
      for n=1:nw
        if dim=="X"
          M(n,1)=Rx(n,1)^2;
          M(n,2)=2*Rx(n,1)*Rx(n,2);
          M(n,3)=Rx(n,2)^2;
        else
          M(n,1)=Ry(n,1)^2;
          M(n,2)=2*Ry(n,1)*Ry(n,2);
          M(n,3)=Ry(n,2)^2;
        end
      end
      
      % Solve least-squares problem, weighted by errors
      for itry=1:2
        zx=x./dx;
        Bx=zeros(nw,3);
        for n=1:nw
          Bx(n,:)=M(n,:)./dx(n);
        end
        Tx=inv(Bx'*Bx);
        u=Tx*Bx'*zx;du=sqrt(diag(Tx));  %#ok<NASGU,MINV>
        if itry==1
          chi2x = sum( (Bx*u - zx).^2 ./ zx ) ;
          dx=dx.*sqrt(chi2x);
        end
      end
      
      % convert fitted input sigma matrix elements to emittance, BMAG, ...
      [px,dpx]=obj.emit_params(u(1),u(2),u(3),Tx,bx0,ax0);
      px(1:3)=abs(px(1:3));
      emitx=px(1);demitx=dpx(1);
      bmagx=px(2);dbmagx=dpx(2);
      betax=px(4);dbetax=dpx(4);
      alphx=px(5);dalphx=dpx(5);
      bcosx=px(6);dbcosx=dpx(6);
      bsinx=px(7);dbsinx=dpx(7);
      emitxn=egamma*emitx;demitxn=egamma*demitx;
      
      % results txt
      fprintf('%s emittance parameters at %s\n',dim,wname{1});
      fprintf('-----------------------------------------------------\n');
      fprintf('energy      = %10.4f              GeV\n',energy);
      fprintf('nemit       = %10.4f +- %9.4f mm-mrad\n',1e6*emitxn,1e6*demitxn);
      fprintf('nemit*bmag  = %10.4f          mm-mrad\n',1e6*emitxn*bmagx);
      fprintf('emit        = %10.4f +- %9.4f nm\n',1e9*emitx,1e9*demitx);
      fprintf('bmag        = %10.4f +- %9.4f      (%9.4f)\n',bmagx,dbmagx,1);
      fprintf('bmag_cos    = %10.4f +- %9.4f      (%9.4f)\n',bcosx,dbcosx,0);
      fprintf('bmag_sin    = %10.4f +- %9.4f      (%9.4f)\n',bsinx,dbsinx,0);
      fprintf('beta        = %10.4f +- %9.4f m    (%9.4f)\n',betax,dbetax,bx0);
      fprintf('alpha       = %10.4f +- %9.4f      (%9.4f)\n',alphx,dalphx,ax0);
      fprintf('chisq/N     = %10.4f\n',chi2x);
      
      % Return data structure
      emitData.emit=emitx;
      emitData.nemit=emitxn;
      emitData.bmag=bmagx;
      emitData.beta=betax;
      emitData.alpha=alphx;
      emitData.demit=demitx;
      emitData.dnemit=demitxn;
      emitData.dbmag=dbmagx;
      emitData.dbeta=dbetax;
      emitData.dalpha=dalphx;
      emitData.dim=dim;
      emitData.sig=data;
      emitData.dsig=data_err;
      emitData.beta0=bx0;
      emitData.alpha0=ax0;
      if dim=="X"
        emitData.nemit0=desemit(1);
      else
        emitData.nemit0=desemit(2);
      end
      emitData.emit0=emitData.nemit0/egamma;
      emitData.R=R;
      emitData.idw=idw;
      emitData.wsel=wsel;
      I = TwissToInitial(obj.LiveModel.DesignTwiss,idw(1),obj.LiveModel.Initial);
      I.Momentum=energy;
      emitData.I_design = I ;
      emitData.section=section;
      emitData.pcols=obj.plotcols;
      emitData.id1=id1;
      
      
      txt_results{1} =  sprintf('%s emittance parameters at %s\n',dim,wname{1});
      txt_results{2} =  sprintf('----\n');
      txt_results{3} =  sprintf('energy      = %10.4f               GeV\n',energy);
      txt_results{4} =  sprintf('nemit       = %10.4f   (%9.4f) mm-mrad\n',1e6*emitxn,1e6*emitData.nemit0);
      txt_results{5} =  sprintf('nemit*bmag  = %10.4f   (%9.4f) mm-mrad\n',1e6*emitxn*bmagx,1e6*emitData.nemit0);
      txt_results{6} =  sprintf('emit        = %10.4f   (%9.4f) nm-rad\n',1e9*emitx,1e9*emitData.nemit0/egamma);
      txt_results{7} =  sprintf('bmag        = %10.4f   (%9.4f)\n',bmagx,1);
      txt_results{8} =  sprintf('beta        = %10.4f   (%9.4f) m\n',betax,bx0);
      txt_results{9} =  sprintf('alpha       = %10.4f   (%9.4f)\n',alphx,ax0);
      txt_results{10} = sprintf('chisq/N     = %10.4f\n',chi2x);
      
      % Get design spot sizes
      emitData.sigma_des=zeros(1,nw) ;
      gamma = BEAMLINE{id1}.P/0.511e-3 ;
      emit = desemit./gamma;
      if dim=="X"
        S0 = emit(1) * ...
          [obj.LiveModel.DesignTwiss.betax(id1)    -obj.LiveModel.DesignTwiss.alphax(id1);
          -obj.LiveModel.DesignTwiss.alphax(id1)  (1+obj.LiveModel.DesignTwiss.alphax(id1)^2)/obj.LiveModel.DesignTwiss.betax(id1)] ;
      else
        S0 = emit(2) * ...
          [obj.LiveModel.DesignTwiss.betay(id1)    -obj.LiveModel.DesignTwiss.alphay(id1);
          -obj.LiveModel.DesignTwiss.alphay(id1)  (1+obj.LiveModel.DesignTwiss.alphay(id1)^2)/obj.LiveModel.DesignTwiss.betay(id1)] ;
      end
      for iw=1:nw
        [~,R]=RmatAtoB(id1,idw(iw));
        if dim=="X"
          S1 = R(1:2,1:2) * S0 * R(1:2,1:2)' ;
        else
          S1 = R(3:4,3:4) * S0 * R(3:4,3:4)' ;
        end
        emitData.sigma_des(iw) = sqrt(S1(1,1)) ;
      end
      txt_results{11} =  sprintf('----\n');
      txt_results{12} = sprintf('SIGMA_DES   = %10.4f %10.4f %10.4f %10.4f um\n',emitData.sigma_des.*1e6);
      
      % Write to local emittance properties
      if dim=="X"
        obj.TwissFitAnalytic(7) = emitxn * 1e6 ;
        obj.TwissFitAnalytic(1) = betax ;
        obj.TwissFitAnalytic(3) = alphx ;
        obj.TwissFitAnalytic(5) = bmagx ;
      else
        obj.TwissFitAnalytic(8) = emitxn * 1e6 ;
        obj.TwissFitAnalytic(2) = betax ;
        obj.TwissFitAnalytic(4) = alphx ;
        obj.TwissFitAnalytic(6) = bmagx ;
      end
      obj.TwissFitModel = obj.TwissFitAnalytic ;
      
    end
    function logplot(obj,whichplot,Pos1,Pos2) %#ok<INUSD>
      %LOGPLOT Write logbook entry of currently selected Tab
      %logplot(whichplot,PosMainTab,PosTwissTab)
      % whichplot: "Magnets" | "QuadScan" | "Optics" | "MultiWire"
      % provide main tab window and Twiss tab window Position properties
      switch whichplot
        case "Magnets" % need to use "exportapp" -> only in 2020b
%           border=2;
%           Position=Pos1; Position(3:4)=Position(3:4)+border*2;
%           Position(3)=Position(3)+border+Pos2(3);
%           fhan = uifigure; fhan.Position=Position;
%           pos=Position;
%           pos(1:2)=pos(1:2)+border; pos(3:4)=Pos1(3:4);
%           uitable(fhan,'Data',obj.MagnetTable,'Position',pos);
%           pos=Position;
%           pos(1:2)=pos(1:2)+border; pos(1)=pos(1)+border+Pos1(3);
%           pos(3:4)=Pos2(3:4);
%           uitable(fhan,'Data',obj.TwissTable,'Position',pos);
%           titletxt = 'Matching Quad Data' ;
%           logtxt='';
        case "QuadScan"
          fhan=obj.PlotQuadScanData(true);
          TA=obj.TwissFitAnalytic; TA(isnan(TA))=0;
          emitx_anal = TA(7) ;
          bmagx_anal = TA(5) ;
          emitx_mdl = obj.TwissFitModel(7) ;
          bmagx_mdl = obj.TwissFitModel(5) ;
          emity_anal = TA(8) ;
          bmagy_anal = TA(6) ;
          emity_mdl = obj.TwissFitModel(8) ;
          bmagy_mdl = obj.TwissFitModel(6) ;
          titletxt = 'Matching Quad Scan Data' ;
          logtxt = sprintf('nemit_x/bmagx : %g / %g (analytic) %g / %g (model) [mm-mrad]\n',emitx_anal,bmagx_anal,emitx_mdl,bmagx_mdl) ;
          logtxt = [logtxt sprintf('nemit_y/bmagy : %g / %g (analytic) %g / %g (model) [mm-mrad]',emity_anal,bmagy_anal,emity_mdl,bmagy_mdl)] ;
        case "Optics"
          fhan=obj.PlotTwiss(true);
          titletxt = 'Matching Twiss Data';
          logtxt='';
        otherwise
          error('Unknown logplot option');
      end
      util_printLog2020(fhan, 'title',sprintf('%s - %s',titletxt,char(obj.ProfName)),'author','F2_Matching.m','text',logtxt);
      delete(fhan);
    end
    % Get/Set
    function set.UseFudge(obj,use)
      obj.LiveModel.UseFudge = use ;
      obj.UseFudge = use ;
      % Update initial Twiss match
      try
        obj.DoMatch("init");
      catch
        fprintf(2,"Failed to match initial Twiss Parameters!\n");
      end
    end
    function twiss=get.TwissFit(obj)
      switch obj.TwissFitSource
        case "Model"
          twiss=obj.TwissFitModel;
        case "Analytic"
          twiss=obj.TwissFitAnalytic;
      end
    end
    function set.ModelDate(obj,val)
      obj.ModelDate=val;
      obj.ArchiveDate=val;
      obj.LiveModel.ArchiveDate=val;
      if ~isempty(obj.guihan) && obj.ModelSource=="Archive"
        obj.guihan.ModelDateEditField.Value=datestr(val);
      end
    end
    function set.UndoAvailable(obj,val)
      obj.UndoAvailable=val;
      if ~isempty(obj.guihan)
        obj.guihan.UndoButton.Enable=val;
      end
    end
    function set.goodmatch(obj,val)
      obj.goodmatch=val;
      if ~isempty(obj.guihan)
        obj.guihan.SetMatchingQuadsButton.Enable=obj.goodmatch;
      end
    end
    function set.ModelSource(obj,src)
      switch string(src)
        case "Live"
          obj.UseArchive=false;
          if ~isempty(obj.guihan)
            obj.guihan.ModelDateEditField.Value="LIVE";
          end
        case "Archive"
          obj.UseArchive=true;
          if ~isempty(obj.guihan)
            obj.guihan.ModelDateEditField.Value=datestr(obj.ModelDate);
          end
        case "Design"
          if ~isempty(obj.guihan)
            obj.guihan.ModelDateEditField.Value="USE DESIGN";
          end
          obj.UseArchive=false;
      end
      obj.LiveModel.ModelSource=src; % causes Live and Design model updates, Archive waits until explicite UpdateModel call
      obj.ModelSource=src;
      obj.goodmatch=false;
      obj.UndoAvailable=false;
    end
    function kdes = get.quadscan_k(obj)
      if isempty(obj.QuadScanData)
        kdes=[];
        return
      end
      id = obj.QuadScanData.QuadInd ;
      LM_mags=obj.LiveModel.LEM.Mags.LM;
      kdes = obj.QuadScanData.QuadVal./LM_mags.ModelBDES_L(id)./LM_mags.ModelP(id)./LucretiaModel.GEV2KGM ;
    end
    function set.ProfName(obj,name)
      if string(name) == "<Select From Below>" % default startup condition for GUI
        return
      end
      obj.LM.ModelClasses=["PROF" "WIRE"];
      pind = obj.LM.ModelID(obj.LM.ControlNames==string(name)) ;
      if isempty(pind)
        error('Profile monitor not found in model')
      else
        obj.ProfModelInd=pind(1);
        obj.ProfName=string(name);
      end
      obj.QuadScanData=[];
      obj.ReadEmitData; % Reads emittance data from PVs if there are any and sets archive date
      if obj.ModelSource ~= "Design"
        obj.LiveModel.UpdateModel ;
      end
      obj.NumMatchQuads=obj.NumMatchQuads; % triggers filling matching quad data and propogation of measured Twiss parameters
    end
    function set.NumMatchQuads(obj,num)
      obj.NumMatchQuads=num;
      obj.QuadScanData=[];
      obj.LM.ModelClasses="QUAD";
      if obj.MatchDir=="BKW"
        iquads = obj.LM.ModelID < obj.ProfModelInd ;
        idq=find(iquads,double(obj.NumMatchQuads),'last') ;
      else
        iquads = obj.LM.ModelID > obj.ProfModelInd ;
        idq=find(iquads,double(obj.NumMatchQuads),'first') ;
      end
      obj.MatchQuadNames = obj.LM.ControlNames(idq) ;
      obj.MatchQuadModelInd = obj.LM.ModelID(idq) ;
      obj.MatchQuadID = idq ;
      obj.UseMatchQuad = true(1,length(obj.MatchQuadNames)) ;
      % Update initial Twiss match
      try
        obj.DoMatch("init");
      catch
        fprintf(2,"Failed to match initial Twiss Parameters!\n");
      end
    end
  end
  methods
     function emitMW_plot(obj,emitData,writeToLog)
      %EMITMW_PLOT Plot data from multi-wire emittance data analysis
      %emitMW_plot(emitData)
      global BEAMLINE
      
      section=emitData.section;
      
      % Get wire names
      switch section
        case "L2"
          wname={'WS11444' 'WS11614' 'WS11744' 'WS12214'};
        case "L3"
          wname={'WS18944' 'WS19144' 'WS19244' 'WS19344'};
        otherwise
          error('Incorrect section selection');
      end
      
      % --- reconstructed normalized phase space plots
      % Unpack data
      e0=emitData.emit0;
      b0=emitData.beta0;
      a0=emitData.alpha0;
      g0=(1+a0^2)/b0;
      e=emitData.emit;
      b=emitData.beta;
      a=emitData.alpha;
      sig=emitData.sig;
      dsig=emitData.dsig;
      R=emitData.R;
%       bmag=emitData.bmag;
      idw=emitData.idw;
      wsel=emitData.wsel;
      nw=sum(wsel);
      pcols=emitData.pcols;
      id1=emitData.id1;
      
      S = e * [b -a; -a (1+a^2)/b] ;
      S0 = e0 * [b0 -a0; -a0 (1+a0^2)/b0] ;
      
      % Target figure axis
      fhan=figure; fhan.Position(3:4)=[900 500]; subplot(1,2,1);
      
      if emitData.dim=="X"
        txt1='Horizontal';
        ixy=0;
      else
        txt1='Vertical';
        ixy=1;
      end
      
      % plot measured beam ellipse (normalized coordinates)
      [xe,ye]=noplot_ellipse(inv(S0),0,0);
      ye = (a0*xe+b0*ye)./sqrt(b0) ;
      xe = xe./sqrt(b0) ;
      xsca = max(xe(:)); ysca=max(ye(:));
      plot(xe./xsca,ye./ysca,'k--');
      hold on
      [xe,ye]=noplot_ellipse(inv(S),0,0);
      ye = (a0*xe+b0*ye)./sqrt(b0) ;
      xe = xe./sqrt(b0) ;
      plot(xe./xsca,ye./ysca,'k-');
      v=axis;
      axis(sqrt(2)*max(abs(v))*[-1,1,-1,1],'square')
      
      % plot mapped OTR measurement phases
      ioff=2*ixy;
      big=1000*sqrt(e0*g0);
      c=pcols(wsel); 
      ho=zeros(1,nw);
      for n=1:nw
        sigw=sig(n);
        dsigw=dsig(n);
        X1 = R{n}(1+ioff:2+ioff,1+ioff:2+ioff) \ [sigw+dsigw; big];
        X2 = R{n}(1+ioff:2+ioff,1+ioff:2+ioff) \ [sigw+dsigw; -big];
        X3 = R{n}(1+ioff:2+ioff,1+ioff:2+ioff) \ [sigw-dsigw; big]; 
        X4 = R{n}(1+ioff:2+ioff,1+ioff:2+ioff) \ [sigw-dsigw; -big];
        x1=[X1(1);X2(1)]; y1=[X1(2);X2(2)];
        x2=[X3(1);X4(1)]; y2=[X3(2);X4(2)];
        y1 = (a0*x1+b0*y1)./sqrt(b0) ; x1 = x1./sqrt(b0) ;
        y2 = (a0*x2+b0*y2)./sqrt(b0) ; x2 = x2./sqrt(b0) ;
        h=plot(x1./xsca,y1./ysca,c(n),x2./xsca,y2./ysca,c(n));
        ho(n)=h(1);
      end
      hold off
      hor_line(0,'k:'),ver_line(0,'k:')
      
      % title, labels, and legend
      title([txt1,' Phase Space @ ',wname{1}])
      ylabel('Angle')
      xlabel('Position')
      
      % propagated fitted beam size plots
      subplot(1,2,2);
      Z=arrayfun(@(x) BEAMLINE{x}.Coordf(3),id1:idw(end));
      
      sigfit = zeros(1,1+idw(end)-id1);
      for iele=id1:idw(end)
        [~,R]=RmatAtoB(id1,iele);
        if emitData.dim=="X"
          S1 = R(1:2,1:2) * S * R(1:2,1:2)' ;
        else
          S1 = R(3:4,3:4) * S * R(3:4,3:4)' ;
        end
        sigfit(1+iele-id1) = sqrt(S1(1,1)) ;
      end
      
      plot(Z,sigfit*1e6,'b--')
      hold on
      idz=arrayfun(@(x) BEAMLINE{x}.Coordi(3),idw);
      errorbar(idz,sig*1e6,dsig*1e6,'kx','MarkerSize',10);
      ax=axis;
      for iw=1:length(idw)
        line([idz(iw) idz(iw)],ax(3:4),'Color',c(iw),'LineStyle','--');
      end
      hold off
      title(sprintf('%s',section));
      ylabel(sprintf('%s Beam Size [um]',txt1))
      xlabel('Z [m]')
      AddMagnetPlotZ(idw(1),idw(end));
      
      if exist('writeToLog','var')
        util_printLog2020(fhan, 'title',sprintf('%s Multi-Wire Emittance Measurement (%c)',section,emitData.dim),'author','F2_Matching.m','text',writeToLog);
        obj.WriteEmitData(emitData.dim);
        close(fhan);
      end
     end
  end
  methods(Static,Hidden)
    function [p,dp] = emit_params(sig11,sig12,sig22,C,b0,a0)
      
      %       [p,dp] = emit_params(sig11,sig12,sig22,C,b0,a0);
      %
      %       Returns emittance, bmag, emit*bmag, beta, alpha, and errors on
      %       all these given fitted sigma11(or 33), sigma12(or 34),
      %       sigma22(or 44) and the 3X3 covariance matrix of this fit.
      %
      %     INPUTS:   sig11:  The 1,1 (3,3) sigma matrix element (in m^2-rad)
      %               sig12:  The 1,2 (3,4) sigma matrix element (in m-rad)
      %               sig22:  The 2,2 (4,4) sigma matrix element (in rad^2)
      %               C:      The 3X3 covariance matrix of the above fitted
      %                       sig11,12,22 (33,34,44) (in squared units of
      %                       the 3 above sigij's)
      %               b0:     The matched (design) beta function (in meters)
      %               a0:     The matched (design) alpha function (unitless)
      %     OUTPUTS:  p:      p(1) = unnormalized emittance in (meter-radians)
      %                       p(2) = bmag (unitless beta beat magnitude)
      %                       p(3) = unnormalized emit*bmag (meter-radians)
      %                       p(4) = beta function (in meters)
      %                       p(5) = alpha function (unitless)
      %                       p(6) = (BMAG^2-1)*cos(phi)/BMAG
      %                       p(7) = (BMAG^2-1)*sin(phi)/BMAG
      %               dp:     Measurement errors on above p(1),p(2),...
      
      %===============================================================================
      
      sig0 = [ b0         -a0
        -a0 (1+a0^2)/b0];
      
      sig  = [sig11 sig12
        sig12 sig22];
      
      ebm = 0.5*trace(inv(sig0)*sig); %#ok<MINV>
      
      e2  =  sig11*sig22 - sig12^2;
      if e2 < 0
        e  = -sqrt(abs(e2));
        b  =  sig11/(-e);
        a  = -sig12/(-e);
      else
        e =  sqrt(e2);
        b =  sig11/e;
        a = -sig12/e;
      end
      
      bm  = ebm/e;
      g0  = (1+a0^2)/b0;
      
      bm_cos = (b/b0 - abs(bm))/abs(bm);
      bm_sin = (a - a0*b/b0)/abs(bm);
      
      grad_e   = [ sig22/(2*e)
        -sig12/e
        sig11/(2*e)];
      
      grad_bm  = [(g0/(2*e) - bm*sig22/(2*e^2))
        (a0/e     + bm*sig12/(e^2))
        (b0/(2*e) - bm*sig11/(2*e^2))];
      
      grad_ebm = [g0/2
        a0
        b0/2];
      
      grad_b   = [-0.5*sig11*sig22/(e^3)+(1/e)
        sig11*sig12/(e^3)
        -0.5*sig11*sig11/(e^3)      ];
      
      grad_a   = [ 0.5*sig12*sig22/(e^3)
        -sig12*sig12/(e^3)-(1/e)
        0.5*sig11*sig12/(e^3)      ];
      
      de       = sqrt(grad_e'*C*grad_e);
      dbm      = sqrt(grad_bm'*C*grad_bm);
      debm     = sqrt(grad_ebm'*C*grad_ebm);
      db       = sqrt(grad_b'*C*grad_b);
      da       = sqrt(grad_a'*C*grad_a);
      
      p  = [ e  bm  ebm  b  a bm_cos bm_sin];
      dp = [de dbm debm db da 0 0];
    end
    function opt =  ModelTwissFitFn(x,dims,Rscan,sigma,sigma_err)
      
      S = x(3) .* [x(1) -x(2);-x(2) (1+x(2)^2)/x(1)] ; % Sigma matrix at entrance of scan quad
      sigma_scan = zeros(size(sigma)) ;
      for iscan=1:length(Rscan)
        S_prof = Rscan{iscan}(dims,dims) * S * Rscan{iscan}(dims,dims)' ; % Sigma matrix at profile monitor location
        sigma_scan(iscan) = sqrt(S_prof(1,1)) ;
      end
      opt = (sigma - sigma_scan) ./ sigma_err ;
      
    end
  end
end