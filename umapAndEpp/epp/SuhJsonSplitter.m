classdef SuhJsonSplitter<SuhModalSplitter
    %   AUTHORSHIP
    %   Primary Developer: Stephen Meehan <swmeehan@stanford.edu>
    %   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
    %   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
    %   Provided by the Herzenberg Lab at Stanford University
    %   License: BSD 3 clause
    
    properties(Constant)
        TEST_JSON_FILE=false;
        JSON_FMT=['{"sigma":%s,"min_relative":%s, "goal":"%s",'...
            '"W":%s,"KLD":{"Normal2D":%s,"Normal1D":%s,"Exponential1D":%s}, ' ...
            '"recursive":%s}'];
        MSGID_NO_JSON='SuhJsonSplitter:Incomplete';
        HOME='.EPP';
        DEFAULT_BRANCH='main';%the most sane version
    end
    
    properties(SetAccess=private)
        json;
        startupCost;
        lostChildren;
    end
    
    methods
        function this=SuhJsonSplitter(varargin)
            this=this@SuhModalSplitter(varargin{:});
            if ~this.argued.contains('use_not_gate')
                this.args.use_not_gate=false;
            end
            this.lostChildren=java.util.HashSet;
            this.type='json';
            this.splitsWithPolygons=true;
            this.args.cpp_branch=...
                Args.GetStartsWith('cpp_branch', ...
                SuhJsonSplitter.DEFAULT_BRANCH, varargin);
        end
        
        function [progressMax, txt]=getProgessMax(this, dataSize)
            this.startupCost=.5*dataSize;
            progressMax=dataSize+this.startupCost;
            txt='Running Wayne''s C++ engine for EPP';
        end
        
        function incrementFirstProgress(this, pu)
            pu.incrementProgress(this.startupCost);
        end
        
        function finish(this, ~)
            if this.lostChildren.size>0
                msgWarning(Html.WrapHr(...
                    sprintf(['<html>The JSON returned by EPP'...
                    '<br>has %d lost children!!</html>'], ...
                    this.lostChildren.size)));
            end
        end
        
        function ok=describeProgress(~) %this
            ok=false;
        end
        
    end
    
    methods(Static)
        function epp=Test(branch, useNotGate)
            if nargin<2
                useNotGate=false;
                if nargin<1
                    branch=SuhJsonSplitter.DEFAULT_BRANCH;
                end
            end
            epp=run_epp('eliver4D.csv', 'create_splitter', 'json',...
                'use_not_gate', useNotGate, 'cpp_branch', branch);
        end

        function epp=TestFlowJoRun(branch, useNotGate)
            if nargin<2
                useNotGate=false;
                if nargin<1
                    branch=SuhJsonSplitter.DEFAULT_BRANCH;
                end
            end
            epp=run_epp(File.Home('.flowjoAG', 'eppIn.csv'), ...
                'create_splitter', 'json',...
                'use_not_gate', useNotGate, 'cpp_branch', branch);
        end

        function [app, eppHome, cloudFldr, branchesFldr, ...
                thirdPartyFldr, versionFile]=Files(branch)
            if nargin<1
                branch=SuhJsonSplitter.DEFAULT_BRANCH;
            end
            eppHome=File.Home(SuhJsonSplitter.HOME);
            File.mkDir(eppHome);
            app=['EPPcli_' branch];
            if ispc
                app=[app '.exe'];
            end
            app=fullfile(eppHome, app);
            cloudFldr='EPP';
            branchesFldr=fullfile(eppHome, 'branches');
            thirdPartyFldr=fullfile(eppHome, 'thirdParty');
            versionFile=['version_' branch '.txt'];
        end
        
        function [ok, extraUrls]=Install(branch)
            if nargin<1
                branch=SuhJsonSplitter.DEFAULT_BRANCH;
            end
            ok=false;
            btnClang=Gui.NewBtn('Get Clang', ...
                @(h,e)getClang,...
                'Go to clang++ website');
            [choice, cancelled]=Gui.Ask([...
                '<html><center>Install EPP (branch=' branch ...
                ') <br>by downloading the:</center><hr></html>'], ...
                {'Executable', ...
                ['<html>Source and then building with <br>Clang ...'...
                '(<i>click install below</i></html>']},...
                'SuhJsonSplitter.Install', 'Confirm...', 1,...
                Gui.Panel(btnClang));
            extraUrls={};
            if cancelled || isempty(choice)
                return;
            end
            if choice==2
                ok=SuhJsonSplitter.Build(branch);
                if ~ok && askYesOrNo(...
                        'Download the executables?', ...
                        'Confirm', 'south')
                    choice=1;
                end
            end
            if choice==1
                [ok,~,extraUrls]=SuhJsonSplitter.Download;
            end
            
            function getClang
                if ismac
                    compiler='apple Xcode C++ compiler';
                elseif ispc
                    compiler='Microsoft Visual Studio C++ compiler';
                else
                    compiler='GNU C++ compiler';
                end
                msgBox(Html.WrapHr(['<b>Note</b>: clang++ is an optimized compiler<br>'...
                    'that works with your installation of<br>' ...
                    compiler]));

                web('https://clang.llvm.org/get_started.html', '-browser');
            end
        end
        
        
        function ok=Build(branch)
            if nargin<1
                branch=SuhJsonSplitter.DEFAULT_BRANCH;
            end
            ok=false;
            [app, eppHome, cloudFldr, branchesFldr, thirdPartyFldr]...
                =SuhJsonSplitter.Files(branch);
            branchesExists=exist(branchesFldr, 'dir');
            thirdPartyExists=exist(thirdPartyFldr, 'dir');
            
            if branchesExists && thirdPartyExists
                [refreshZip, cancelled]=askYesOrNo(...
                    'Refresh source from cloud?', ...
                    'Prior download exists...', 'center', true,...
                    '','SuhJsonSplitter.Refresh');
                if cancelled
                    return;
                end
            else
                refreshZip=true;
            end
            if refreshZip 
                zipFile=fullfile(eppHome, 'build.zip');
                delete(zipFile);
                [zipFile, exists, downloaded]=WebDownload.GetFile(...
                    'build.zip', eppHome, cloudFldr);
                if ~exists
                    return;
                end
                if ~ispc
                    fl=String.ToSystem(fullfile(eppHome, 'build'));
                    system(['chmod 777 ' fl]);
                end
                if thirdPartyExists
                    File.rmDir(fullfile(eppHome, 'thirdParty'));
                end
                if branchesExists
                    File.rmDir(fullfile(eppHome, 'branches'));
                end
                unzip(zipFile, eppHome)
            else
                zipFile=fullfile(eppHome, 'build.zip');
            end
            if exist(app, 'file')
                old=[app '.old'];
                File.moveFile(app, old, true);
            else
                old=[];
            end
            script=fullfile(eppHome, 'buildEppJson.cmd');
            
            cmds={['cd ' String.ToSystem(eppHome)]};
            if ispc
                cmds{2}=['call build ' branch];
            else
                cmds{2}=['./build ' branch];
            end
            File.Spawn(cmds, script, 'Building EPP with clang++', false, true);
            %script puts build in branch subfolder
            if ~exist(app, 'file')
                msgWarning('Build failed!!', 12, 'north');
                if askYesOrNo(['<html>Open windows to<br>'...
                        'edit build script?<hr></html>'])
                    File.OpenFolderWindow(eppHome, '', false)
                    if ispc
                        system('start cmd');
                        buildFile='build.cmd';
                    else
                        system('open -b com.apple.terminal');
                        buildFile='build';
                    end
                    msg(Html.Wrap(['The build folder is ' ...
                        Html.FileTree(eppHome) ...
                        '<br><br>Edit the file <b>' buildFile '</b> according' ...
                        '<br>to your C++ compiler''a needs<hr>']));
                end
                if ~isempty(old)
                    File.moveFile(old, app);
                end
            else
                ok=true;
                if ~isempty(old)
                    delete(old);
                end
                if exist(zipFile, 'file')
                    delete(zipFile);
                end
            end
        end
        
        function ok=GetUpdate(branch)
            pu=PopUp(...
                '<html>Checking Google Cloud<br>for EPP updates!</html>');
            if nargin<1
                branch=SuhJsonSplitter.DEFAULT_BRANCH;
            end
            ok=false;
            [app, eppHome, cloudFldr, ~, ~, vf]...
                =SuhJsonSplitter.Files(branch);
            f=fullfile(eppHome, vf);
            localVersion=str2double(strtrim(File.ReadTextFile(f)));
            if isnan(localVersion)
                WebDownload.GetFile(vf, eppHome, cloudFldr, true);
                localVersion=str2double(strtrim(File.ReadTextFile(f)));
            end
            if isnan(localVersion)
                pu.close;
                return;
            end
            url=WebDownload.ResolveUrl(vf, cloudFldr);
            host=WebDownload.UrlParts(url);
            if isempty(host)
                warning('Google Cloud can''t be reached for updates?');
                pu.close;
                return;
            end
            serverVersion=str2double(strtrim(WebDownload.ReadText(url)));
            if serverVersion>localVersion
                if askYesOrNo(struct(...
                        'javaWindow', pu.dlg,...
                        'msg', ['<html><center>The Google Cloud has EPP version '...
                        '<b>' num2str(serverVersion) '</b> available'...
                        '<br>... mean while you currently are running version <b>'...
                        num2str(localVersion) '</b><font color="red">'...
                        '&nbsp;&nbsp;:-( </font> ...<br><br><b>Download'...
                        ' the newer version</b>???<hr>'...
                        '</center></html>']), 'Good NEWS!!', 'south++')
                    ok=true;
                    delete(f)
                    delete(app);
                    SuhJsonSplitter.Download;
                end
            end
            pu.close;
        end
        
        function [ok, cancelled, extraUrls]=Download
            urls={};
            extraUrls={};
            [~, eppHome, cloudFldr]=SuhJsonSplitter.Files;
            if ispc
                urls{end+1}=WebDownload.ResolveUrl('EPPcli_main.exe', cloudFldr);
                urls{end+1}=WebDownload.ResolveUrl('EPPcli_bleeding.exe', cloudFldr);
                urls{end+1}=WebDownload.ResolveUrl('EPPcli_builds.exe', cloudFldr);
                urls{end+1}=WebDownload.ResolveUrl('libfftw3-3.dll', cloudFldr);
                urls{end+1}=WebDownload.ResolveUrl('libfftw3f-3.dll', cloudFldr);
                urls{end+1}=WebDownload.ResolveUrl('libfftw3l-3.dll', cloudFldr);   
                extraUrls{1}=urls{3};
                extraUrls{2}=urls{4};
                extraUrls{3}=urls{5};
            else
                urls{end+1}=WebDownload.ResolveUrl('EPPcli_main', cloudFldr);
                urls{end+1}=WebDownload.ResolveUrl('EPPcli_bleeding', cloudFldr);
                urls{end+1}=WebDownload.ResolveUrl('EPPcli_builds', cloudFldr);
            end
            urls{end+1}=WebDownload.ResolveUrl('version_main.txt', cloudFldr);
            urls{end+1}=WebDownload.ResolveUrl('version_bleeding.txt', cloudFldr);
            urls{end+1}=WebDownload.ResolveUrl('version_builds.txt', cloudFldr);
            [ok, cancelled]=WebDownload.Many(urls, eppHome);
            if ok
                if ismac
                    fl=String.ToSystem(fullfile(...
                        eppHome, 'EPPcli_builds'));
                    system(['chmod 777 ' fl]);
                    system(['xattr -r -d com.apple.quarantine ' fl]);
                    fl=String.ToSystem(fullfile(...
                        eppHome, 'EPPcli_bleeding'));
                    system(['chmod 777 ' fl]);
                    system(['xattr -r -d com.apple.quarantine ' fl]);
                    fl=String.ToSystem(fullfile(...
                        eppHome, 'EPPcli_main'));
                    system(['chmod 777 ' fl]);
                    system(['xattr -r -d com.apple.quarantine ' fl]);
                end
            end
        end
        
        function [app, reason]=GetTheApp(branch)
            if nargin<1
                branch=SuhJsonSplitter.DEFAULT_BRANCH;%get MOST tested
            end
            reason='';
            if ismac &&  MatBasics.OsVerCmp('10.12')<0
                reason=['<center>EPP requires <br>'...
                    'macOS Sierra (10.15) or later</center>'];
                app='';
            else
                app=['EPPcli_' branch];
                if ispc
                    app=[app '.exe'];
                end
                app=SuhJsonSplitter.Files(branch);
                if ~exist(app, 'file')
                    [~, extraUrls]=SuhJsonSplitter.Install(branch);
                else
                    extraUrls={};
                end
                if ~exist(app, 'file')
                    reason=['Need the app ',  Html.FileTree(app) ];
                    if ~isempty(extraUrls)
                        reason=[reason ' plus these :' ...
                            Html.ToList(extraUrls, 'ul') ];
                    end
                    app='';
                end
            end
            if isempty(app)
                msgError(Html.WrapHr(reason));
            else
                if Gui.UnderConstruction('EPP executable')
                    SuhJsonSplitter.GetUpdate(branch);
                end
            end
        end
        
        function json=RunFromFileSystem(columns, rows, csv, json, branch)
            if nargin<5
                branch=SuhJsonSplitter.DEFAULT_BRANCH;
            end
            app=SuhJsonSplitter.GetTheApp(branch);
            if isempty(app)
                json='';
                return;
            end
            fldr=fileparts(csv);
            outJson=[csv '_out.json'];
            if ismac
                macStdin=File.Home(SuhJsonSplitter.HOME, 'macStdin.txt');
                if ~exist(macStdin, 'file')
                    File.WriteTextFile(macStdin, {'1', '2'});
                end
                suffix=['<' macStdin];
            else
                suffix='';
            end
            fullCmd=[String.ToSystem(app) ' ' num2str(columns) ' ' ...
                num2str(rows) ' ' String.ToSystem(csv) ' ' ...
                String.ToSystem(json)  ' ' ...
                String.ToSystem(outJson) suffix];
            fprintf('This command is on the clipboard:\n  %s\n', fullCmd);
            clipboard('copy', fullCmd);
            script=fullfile(fldr, 'eppJson.cmd');
            [status, ~, isDone]=File.Spawn(fullCmd, script, ...
                ['<html><center>EPP on ' String.encodeInteger(rows) ...
                ' events X ' num2str(columns) ' measurements<br>' ...
                Html.WrapSmall(...
                '(<i>It is running in a separate command window</i>)')...
                '</center></html>'],...
                false, true);
            if isequal(File.Home('.flowjoAG'), fileparts(outJson)) ...
                    && exist(File.Home('.flowjoAG', 'gating.json'), 'file')
                if askYesOrNo('Use last FlowJo output?')
                    json=File.ReadTextFile(...
                        File.Home('.flowjoAG', 'gating.json'));
                else
                    json=File.ReadTextFile(outJson);
                end
            else
                json=File.ReadTextFile(outJson);
            end
            if status ~= 0 || ~isDone || isempty(json)
                msgError('EPP executable did not finish');
            else
                try
                    json=jsondecode(json);
                catch ex
                    json=[];
                    BasicMap.Global.reportProblem(ex);
                    throw(ex);
                end
            end
        end
        
        function json=EncodeJSON(isBalanced, W, sigma, ...
                kld1, kld2, kldExp, minRel, recursive)
            if nargin<8
                recursive=true;
                if nargin<7
                    minRel=.005;
                end
            end
            if recursive
                recursive='true';
            else
                recursive='false';
            end
            if isBalanced
                goal='best_balance';
            else
                goal='best_separation';
            end
            sigma=String.encodeRounded(sigma,1,true);
            W=String.encodeRounded(W,2);
            minRel=String.encodeRounded(minRel,3);
            kld1=String.encodeRounded(kld1,2);
            kldExp=String.encodeRounded(kldExp,2);
            kld2=String.encodeRounded(kld2,2);
            json=sprintf(SuhJsonSplitter.JSON_FMT, sigma, minRel, goal,...
                W, kld2, kld1, kldExp, recursive);
        end
        
        function json=Run(dataSet,args)
            if ~isempty(dataSet.file) && exist(dataSet.file, 'file') ...
                    && isempty(dataSet.labels)
                csv=dataSet.file;
            else
                csv=[tempname '.csv'];
                if ~isempty(dataSet.columnNames)
                    File.SaveMatrix(csv, dataSet.data, dataSet.columnNames);
                else
                    File.SaveMatrix(csv, dataSet.data, true);
                end
            end
            [p,f]=fileparts(csv);
            configJson=fullfile(p, [f '_config.json']);
            if isfield(args, 'minRelative')
                minR=args.minRelative;
            else
                minR=0.005;
            end
            jsonIn=SuhJsonSplitter.EncodeJSON(args.balanced, args.W,...
                args.sigma, args.KLD_normal_1D, args.KLD_exponential_1D,...
                args.KLD_normal_2D, minR);
            File.WriteTextFile(configJson,jsonIn);
            json=SuhJsonSplitter.RunFromFileSystem(...
                dataSet.C, dataSet.R, csv, configJson, args.cpp_branch);
        end
        
        
        function json=FindNode(json, key, lostChildren)
            try
                N=length(key);
                for i=2:N
                    idx=str2double(key(i));
                    if isstruct(json.children) %leaf
                        json=json.children;
                        if length(json)>1
                            json=json(idx);
                        else
                            if iscell(json.children)
                                json=json.children{idx};
                            else
                                json=json.children(idx);
                            end
                            if nargin>2
                                lostChild=key(1:idx);
                                if lostChildren.add(...
                                        java.lang.String(lostChild))
                                    disp(['Lost child:  ' lostChild]);
                                end
                            end
                        end
                    else
                        json=json.children{idx};
                    end
                end
            catch
                json=[];
            end
        end
        
    end
    
    methods(Access=protected)
        function [X, Y, polygonA, polygonB, leafCause]...
                =split(this, subset, key)
            if isempty(this.json)
                this.json=SuhJsonSplitter.Run(...
                    subset.dataSet, this.args);
                if isempty(this.json)
                    ex=MException(SuhJsonSplitter.MSGID_NO_JSON, ...
                        'EPP did not complete.');
                    throw(ex);
                end
            end
            leafCause='';
            try
                A=SuhJsonSplitter.FindNode(this.json, ...
                    [key '1'], this.lostChildren);
                if isempty(A)
                    X=0; Y=0; polygonA=[]; polygonB=[];
                    return;
                end
                if this.args.use_not_gate
                    polygonB='';
                else
                    B=SuhJsonSplitter.FindNode(this.json, ...
                        [key '2'], this.lostChildren);
                    if isempty(B)
                        X=0; Y=0; polygonA=[]; polygonB=[];
                        return;
                    end
                    polygonB=B.polygon;
                    assert(A.X==B.X && A.Y==B.Y);
                end
                X=A.X+1;
                Y=A.Y+1;
                polygonA=A.polygon;
            catch ex
                X=0; Y=0; polygonA=[]; polygonB=[];
                ex.getReport
                return;
            end
        end
    end
end