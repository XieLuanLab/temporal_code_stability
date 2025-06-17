classdef RoiUtil < handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
    properties(Constant)
        V='9.6';
        POLYGON='impoly';
        RECTANGLE='imrect';
        ELLIPSE='imellipse';
        NEW_COLOR=[.1 .66 .9];
        EDIT_COLOR=[.71 .839 .06];
    end
    
    methods(Static)
        function ok=CanDoNew
            ok=~verLessThan('matlab',  RoiUtil.V);
        end
        
        function [roi, str]=NewForXy(ax, xy, tolerance, cbMoved)
            if nargin<4
                cbMoved=[];
                if nargin<3
                    tolerance=1;
                end
            end
            xy=xy(boundary(xy(:,1), xy(:,2), .1241), :);
            xy=double(edu.stanford.facs.swing.CurveSimplifier.ToArray(...
                xy, tolerance));
            [roi, str]=RoiUtil.NewPolygon(ax, xy, cbMoved);
        end
        
        function oldClr=SetColor(roi, clr)
            oldRoi=~RoiUtil.CanDoNew;
            if ~oldRoi
                oldClr=roi.Color;
                roi.Color=clr;
            else
                oldClr=roi.getColor;
                roi.setColor(clr);
            end
        end
        
        function oldLbl=SetLabel(roi, lbl)
            oldRoi=~RoiUtil.CanDoNew;
            try
                if ~oldRoi
                    oldLbl=roi.Label;
                    roi.Label=lbl;
                end
            catch
            end
        end
        
        
        function [roi, str]=NewPolygon(ax, xy, cbMoved)
            str=MatBasics.XyToString(xy);            
            oldRoi=~RoiUtil.CanDoNew;
            if ~oldRoi
                roi=drawpolygon(ax, ...
                    'Position', xy,...
                    'ContextMenu', [], 'SelectedColor', ...
                    RoiUtil.NEW_COLOR);
            else
                roi=impoly(ax, xy);
                roi.setColor(RoiUtil.NEW_COLOR);
            end
            if ~isempty(cbMoved)
                if oldRoi
                    roi.addNewPositionCallback(...
                        @(pos)RoiUtil.OldRoiMoved(cbMoved, roi));
                else
                    addlistener(roi, 'ROIMoved', ...
                        @(src,evt)RoiUtil.Moved(cbMoved, src));
                    
                end
                feval(cbMoved,roi);
            end
        end
        
        function roi=New(ax, type, cbMoved, cbMoving)
            oldRoi=~RoiUtil.CanDoNew;
            if strcmp(type, RoiUtil.RECTANGLE)
                if ~oldRoi
                    roi=drawrectangle(ax, 'ContextMenu', [],...
                        'Rotatable', true, 'SelectedColor', ...
                        RoiUtil.NEW_COLOR);
                else
                    roi=imrect(ax);
                end
            elseif strcmp(type, RoiUtil.ELLIPSE)
                if ~oldRoi
                    roi=drawellipse(ax, 'RotationAngle', 0, ...
                        'ContextMenu', [], 'SelectedColor', ...
                        RoiUtil.NEW_COLOR);
                else
                    roi=imellipse(ax);
                end
            else
                if ~oldRoi
                    roi=drawpolygon(ax, ...
                        'ContextMenu', [], 'SelectedColor', ...
                        RoiUtil.NEW_COLOR);
                else
                    roi=impoly(ax);
                end
            end
            if oldRoi
                roi.addNewPositionCallback(...
                    @(pos)RoiUtil.OldRoiMoved(cbMoved, roi));
            else
                if nargin>3
                    addlistener(roi, 'MovingROI', ...
                        @(src,evt)RoiUtil.Moving(cbMoving, src));
                end
                addlistener(roi, 'ROIMoved', ...
                    @(src,evt)RoiUtil.Moved(cbMoved, src));
                
            end
            feval(cbMoved,roi);
        end
        
        function OldRoiMoved(cb, roi)
            feval(cb, roi);
        end
        
        function Moving(cb, roi)
            feval(cb, roi);
        end

        function Moved(cb, roi)
            feval(cb, roi);
        end

        function Clicked(cb, roi, evt)
            feval(cb, roi, evt);
        end

        function ok=IsNewRoi(roi)
            ok=isa(roi, 'images.roi.Ellipse') ||...
                isa(roi, 'images.roi.Rectangle') ||...
                isa(roi, 'images.roi.Polygon');
        end
        
        function position=Position(roi)
            if ~RoiUtil.IsNewRoi(roi)
                position=roi.getPosition();
            else
                if isa(roi, 'images.roi.Ellipse')
                    c=get(roi, 'Center');
                    sa=get(roi, 'SemiAxes');
                    position=[c(1)-sa(1) c(2)-sa(2) sa(1)*2 sa(2)*2];
                    if roi.RotationAngle ~=0
                        position(end+1)=roi.RotationAngle;
                    end
                else
                    position=get(roi, 'Position');
                    if isa(roi, 'images.roi.Rectangle')
                        if roi.RotationAngle ~=0
                            position(end+1)=roi.RotationAngle;
                        end
                    end
                end
            end
        end
        
                
        function rows=GetRows(roi, data2D)
            if RoiUtil.IsNewRoi(roi)
                rows=roi.inROI(data2D(:,1), data2D(:,2));
                return;
            end
            pos=RoiUtil.Position(roi);
            typeT=class(roi);
            if strcmpi(typeT,RoiUtil.ELLIPSE)
                rows=RoiUtil.InEllipseUnrotated(data2D(:,1), data2D(:,2), pos);
            elseif strcmpi(typeT, RoiUtil.RECTANGLE)
                rows=RoiUtil.InRectUnrotated(data2D(:,1), data2D(:,2), pos);
            elseif strcmp(typeT, RoiUtil.POLYGON)
                rows=inpolygon(data2D(:,1), data2D(:,2),pos(:,1),pos(:,2));
            else
                rows=[];
            end
        end
        

        function inside=InEllipseUnrotated(X, Y, pos)
            xmin = pos(1);
            ymin = pos(2);
            width = pos(3);
            height = pos(4);
            a = width/2;
            b = height/2;
            center = [xmin+a, ymin + b];
            inside = (X - center(1)*ones(size(X))).^2./a^2 + ...
                (Y - center(2)*ones(size(Y))).^2./b^2 <= ones(size(X));
        end
        
        function inside=InRectUnrotated(X, Y, pos)
            xmin = pos(1);
            ymin = pos(2);
            width = pos(3);
            height = pos(4);
            a = width/2;
            b = height/2;
            center = [xmin+a, ymin + b];
            inside = abs(X - center(1)*ones(size(X))) <= a*ones(size(X))...
                & abs(Y - center(2)*ones(size(Y))) <= b*ones(size(Y));
        end
        
        function ok=IsHandleOk(roi)
            try
                ok=ishandle(roi.Parent);
                if isempty(ok) %odd .... does NOT look ok
                    ok=false;
                end
            catch
                ok=false;
            end
        end
        
    end
    
    properties
        position;
        label;
        roi;
        type;
        newRoi;%r2018b or later
        cbMoved;
    end
    
    methods
        function select(this, yes)
            if nargin<2
                yes=true;
            end
            if this.newRoi
                this.roi.Selected=yes;
                if yes
                    bringToFront(this.roi)
                end
            end
        end
        
        function setLabel(this, label)
            this.label=label;
            if this.newRoi
                this.roi.Label=label;
            end
        end
        
        function this=RoiUtil(type, pos)
            if isa(type, 'RoiUtil')
                prior=type;
                type=prior.type;
                if nargin<2
                    pos=prior.position;
                end
            else
                prior=[];
            end
            this.newRoi=RoiUtil.CanDoNew;
            this.type=type;
            this.position=pos;
            if this.newRoi
                if strcmp(type, RoiUtil.RECTANGLE)
                    this.roi=images.roi.Rectangle('Position', pos);
                elseif strcmp(type, RoiUtil.ELLIPSE)
                    semiAxes=[pos(3)/2 pos(4)/2];
                    center=[pos(1)+semiAxes(1) pos(2)+semiAxes(2)];
                    if length(pos)==5
                        this.roi=images.roi.Ellipse('RotationAngle', pos(5),...
                            'Center', center, 'SemiAxes', semiAxes);
                    else
                        this.roi=images.roi.Ellipse(...
                            'Center', center, 'SemiAxes', semiAxes);
                    end
                else
                    this.roi=images.roi.Polygon('Position', pos);
                end                
                if this.newRoi && ~isempty(prior)
                    this.roi.Label=prior.roi.Label;
                    this.roi.UserData=prior.roi.UserData;
                    this.roi.Color=prior.roi.Color;
                    this.roi.SelectedColor=prior.roi.SelectedColor;
                    this.setMovedCallback(prior.cbMoved);
                end
            end
        end
        
        function H=setParent(this, parent)
            if ~this.newRoi
                if strcmp(this.type, RoiUtil.RECTANGLE)
                    H=imrect(parent, this.position);
                elseif strcmp(this.type, RoiUtil.ELLIPSE)
                    H=imellipse(parent, this.position);
                else
                    H=impoly(parent, this.position);
                end
                this.roi=H;
            else
                try
                    H=this.roi;
                    H.Parent=parent;
                catch
                end
            end
        end
        
        function ok=isHandleOk(this)
            if ~this.newRoi
                ok=true;
            else
                try
                    ok=ishandle(this.roi.Parent);
                    if isempty(ok) %odd .... does NOT look ok
                        ok=false;
                    end
                catch
                    ok=false;
                end
            end
        end
        
        function pos=getPosition(this)
            pos=RoiUtil.Position(this.roi);
        end
        
        function setPosition(this, pos)
            this.position=pos;
            if this.newRoi
                if strcmp(this.type, RoiUtil.RECTANGLE)
                    this.roi.Position=pos;
                elseif strcmp(this.type, RoiUtil.ELLIPSE)
                    center=[pos(1)+pos(3)/2 pos(2)+pos(4)/2];
                    semiAxes=[pos(3)/2 pos(4)/2];
                    if length(pos)==5
                        this.roi.RotationAngle=pos(5);
                    end
                    this.roi.Center=center;
                    this.roi.SemiAxes=semiAxes;
                else
                    this.roi.Position=pos;
                end
            end
        end
        
        function rows=inROI(this, XY)
            if this.newRoi
                rows=this.roi.inROI(XY{1}, XY{2});
            else
                if strcmpi(this.type,RoiUtil.ELLIPSE)
                    rows=RoiUtil.InEllipseUnrotated(XY{1}, XY{2}, ...
                        this.position);
                elseif strcmpi(this.type, RoiUtil.RECTANGLE)
                    rows=RoiUtil.InRectUnrotated(XY{1}, XY{2}, ...
                        this.position);
                elseif strcmp(this.type, RoiUtil.POLYGON)
                    rows=inpolygon(XY{1}, XY{2},...
                        this.position(:,1), this.position(:,2));
                else
                    rows=[];
                end
            end
        end
        
        function setMovedCallback(this, cbMoved)
            if ~isempty(cbMoved)
                this.cbMoved=cbMoved;
                if ~this.newRoi
                    this.roi.addNewPositionCallback(...
                        @(pos)RoiUtil.OldRoiMoved(cbMoved, ...
                        this.roi));
                else
                    addlistener(this.roi, 'ROIMoved', ...
                        @(src,evt)RoiUtil.Moved(cbMoved, this.roi));
                end
            end
        end
        
        function setClickedCallback(this, cbClicked)
            if ~isempty(cbClicked)
                this.cbMoved=cbClicked;
                if ~this.newRoi
                    disp('No clicked callback on ROI BEFORE r2018b');
                else
                    addlistener(this.roi, 'ROIClicked', ...
                        @(src,evt)RoiUtil.Clicked(cbClicked, src, evt));
                end
            end
        end
        
    end
end