function status = solver_output(t,y,flag,varargin)
%PRINT
% Prints the current time to STDOUT and write time and vertices to log
% files. Plots area v. time using ODEPLOT.
%
% For use with an odesolver.
%
% See also: ODE23, ODEPLOT, ODEPRINT

% Status -- 0 is no problem / stop
tis = varargin{1};

% Write to log file - used for assembling model
if numel(t) < 2 % Skip the first line
    dlmwrite([tis.dir '/times.csv'],t,'delimiter',',','-append');
    dlmwrite([tis.dir '/vertices.csv'],y','delimiter',',','-append');
end

display( ['Time = ' num2str(t)] );

A = tis.getArea;

% --- ODEPLOT ---

persistent TARGET_FIGURE TARGET_AXIS TARGET_HGCLASS

status = 0;                             % Assume stop button wasn't pushed.
chunk = 128;                            % Memory is allocated in chunks.

drawnowDelay = 3;  % HGUsingMATLABClasses - postpone drawnow for performance reasons
doDrawnow = true;  % default

if nargin < 3 || isempty(flag) % odeplot(t,y) [v5 syntax] or odeplot(t,y,'')
    
    if (isempty(TARGET_FIGURE) || isempty(TARGET_AXIS))
        
        error(message('MATLAB:odeplot:NotCalledWithInit'));
        
    elseif (ishghandle(TARGET_FIGURE) && ishghandle(TARGET_AXIS))  % figure still open
        
        try
            ud = get(TARGET_FIGURE,'UserData');
            % Append t and y to ud.t and ud.y, allocating if necessary.
            nt = length(t);
            chunk = max(chunk,nt);
            [rows,cols] = size(ud.y);
            oldi = ud.i;
            newi = oldi + nt;
            if newi > rows
                ud.t = [ud.t; zeros(chunk,1)];
                ud.y = [ud.y; zeros(chunk,cols)];
            end
            ud.t(oldi+1:newi) = t;
            ud.y(oldi+1:newi,:) = A.'; % Plot area instead of position
            ud.i = newi;
            
            if TARGET_HGCLASS
                ploti = ud.ploti;
                doDrawnow = (ud.drawnowSteps > drawnowDelay);
                if doDrawnow
                    ud.ploti = newi;
                    ud.drawnowSteps = 1;
                else
                    ud.drawnowSteps = ud.drawnowSteps + 1;
                end
            end
            
            set(TARGET_FIGURE,'UserData',ud);
            
            if ud.stop == 1                       % Has stop button been pushed?
                status = 1;
            else
                % Rather than redraw all of the data every timestep, we will simply move
                % the line segments for the new data, not erasing.  But if the data has
                % moved out of the axis range, we redraw everything.
                ylim = get(TARGET_AXIS,'ylim');
                
                % Replot everything if out of axis range or if just initialized.
                if (oldi == 1) || (min(y(:)) < ylim(1)) || (ylim(2) < max(y(:)))
                    for j = 1:cols
                        set(ud.lines(j),'Xdata',ud.t(1:newi),'Ydata',ud.y(1:newi,j));
                    end
                else
                    % Plot only the new data.
                    if doDrawnow
                        if TARGET_HGCLASS  % start new segment
                            if ~ishold
                                hold on
                                plot(ud.t(ploti:newi),ud.y(ploti:newi,:),'-o');
                                hold off
                            else
                                plot(ud.t(ploti:newi),ud.y(ploti:newi,:),'-o');
                            end
                        else
                            for j = 1:cols
                                set(ud.line(j),'Xdata',ud.t(oldi:newi),'Ydata',ud.y(oldi:newi,j));
                            end
                        end
                    end
                end
            end
            
        catch ME
            error(message('MATLAB:odeplot:ErrorUpdatingWindow', ME.message));
        end
    end
    
else
    
    switch(flag)
        case 'init'                           % odeplot(tspan,y0,'init')
            TARGET_HGCLASS = feature('HGUsingMATLABClasses');
            
            ud = [];
            cols = length(A);
            ud.t = zeros(chunk,1);
            ud.y = zeros(chunk,cols);
            ud.i = 1;
            ud.t(1) = t(1);
            ud.y(1,:) = A.';
            
            if TARGET_HGCLASS
                EraseMode = {};
                ud.ploti = 1;
                ud.drawnowSteps = 1;
            else
                EraseMode = {'EraseMode','none'};
            end
            
            % Rather than redraw all data at every timestep, we will simply move
            % the last line segment along, not erasing it.
            f = figure(gcf);
            
            TARGET_FIGURE = f;
            TARGET_AXIS = gca;
            
            if ~ishold
                ud.lines = plot(ud.t(1),ud.y(1,:),'-o');
                hold on
                ud.line = plot(ud.t(1),ud.y(1,:),'-o',EraseMode{:});
                hold off
                set(TARGET_AXIS,'XLim',[min(t) max(t)]);
            else
                ud.lines = plot(ud.t(1),ud.y(1,:),'-o',EraseMode{:});
                ud.line = plot(ud.t(1),ud.y(1,:),'-o',EraseMode{:});
            end
            
            % The STOP button.
            h = findobj(f,'Tag','stop');
            if isempty(h)
                ud.stop = 0;
                pos = get(0,'DefaultUicontrolPosition');
                pos(1) = pos(1) - 15;
                pos(2) = pos(2) - 15;
                uicontrol( ...
                    'Style','pushbutton', ...
                    'String',getString(message('MATLAB:odeplot:ButtonStop')), ...
                    'Position',pos, ...
                    'Callback',@StopButtonCallback, ...
                    'Tag','stop');
            else
                set(h,'Visible','on');            % make sure it's visible
                if ishold
                    oud = get(f,'UserData');
                    ud.stop = oud.stop;             % don't change old ud.stop status
                else
                    ud.stop = 0;
                end
            end
            set(f,'UserData',ud);
            
        case 'done'                           % odeplot([],[],'done')
            
            f = TARGET_FIGURE;
            TARGET_FIGURE = [];
            ta = TARGET_AXIS;
            TARGET_AXIS = [];
            
            if ishghandle(f)
                ud = get(f,'UserData');
                ud.t = ud.t(1:ud.i);
                ud.y = ud.y(1:ud.i,:);
                set(f,'UserData',ud);
                cols = size(ud.y,2);
                for j = 1:cols
                    set(ud.lines(j),'Xdata',ud.t,'Ydata',ud.y(:,j));
                end
                if ~ishold
                    set(findobj(f,'Tag','stop'),'Visible','off');
                    
                    if ishghandle(ta)
                        set(ta,'XLimMode','auto');
                    end
                    
                    refresh;                          % redraw figure to remove marker frags
                end
            end
    end
end

if doDrawnow
    xlabel('Simulation time (sec)')
    ylabel('Cell area (\mum^2)')
    drawnow;
end

end

% --------------------------------------------------------------------------
% Sub-function
%

function StopButtonCallback(src,eventdata)
ud = get(gcbf,'UserData');
ud.stop = 1;
set(gcbf,'UserData',ud);
end  % StopButtonCallback

% --------------------------------------------------------------------------