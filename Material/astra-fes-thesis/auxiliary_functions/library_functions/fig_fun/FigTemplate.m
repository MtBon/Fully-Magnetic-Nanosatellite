function FigTemplate(FIG)
% Apply common template to selected Figure
%

if nargin<1
    flist=findobj('Type', 'figure');
else
    flist.Number=FIG;
end
%% Parameters (may become input)
FontSize=16;
Font='Courier';
AxesLineWidth=0.75;
LineLineWidth=2;
Interpreter='tex';

%% Apply parameters to figure


for i=1:length(flist)
    FIG=flist(i).Number;
    % Axes
    h.ax=findall(FIG,'Type','Axes');
    if ~isempty(h.ax) % Check if you actually have a figure with axis... (why shouldn't you?)
        h.ax(1).FontSize=FontSize;
        h.ax(1).FontName=Font;
        h.ax(1).LineWidth=AxesLineWidth;
        h.ax(1).XLabel.Interpreter=Interpreter;
        h.ax(1).YLabel.Interpreter=Interpreter;
        h.ax(1).ZLabel.Interpreter=Interpreter;
        h.ax(1).TickLabelInterpreter=Interpreter;
        h.ax(1).Title.Interpreter=Interpreter;
        
        % Legend
        h.legend=findall(FIG,'Type','Legend');
        if ~isempty(h.legend)
            h.legend.FontSize=FontSize;
            h.legend.Interpreter=Interpreter;
            h.legend.Location='best';
        end
        
        % Lines
        h.lines=findall(FIG,'Type','Line');
        if ~isempty(h.lines)
            for j=1:length(h.lines)
                h.lines(j).LineWidth=LineLineWidth;
            end
        end
    end
end