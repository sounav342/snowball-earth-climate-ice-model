function set_axis_limits(h,x,y)
%% Set axis limits tighter around min/max of y data than the
%% typical Matlab default.
%% If plotting more than one array using more than one plot call, use:
%% set_axis_limits(gca,x,[y1,y2])
%% If y1,y2 are column vectors, use
%% set_axis_limits(gca,x,[y1',y2'])
%% If wish to add min/max values in addition to the data (e.g. if
%% data are a constant, and a wider range is desired):
%% set_axis_limits(gca,x,[ymin,ymax,y1,y2])

%% set y limits:
y=reshape(y,1,[]); % returns a 1xN vector of y
V=axis; % returns x-axis and y-axis limits for current axes
a1=min(y);
a2=max(y);
V(3)=a1-(a2-a1)*0.15;
V(4)=a2+(a2-a1)*0.15;
if V(3)==V(4)
    V(3)=a1*0.9;
    V(4)=a2*1.1;
    if V(3)==0.0
        V(3)=-1;
        V(4)=1;
    end
end
%% set x limits:
V(1)=x(1);
V(2)=x(length(x));

axis(h,V);
set(h,'YTickMode','auto');
