function TimeAbs=RelTimeSec2AbsTime(TStart,TimeSec,formatOut)

if nargin==2
    formatOut = 'HH:MM:SS.FFF';
end
TimeAbsSec=TimeSec/24/3600+TStart;
TimeAbs=datestr(TimeAbsSec,formatOut);
