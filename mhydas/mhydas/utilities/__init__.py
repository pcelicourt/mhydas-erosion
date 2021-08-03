
from datetime import datetime as dt

def datenum(dates):
    return list(map(lambda _date: 366 + _date.toordinal() + (_date - dt.fromordinal(_date.toordinal())).total_seconds()/(24*60*60), dates))
