function res = az_equalize(ts)
% AZ_EQUALIZE  normalizes magnitude/phase response of sensor data based on
% characterized BUMP rev B transfer functions

res.fs = ts.fs;
res.data = filtfilt(b,a,ts.data);
