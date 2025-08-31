function [start_time,end_time] = timegreg_2_timeCMIP(start_date,end_date,base_date,model_name)
% this code converts the start and end date for the considered period into
% CMIP time, i.e. number of days starting from a certain base date (usually
% 1850-01-01). Some models have leap years, some others have a fixed number
% of days (either 365 or 360)

% calculate the number of extra days over all the leap years of the time
% interval
if base_date.Year==start_date.Year
    n_leap_start = 0;
else
    n_leap_start = sum(leapyear(base_date.Year:1:start_date.Year));
end
n_leap_end = sum(leapyear(base_date.Year:1:(end_date.Year-1)));

switch model_name
    case 'ACCESS-ESM1-5'
        start_time  = 1 + days(start_date-base_date);
        end_time    = 1 + days(end_date-base_date);
    case 'INM-CM4-8'
        start_time  = 1 + days(start_date-base_date) - n_leap_start;
        end_time    = 1 + days(end_date-base_date) - n_leap_end;
    case 'MIROC-ES2L'
        start_time  = 1 + days(start_date-base_date);
        end_time    = 1 + days(end_date-base_date);
    case 'MRI-ESM2-0'
        start_time  = 1 + days(start_date-base_date);
        end_time    = 1 + days(end_date-base_date);
    case 'MPI-ESM1-2-LR'
        start_time  = 1 + days(start_date-base_date);
        end_time    = 1 + days(end_date-base_date);
    case 'MPI-ESM-1-0-HAM'
        start_time  = 1 + days(start_date-base_date);
        end_time    = 1 + days(end_date-base_date);
    case 'UKESM1-0-LL'
        start_time  = 1 + days(start_date-base_date) - n_leap_start - 5*(start_date.Year-base_date.Year);
        end_time    = 1 + days(end_date-base_date) - n_leap_end - 5*(end_date.Year-base_date.Year);        
    case 'GFDL-ESM4'
        start_time  = 1 + days(start_date-base_date) - n_leap_start;
        end_time    = 1 + days(end_date-base_date) - n_leap_end;
    case 'MIROC-ES2H'
        start_time  = 1 + days(start_date-base_date);
        end_time    = 1 + days(end_date-base_date);
    case 'proleptic_gregorian'
        start_time  = 1 + days(start_date-base_date);
        end_time    = 1 + days(end_date-base_date);        
    otherwise
        disp('Invalid model name: quitting')
        return
end

end