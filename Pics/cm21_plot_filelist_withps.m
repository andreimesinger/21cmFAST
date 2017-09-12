% 21cm_plot_filelist_withps

% open the filelists
[input_filename] = textread('input_box_filenames','%s');
[input_ps_filename] = textread('input_ps_filenames','%s');
color_bar = 1;
expand_factor = 1;
MIN = -52;
MAX = 52;
dim = 300;
for file_ct=1:length(input_filename)
    curr_filename = input_filename{file_ct};
    [params] = regexp(curr_filename, '\d+\.\d\d+', 'match');
    xH = params{1};
    xH_str = num2str(round(str2num(xH)*100)/100);
    z = params{2};
    params = regexp(curr_filename, 'Tb\S+\.\d\d', 'match');
    [aveTb] = regexp(params{1}, 'b', 'split');
    aveTb_str = num2str(round(str2num(aveTb{2})*10)/10);
%    outfilename = num2str(1000 - str2num(z) * 100);
    if (str2num(z) > 99.9999)
        outfilename = strcat(z, '_withps.tiff');
    else
         if (str2num(z) > 9.9999)
             outfilename = strcat('0', z, '_withps.tiff');
         else
             outfilename = strcat('00', z, '_withps.tiff');
         end
    end
    cm21_image_slice_withps(curr_filename, outfilename, dim, color_bar, expand_factor, z, xH_str, aveTb_str, MIN, MAX, input_ps_filename{file_ct});
end